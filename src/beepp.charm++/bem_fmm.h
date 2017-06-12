#ifndef __BEM_FMM_H
#define __BEM_FMM_H

#include "prerequisites.h"
#include <string>
#include <vector>
#include "../fmm/fmm_bem_hybrid.h"

#include "fh_values_nodegroup.decl.h"
#include "fh_values_nodegroup.h"
#include "opencl_nodegroup.decl.h"
#include "opencl_nodegroup.h"

extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

template<int NTERMS, int NLAMBS, int NWAVES, typename PatchType> 
         class BEM_FMM : public CBase_BEM_FMM<NTERMS, NLAMBS, NWAVES, PatchType>, 
                                   public fmm::FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,PatchType>
{

public:

    typedef CBase_BEM_FMM<NTERMS, NLAMBS, NWAVES, PatchType> super_charm;
    typedef fmm::FMM_BEM_Octree<NTERMS,NLAMBS,NWAVES,PatchType> super;
    typedef typename super::NodeT NodeT;

    /// Constructors ///
    BEM_FMM(unsigned int max_items_per_node,
                        Vector universe_centre,
                        double universe_edge_length, 
                        unsigned int expected_total=0,
                        CkCallback cb_completed=CkCallback(CkCallback::ignore)) : super(max_items_per_node, universe_centre, universe_edge_length), expect(expected_total), stashed_cb(cb_completed)
    {
        // for BEM-FMM we optimize the tree to have average number of neighbours,
        // so we relax the subdivision policy on octree nodes.
        // The usual octree divides an octree node as soon as the contents
        // reach the max_items_per_node limit.
        // NB: since we unset_strict_divisions we *must*
        // call optimize_neighbourhoods once we're done inserting
        // items.
        super::unset_strict_divisions();
        
        done_precalc = false;
        if (expected_total == 0) {
            super_charm::contribute(cb_completed);
        }
    }
                        
    // should never migrate- it's a group object.
    BEM_FMM(CkMigrateMessage *msg) : super(0, Vector(0,0,0), 0.0) { throw std::exception(); }
    ~BEM_FMM() { }

    /// Entry methods
    void insert(std::vector<PatchType> insertions) {
        CkGroupID gid = super_charm::thisgroup;

        for (typename std::vector<PatchType>::const_iterator it=insertions.begin(), end=insertions.end();
                it != end;
                ++it)
        {
            super::add(*it);
        }
        
        // check filledness
        if (super::size() == expect) {
            finalize(stashed_cb);
        }
    }
    
    void finalize(CkCallback cb)
    {
        // delete nodes which this compute group is not interested in
        decision_level = level_with_more_nodes_than_procs();
        super::linearize(decision_level);
        
        typename super::NodeList& nlist = super::get_node_list(decision_level);
        size_t num_nodes = nlist.size();
        int num_nodes_per_pe = static_cast<int>(floor(static_cast<double>(num_nodes) / CkNumPes()));
        int pe_ctr=0, nctr=0;
        for (typename super::NodeList::iterator it=nlist.begin(), end=nlist.end(); it != end; ++it)
        {
            if (pe_ctr != CkMyPe())
            {
                super::erase_children(*(it->second));
                it->second->unset_is_leaf();
                it->second->set_shadow();
            }
            else
            {
                my_responsibilities.push_back(it->second->get_idx());
            }
            
            if (++nctr >= num_nodes_per_pe) {
                ++pe_ctr;
                if (pe_ctr >= CkNumPes()) { pe_ctr = 0; }
                nctr=0;
            }
        }
        
        // finalize the tree
        super::build_neighbourhoods();
        super::optimize_neighbourhoods();
        super::remove_empty_nodes();
        
        super_charm::contribute(cb);
    }
    
    void solve(double kappa, double kappa0, CkCallback cb)
    {
        FH_Values_NodeGroupProxy.ckLocalBranch();
        FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());
        unsigned int num_patches = fh_vals_group.get_num_patches();
        const double *f_lhs = fh_vals_group.fvals();
        const double *h_lhs = fh_vals_group.hvals();
        boost::scoped_array<double> results(new double[num_patches*2]);
        memset(results.get(), 0, sizeof(double)*num_patches*2);
        
        super::solve(kappa, kappa0, f_lhs, h_lhs, &(results[0]), &(results[num_patches]));
        
        // add the results into the group object
        fh_vals_group.add_results(results.get());
        
        // signal done fmm evals
        cb.send();
    }
    
    void precalc(double kappa);
    void self_patch_integrals(double kappa);
    void near_field_integration(double kappa, CkCallback cb);
       
private:

    unsigned int expect;
    CkCallback stashed_cb;
    
    double kappa;
    unsigned short decision_level;
    std::vector<OctreeIndexer> my_responsibilities;

    std::vector<LintArray_Size> local_integrations;
    bool done_precalc;
    
    unsigned short level_with_more_nodes_than_procs() const
    {
        for (unsigned short level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
        {
            if ((super::get_num_nodes_on_level(level)) >= CkNumPes()) 
            {
                return level;
            }
            
        }
        return super::get_bottom_level();
    }
    
};

template <int NTERMS, int NLAMBS, int NWAVES, typename PatchType>
void BEM_FMM<NTERMS, NLAMBS, NWAVES, PatchType>::self_patch_integrals(double kappa)
{

    long start_clock = myclock();
    typedef typename super::NodeT NodeT;
    typedef typename NodeT::ContentList ContentList;
    unsigned int total_singular_integrals=0;
    
    for (std::vector<OctreeIndexer>::iterator it=my_responsibilities.begin(), end=my_responsibilities.end();
          it != end; ++it)
    {
        // get node of responsibility
        const NodeT* node_ptr;
        try {
            node_ptr = &(super::get_node(*it));
        }
        catch (BadIndexer) {
            std::cerr << "A BEM_FMM compute group caught a bad node in self-patch integrals." << std::endl;
            CkExit(); // bail
        }
        const NodeT& node = *node_ptr;
        assert(node.isShadow() == false);
        assert(node.empty() == false);
        const ContentList& contents = node.get_contents();
        unsigned int num_patches = contents.size();
        total_singular_integrals += num_patches;

        // diagonal values
        LintArray diagonal_values(new LocalIntegrations[num_patches]);
        local_integrations.push_back(LintArray_Size(diagonal_values, num_patches));

        #pragma omp parallel for
        for (int ii=0; ii < num_patches; ++ii)
        {
            const BasicNodePatch& np = *(contents[ii]);
            double dielectric_ratio = np.get_dielectric_ratio();

            diagonal_values[ii].set(np.get_idx(),
                                    np.get_idx(),
                                    0,
                                    +fGeometricCorrection(np.gc, dielectric_ratio),
                                    -hGeometricCorrection(np.gc, dielectric_ratio),
                                    0);
                            
        }
        
#ifdef OPENCL

        // need this for OpenCL
        OpenCL_Handler& global_ocl_handler = OpenCL_NodeGroupProxy.ckLocalBranch()->get_ocl_handler();

        // create local integrations for the self-geometric interactions
        LintArray gpu_local_ints(new LocalIntegrations[num_patches]);
        local_integrations.push_back(LintArray_Size(gpu_local_ints, num_patches));
        
        unsigned int chunksize = 512;
        unsigned int ii=0;
        while (ii < num_patches)
        {
    
            PatchPtrList patch_ptrs(new PPList);
            
            unsigned chunk_ctr=0;
            for ( ; ii+chunk_ctr < num_patches && chunk_ctr < chunksize; ++chunk_ctr)
            {
                const BasicNodePatch& np = *(contents[ii]);
                patch_ptrs->push_back(&np);
            }
            
            // singular part
            SingularBEM* ocl_singular_bem = new SingularBEM(patch_ptrs, kappa, &(gpu_local_ints[ii]));
            global_ocl_handler.add_work_to_queue(ocl_singular_bem);
            
            ii += chunk_ctr;
        }
        
#else
        // create local integrations for the self-geometric interactions
        LintArray local_ints(new LocalIntegrations[num_patches]);
        local_integrations.push_back(LintArray_Size(local_ints, num_patches));
        
        for (int ii=0; ii < num_patches; ++ii)
        {
            const BasicNodePatch& np = *(contents[ii]);
            double dielectric_ratio = np.get_dielectric_ratio();
            
            float A=0,B=0,C=0,D=0;
            singular_BEM_kernels(kappa, np, A, B, C, D);
            local_ints[ii].set(np.get_idx(), np.get_idx(), A, B, C, D);
        }
#endif
    }

    std::cout << total_singular_integrals << " singular integrals took: " << (myclock() - start_clock) / 1000. << " ms" << std::endl;
       
    return;
}

template <int NTERMS, int NLAMBS, int NWAVES, typename PatchType>
void BEM_FMM<NTERMS, NLAMBS, NWAVES, PatchType>::precalc(double kappa)
{
    // only do precalcs once
    if (done_precalc) { return; }
    done_precalc = true;
    
    long start_clock = myclock();
    
    // clear any previously stored local integration results
    local_integrations.clear();

    // singular self-patch integrations
    self_patch_integrals(kappa);

    size_t total_interactions=0;
    
#ifdef OPENCL
    // need this for OpenCL
    OpenCL_Handler& global_ocl_handler = OpenCL_NodeGroupProxy.ckLocalBranch()->get_ocl_handler();
#endif        
    
#ifdef CACHE_GPU_RESULTS
    
    typedef typename super::NodeList NodeList;
    typedef typename super::NodeT NodeT;
    typedef typename NodeT::ContentList ContentList;
    
    // traverse all leaf nodes
    for (int level=decision_level; level <= super::get_bottom_level(); ++level)
    {
        // quadrature_point cache
        std::vector<boost::shared_ptr<QuadList> > qp_cache;

        for (typename NodeList::const_iterator node_it=super::get_node_list(level).begin(), 
                                              node_end=super::get_node_list(level).end(); 
               node_it != node_end; 
               ++node_it)
        {
            const NodeT& node = *(node_it->second);
            if (node.isShadow()==true || node.isLeaf()==false || node.empty()==true || node.isDeleted()==true) { continue; }

            // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
            boost::shared_ptr<std::vector<CharmNodePatch*> > neighbourlist = super::get_neighbourhood_contents(node);
            size_t num_neighbours = neighbourlist->size();
            
            size_t chunksize = static_cast<size_t>(ceil(BEM_EXPLICIT_CHUNKSIZE/static_cast<double>(node.size())));
            chunksize = (num_neighbours > chunksize) ? chunksize : num_neighbours;
            
            size_t ctr=0;
            while (ctr < num_neighbours)
            {
                boost::shared_ptr<std::vector<const BasicNodePatch*> > list_chunk(new std::vector<const BasicNodePatch*>());
                list_chunk->reserve(chunksize);

                for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
                {
                    ((*neighbourlist)[ctr])->obtain_shared_quad_ptrs(qp_cache);
                    list_chunk->push_back((*neighbourlist)[ctr++]);
                }

                size_t size = list_chunk->size() * node.size();
                LintArray near_field_integrals(new LocalIntegrations[size]);
                LintArray_Size arr_siz(near_field_integrals, size);
                local_integrations.push_back( arr_siz );
                
                // running total
                total_interactions += size;

#ifdef OPENCL
                //std::cout << "Size of neighbourlist: " << neighbourlist->size() << std::endl;
                typedef std::vector<const BasicNodePatch*> PatchList;
                boost::shared_ptr< PatchList > xlated_list(new PatchList);
                for (typename ContentList::const_iterator xx_it=node.get_contents().begin(), xx_end=node.get_contents().end(); 
                      xx_it != xx_end; ++xx_it)
                {
                    xlated_list->push_back(static_cast<const BasicNodePatch*>(*xx_it));
                }
                BEM_Resources* res_ptr = new BEM_Resources(xlated_list, list_chunk, kappa, near_field_integrals.get());
                global_ocl_handler.add_work_to_queue(res_ptr);
#else
                #pragma omp parallel
                for (int src_ctr=0; src_ctr < node.get_contents().size(); ++src_ctr)
                {
                    const BasicNodePatch& src_patch = *(node.get_contents()[src_ctr]);
                    
                    #pragma omp for
                    for (int targ_ctr=0; targ_ctr < list_chunk->size(); ++targ_ctr)
                    {
                        const BasicNodePatch& targ_patch = *((*list_chunk)[targ_ctr]); // that's a bit ugly
                        
                        size_t idx = (src_ctr * list_chunk->size()) + targ_ctr;
                        near_field_integrals[idx].init(kappa, src_patch, targ_patch);
                    }
                }
#endif

            }
        }

#ifdef OPENCL
        // ensure that QuadPoints stay in scope for OpenCL
        OpenCL_WorkBlob* ocl_cache_ptr = new QuadPointCache(qp_cache);
        global_ocl_handler.add_work_to_queue(ocl_cache_ptr);
#endif
 
    }

#endif

#ifdef OPENCL
    // need this for OpenCL
    global_ocl_handler.wait_until_idle();
#endif

    //std::cout << "Done " << local_integrations.size() << " groups of pair-wise local neighbour interactions" << std::endl;
    std::cout << "Precalc integrals (incl. singulars) took: " << (myclock() - start_clock) / 1000. << " ms (" << total_interactions << " interactions)" << std::endl;

    return;

}

template <int NTERMS, int NLAMBS, int NWAVES, typename PatchType>
void BEM_FMM<NTERMS, NLAMBS, NWAVES, PatchType>::near_field_integration(double kappa, CkCallback cb)
{
    precalc(kappa);

    // where to put results
    FH_Values_NodeGroupProxy.ckLocalBranch();
    FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());
    unsigned int num_patches = fh_vals_group.get_num_patches();
    const double *f_lhs = fh_vals_group.fvals();
    const double *h_lhs = fh_vals_group.hvals();
    boost::scoped_array<double> results(new double[num_patches*2]);
    memset(results.get(), 0, sizeof(double)*num_patches*2);

    // loop over the precalculated local integrations and multiply f/hvals by the matrix elements
    // (There are always at least the singular/cauchy principle values precalculated, if not entire
    // massive chunks of the BEM matrix... (unless using GPU in which case we generally do those
    // integrals on the fly as that's faster than shunting GB of data off the GPU.)
    for (std::vector<LintArray_Size>::const_iterator it=local_integrations.begin(), 
                                                    end=local_integrations.end(); 
          it != end; ++it)
    {
        const LintArray& array = it->first;
        const size_t& sz = it->second;
        //std::cout << "Evaluating: " << sz << " local ints, from mem: " << array.get() << std::endl;
        for (size_t ii=0; ii < sz; ++ii)
        {
            //std::cout << ii << " of " << sz << ": " << array[ii] << "\n";
            array[ii].evaluate_local_contributions(f_lhs, h_lhs, &(results[0]), &(results[num_patches]));
        }
    }

#ifndef CACHE_GPU_RESULTS

#ifdef OPENCL
    OpenCL_Handler& global_ocl_handler = OpenCL_NodeGroupProxy.ckLocalBranch()->get_ocl_handler();
#endif
    
    // temp holder for integrations
    LocalIntegrations bem_kernels;

    typedef typename super::NodeList NodeList;
    typedef typename super::NodeT NodeT;
    typedef typename NodeT::ContentList ContentList;
    
    // traverse all leaf nodes
    for (int level=decision_level; level <= super::get_bottom_level(); ++level)
    {
        // quadrature_point cache
        std::vector<boost::shared_ptr<QuadList> > qp_cache;

        for (typename NodeList::iterator node_it=super::get_node_list(level).begin(), 
                                        node_end=super::get_node_list(level).end(); 
               node_it != node_end; 
               ++node_it)
        {
            NodeT& node = *(node_it->second);
            if (node.isShadow()==true || node.isLeaf()==false || node.empty()==true || node.isDeleted()==true) { continue; }

            // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
            boost::shared_ptr<std::vector<CharmNodePatch*> > neighbourlist = super::get_neighbourhood_contents(node);
            size_t num_neighbours = neighbourlist->size();
            
            size_t chunksize = static_cast<size_t>(ceil(BEM_EXPLICIT_CHUNKSIZE/static_cast<double>(node.size())));
            chunksize = (num_neighbours > chunksize) ? chunksize : num_neighbours;
            
            size_t ctr=0;
            while (ctr < num_neighbours)
            {
                boost::shared_ptr<std::vector<const BasicNodePatch*> > list_chunk(new std::vector<const BasicNodePatch*>());
                list_chunk->reserve(chunksize);

                for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
                {
                    ((*neighbourlist)[ctr])->obtain_shared_quad_ptrs(qp_cache);
                    list_chunk->push_back((*neighbourlist)[ctr++]);
                }
            
#ifdef OPENCL
                //std::cout << "Size of neighbourlist: " << neighbourlist->size() << std::endl;
                PatchPtrList xlated_list(new PPList);
                for (typename ContentList::iterator xx_it=node.get_contents().begin(), xx_end=node.get_contents().end(); 
                      xx_it != xx_end; ++xx_it)
                {
                    xlated_list->push_back(static_cast<const BasicNodePatch*>(*xx_it));
                }
                BEM_OnDemand_Resources* res_ptr = new BEM_OnDemand_Resources(xlated_list, 
                                                                             list_chunk, 
                                                                             kappa, 
                                                                             f_lhs,
                                                                             h_lhs,
                                                                             &(results[0]),
                                                                             &(results[num_patches]));
                global_ocl_handler.add_work_to_queue(res_ptr);
#else
                size_t another_ctr=0;
                for (std::vector<CharmNodePatch*>::const_iterator src_it=node.get_contents().begin(), src_end=node.get_contents().end(); src_it != src_end; ++src_it)
                {
                    const BasicNodePatch& src_patch = **src_it;
                    for (std::vector<const BasicNodePatch*>::const_iterator targ_it=list_chunk->begin(), targ_end=list_chunk->end(); targ_it != targ_end; ++targ_it)
                    {
                        const BasicNodePatch& targ_patch = **targ_it;
                        bem_kernels.init(kappa, src_patch, targ_patch);
                        bem_kernels.evaluate_local_contributions(f_lhs, h_lhs, &(results[0]), &(results[num_patches]));

                    }
                }
#endif
            }
        }
        
#ifdef OPENCL
        // ensure that QuadPoints stay in scope for OpenCL
        OpenCL_WorkBlob* ocl_cache_ptr = new QuadPointCache(qp_cache);
        global_ocl_handler.add_work_to_queue(ocl_cache_ptr);
#endif

    }

#ifdef OPENCL
    global_ocl_handler.wait_until_idle();
#endif
    
#endif

    // add the results into the group object
    fh_vals_group.add_results(results.get());
    
    // signal completion
    cb.send();

}

#endif // __BEM_FMM_H
