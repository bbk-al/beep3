#include "prerequisites.h"

#include "explicit_worker.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "rhs_handler.decl.h"
#include "fmm_globals_nodegroup.decl.h"
#include "opencl_nodegroup.decl.h"
#include "iteration_handler.decl.h"
#include "fmm_worker.decl.h"

#include "explicit_worker.h"
#include "fmm_worker.h"
#include "rhs_handler.h"
#include "parallel_fmm_octree.h"
#include "fh_values_nodegroup.h"
#include "opencl_nodegroup.h"
#include "iteration_handler.h"
#include "mesh_library.decl.h"
#include "mesh_library.h"

#include "main.decl.h"
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

typedef ParallelFMMOctree<CharmNodePatch, CProxy_FMMWorker> ParallelTree;
//#define USE_KAHAN // we do not need to use compensated addition-- uncomment this if you really really want to

extern /* readonly */ CProxy_Main MainProxy;
extern /* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
extern /* readonly */ CProxy_MeshLibrary MeshLibraryProxy;
extern /* readonly */ beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
extern /* readonly */ CProxy_ExplicitWorker ExplicitWorkerProxy;
extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
extern /* readonly */ CProxy_IterationHandler IterationHandlerProxy;
extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

ExplicitWorker::ExplicitWorker() : possibly_migrating(false), done_precalcs(false)
{
    //std::cout << "Creating ExplicitWorker: " << thisIndex << std::endl;
}

void ExplicitWorker::precalc(double kappa, unsigned int total_num_patches)
{
    if (done_precalcs) { return; }
    
    // this should only get called once, so explicit_interactions list should be empty
    assert(explicit_integrations.size() == 0);

    // only do this function once
    done_precalcs = true;
    
    // get the corresponding node from the FMM tree
    const ParallelTree& tree = *(ParallelFMMOctreeProxy.ckLocalBranch());
    const ParallelTree::NodeT &node = tree.get_node(thisIndex);
    assert(node.isLeaf());
    assert(!node.empty());

    typedef ParallelTree::NodeT::ContentList ContentList;
    const ContentList& contents = node.get_contents();
    size_t num_patches = contents.size();
    
    // create local integrations for the self-geometric interactions (Cauchy terms)
    LintArray local_ints(new LocalIntegrations[num_patches]);
    explicit_integrations.push_back(LintArray_Size(local_ints, num_patches));

    for (size_t ii=0; ii < num_patches; ++ii)
    {
        const CharmNodePatch& np = *(contents[ii]);
     
        double dielectric_ratio = np.get_dielectric_ratio();
        float A=0,B=0,C=0,D=0;
        singular_BEM_kernels(kappa, np, A, B, C, D);
        
        local_ints[ii].set(np.get_idx(),
                        np.get_idx(),
                        A,
                        B + fGeometricCorrection(np.gc, dielectric_ratio),
                        C - hGeometricCorrection(np.gc, dielectric_ratio),
                        D);
    }

#ifdef CACHE_GPU_RESULTS
#ifdef OPENCL

    OpenCL_NodeGroupProxy[CkMyNode()].precalc_bem(thisIndex, kappa);
    
#else 

    // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
    boost::shared_ptr<ContentList> neighbourlist = tree.get_neighbourhood_contents(node);
    size_t num_neighbours = neighbourlist->size();
    size_t ctr=0;
    size_t chunksize=static_cast<size_t>(ceil(static_cast<double>(BEM_EXPLICIT_CHUNKSIZE) / node.size()));
    chunksize = (num_neighbours > chunksize) ? chunksize : num_neighbours;

    std::vector<boost::shared_ptr<QuadList> > qp_cache;
    
    while (ctr < num_neighbours)
    {
        boost::shared_ptr<std::vector<CharmNodePatch*> > list_chunk(new std::vector<CharmNodePatch*>());
        list_chunk->reserve(chunksize);

        for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
        {
	        qp_cache.push_back(((*neighbourlist)[ctr])->get_quadrature_points());
            list_chunk->push_back((*neighbourlist)[ctr++]);
        }

        size_t size = list_chunk->size() * node.size();
        LintArray results(new LocalIntegrations[size]);
        LintArray_Size arr_siz(results, size);
        explicit_integrations.push_back( arr_siz );

        size_t chunk_ctr=0;
        for (std::vector<CharmNodePatch*>::const_iterator src_it=node.get_contents().begin(), src_end=node.get_contents().end(); src_it != src_end; ++src_it)
        {
            const CharmNodePatch& src_patch = **src_it;
            qp_cache.push_back(src_patch.get_qualocation_points());
            for (std::vector<CharmNodePatch*>::const_iterator targ_it=list_chunk->begin(), targ_end=list_chunk->end(); targ_it != targ_end; ++targ_it)
            {
                const CharmNodePatch& targ_patch = **targ_it;
                assert(chunk_ctr < size);
                results[chunk_ctr++].init(kappa, src_patch, targ_patch);
            }
        }

    }
    
#endif
#endif

#ifdef USE_JUFFER_PEAKS
    // One-off calculation: surface integration of Juffer-style peaks
    // Use the f/h values embedded within the quadrature points
    juffer_peaks.reset(new double[total_num_patches*2]);
    memset(juffer_peaks.get(), 0, sizeof(double)*total_num_patches*2);
    for (ContentList::const_iterator np_it=contents.begin(), np_end=contents.end();
            np_it != np_end;
            ++np_it)
    {
        const CharmNodePatch& src_np = **np_it;
        unsigned int src_idx = src_np.get_idx();

        for (ContentList::const_iterator neigh_idx_it=neighbourhood->begin(), neigh_idx_end=neighbourhood->end();
                neigh_idx_it != neigh_idx_end;
                ++neigh_idx_it)
        {

            const CharmNodePatch& targ_np = **neigh_idx_it;
            if (targ_np.get_idx() != (**np_it).get_idx())
            {
                peak_integral(src_np, targ_np, juffer_peaks[src_idx], juffer_peaks[src_idx+total_num_patches], kappa);
            }
        }
    }
#endif // USE_JUFFER_PEAKS
}

void ExplicitWorker::peak_integral(const CharmNodePatch& src_patch, const CharmNodePatch& targ, double& fpeak, double& hpeak, double kappa)
{
    static const double sigm = 1e-20;
    float A=0,B=0,C=0,D=0;
    const float inv_epsilon = 1.0 / targ.get_dielectric_ratio();

    const float wt = targ.get_weighted_area() / 3.0;
    const Vector dx = targ.get_node() - src_patch.get_node();
    const Vector& n = targ.get_alt_normal();
    const Vector& n0 = src_patch.get_normal();

    // use inverse multiplication instead of division
    float r2 = (dx.length2() + sigm);
    float r = sqrt(r2);
    float ir2 = 1.0 / r2;
    float ir = 1.0 / r;
    float ir3 = ir2*ir;
    float ir4 = ir2*ir2;
    float ir5 = ir4*ir;

    float gr0 = ONE_OVER_4PI;
    float gr1 = gr0* ir;
    //float gr2 = gr0/r2;
    float gr3 = gr0*ir3;
    //float gr4 = gr0/r4;
    float gr5 = gr0*ir5;
    float ur0=exp(-kappa*r) * ONE_OVER_4PI;
    float ur1=ur0*ir;
    float ur2=ur0*ir2;
    float ur3=ur0*ir3;
    float ur4=ur0*ir4;
    float ur5=ur0*ir5;
    float ur3ur2=ur3+kappa*ur2;
    float pur4ur3=kappa*(3.0*ur4+kappa*ur3)+3.0*ur5;

    Vector gd = -dx*gr3;
    Vector ud = -dx*(ur3 + kappa*ur2);

    A += targ.h*(gr1 - ur1);
    B += targ.f*n.dot(inv_epsilon*gd - ud);
    C += -targ.h*n0.dot(gd - ud*inv_epsilon);
    {
        Vector dx3jk = dx * dx.x;
        Vector gdd_k = 3.0 * dx3jk * gr5;
        Vector udd_k = dx3jk * pur4ur3;
        gdd_k.x -= gr3;
        udd_k.x -= ur3ur2;
        D += - targ.f * n0.x * (n.dot(gdd_k-udd_k));
    }
    {
        Vector dx3jk = dx * dx.y;
        Vector gdd_k = 3.0 * dx3jk * gr5;
        Vector udd_k = dx3jk * pur4ur3;
        gdd_k.y -= gr3;
        udd_k.y -= ur3ur2;
        D += - targ.f * n0.y * (n.dot(gdd_k-udd_k));
    }
    {
        Vector dx3jk = dx * dx.z;

        Vector gdd_k = 3.0 * dx3jk * gr5;
        Vector udd_k = dx3jk * pur4ur3;
        gdd_k.z -= gr3;
        udd_k.z -= ur3ur2;
        D += - targ.f * n0.z * (n.dot(gdd_k-udd_k));
    }

    D *= inv_epsilon;

    fpeak += wt*(B - A);
    hpeak += wt*(D - C);

    return;
}

void ExplicitWorker::evaluate(double kappa)
{
    //std::cout << "Evaluating explicit neighbours for node " << thisIndex << std::endl;
    // get fvals & hvals from the node group
    FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());
    const double *fvals = fh_vals_group.fvals();
    const double *hvals = fh_vals_group.hvals();
    unsigned int total_num_patches = fh_vals_group.get_num_patches();
    boost::scoped_array<double> results(new double[total_num_patches*2]);
    memset(results.get(), 0, sizeof(double)*total_num_patches*2);
#ifdef USE_KAHAN
    boost::scoped_array<double> kahan(new double[total_num_patches*2]);
    memset(kahan.get(), 0, sizeof(double)*total_num_patches*2);
#endif

    if (done_precalcs == false)
    {
        precalc(kappa, total_num_patches);
#ifdef OPENCL
        OpenCL_NodeGroupProxy.ckLocalBranch()->wait_for_completion();
#endif
    }

#ifdef USE_JUFFER_PEAKS
    for (unsigned int ii=0; ii < total_num_patches*2; ++ii)
    {
        results[ii] += juffer_peaks[ii];
    }
#endif

    // loop over the precalculated local integrations and multiply f/hvals by the matrix elements
    for (std::vector<LintArray_Size>::const_iterator it=explicit_integrations.begin(), end=explicit_integrations.end(); it != end; ++it)
    {
        const LintArray& array = it->first;
        const size_t& sz = it->second;
        //std::cout << "Evaluating: " << sz << " local ints, from mem: " << array.get() << std::endl;
        for (size_t ii=0; ii < sz; ++ii)
        {
            //std::cout << ii << " of " << sz << ": " << array[ii] << "\n";
#ifdef USE_KAHAN
        array[ii].evaluate_local_contributions(fvals, hvals, results.get(), &(results[total_num_patches]), kahan.get(), &(kahan[total_num_patches]));
#else
        array[ii].evaluate_local_contributions(fvals, hvals, results.get(), &(results[total_num_patches]));
#endif
        }
    }

    // add them into the node group holder
    fh_vals_group.add_results(results.get());

#ifndef CACHE_GPU_RESULTS
#ifdef OPENCL
    // In this clause we are not stashing the precalc'ed values, but 
    // are generating them on the fly.  So call the Nodegroup Proxy and
    // tell it to work.
    OpenCL_NodeGroupProxy[CkMyNode()].run_bem(thisIndex, kappa);
    return;
#else

    // get the corresponding node from the FMM tree
    const ParallelTree tree = *(ParallelFMMOctreeProxy.ckLocalBranch());
    const ParallelTree::NodeT &node = tree.get_node(thisIndex);
    assert(node.isLeaf());

    typedef ParallelTree::NodeT::ContentList ContentList;
    const ContentList& contents = node.get_contents();
    size_t num_patches = contents.size();
    boost::shared_ptr<ContentList> neighbourhood = tree.get_neighbourhood_contents(node);
    
    // cacheing of quadrature points (compromise between massive memory usage and 
    // inefficiency of recalculating quad points)
    std::vector<boost::shared_ptr<QuadList> > qp_cache;
    for (ContentList::const_iterator targ_it=neighbourhood->begin(), 
				      targ_end=neighbourhood->end(); 
				      targ_it != targ_end; ++targ_it)
    { 
	(**targ_it).obtain_shared_quad_pointers(qp_cache);
    }
    

    LocalIntegrations bem_kernels;
    for (ContentList::const_iterator src_it=node.get_contents().begin(), 
                                     src_end=node.get_contents().end(); 
                                     src_it != src_end; ++src_it)
    {
        const CharmNodePatch& src_patch = **src_it;
        for (ContentList::const_iterator targ_it=neighbourhood->begin(), 
                                         targ_end=neighbourhood->end(); 
                                         targ_it != targ_end; ++targ_it)
        { 
            const CharmNodePatch& targ_patch = **targ_it;
            bem_kernels.init(kappa, src_patch, targ_patch);
            bem_kernels.evaluate_local_contributions(f_lhs.get(), h_lhs.get(), f_results, h_results);
        }
    }
    
#endif
#endif

    IterationHandlerProxy.done_evals();

}

void ExplicitWorker::pup(PUP::er &p) {
    CBase_ExplicitWorker::pup(p);
    //std::cout << "AARGH! ExplicitWorker " << thisIndex << " is pup'ing! - " << usesAtSync << std::endl;

    p | possibly_migrating;
    p | done_precalcs;
   
    // packing/unpacking the local integrations is a bit tricky because it's in 
    // multiple chunks of dynamically allocated memory wrapped in boost 
    // shared array pointers.
    size_t num_arrays=0;
    if (p.isPacking() || p.isSizing())
    {
        num_arrays = explicit_integrations.size();
    }
    p | num_arrays;
    
    for (size_t ii=0; ii < num_arrays; ++ii)
    {
        size_t chunk_size;
        if (p.isPacking() || p.isSizing())
        {
            chunk_size = explicit_integrations[ii].second;
        }
        p | chunk_size;
        
        LintArray lints;
        if (p.isPacking() || p.isSizing())
        {
            lints = explicit_integrations[ii].first;
        }
        if (p.isUnpacking())
        {
            lints.reset(new LocalIntegrations[chunk_size]);
            explicit_integrations.push_back(LintArray_Size(lints, chunk_size));
        }
        
        for (size_t jj=0; jj < chunk_size; ++jj)
        {
            p | lints[jj];
        }
    }

}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
ExplicitWorker::ExplicitWorker(CkMigrateMessage *msg) {
    //std::cout << "AARGH!  ExplicitWorker " << thisIndex << " is migrating!" << std::endl;
}

#include "explicit_worker.def.h"
