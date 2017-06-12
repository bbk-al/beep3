/*
* fmm_octree.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef FMM_OCTREE_H_
#define FMM_OCTREE_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_ptr.hpp>

#include "octree.h"
#include "fmm_globals.h"
#include "fmm.h"
#include "../common/multipole_holder.h"
#include "interaction_list.h"
#include "eval_pt.h"
#include "../common/charge.h"
#include "fmm_time_info.h"

#ifdef OPENCL
class OpenCL_Handler;
#include "opencl_fmm.h"
#endif

namespace fmm
{

template <typename ContentType>
class FMM_NodeT : public Node<ContentType>
{
    typedef Node<ContentType> super;
public:
    FMM_NodeT(const Vector& _centre, double _edge_length, const OctreeIndexer& root_idxer) : super(_centre, _edge_length, root_idxer), pw_idx(-1), mp_idx(-1), lc_idx(-1) {}
    FMM_NodeT(FMM_NodeT& other, OctreeIndexer idxer, unsigned short child_id) :
        super(other, idxer, child_id),
        pw_idx(-1),
        mp_idx(-1),
        lc_idx(-1) {}
    virtual ~FMM_NodeT() {}

    inline long get_pw_idx() const { assert(pw_idx != -1); return pw_idx;}
    inline void set_pw_idx(long val) {pw_idx = val;}
    inline long get_mp_idx() const {assert(mp_idx != -1); return mp_idx;}
    inline void set_mp_idx(long val) {mp_idx = val;}
    inline long get_lc_idx() const {assert(lc_idx != -1); return lc_idx;}
    inline void set_lc_idx(long val) {lc_idx = val;}

protected:

    long pw_idx;
    long mp_idx;
    long lc_idx;
};

template<typename NodeT, typename NodeListT, int NTERMS, int NLAMBS, int NWAVES, int multiplicity=1>
class InteractionListManager
{

    typedef MultiHolder<multiplicity, BasePlaneWaveHolder<NLAMBS,NWAVES> > WaveT;

public:

    inline WaveT* wave(unsigned char nb_id, unsigned char src_id, DIRECTION dir) const
    {
        assert(nb_id < 27);
        assert(src_id >=1 && src_id <=8);
        return wave_ptrs[nb_id*8*6 + (src_id-1)*6 + dir];
    }

    typedef WaveT* PW_Ptr;
    inline PW_Ptr& wave(unsigned char nb_id, unsigned char src_id, DIRECTION dir)
    {
        assert(nb_id < 27);
        assert(src_id >=1 && src_id <=8);
        return wave_ptrs[nb_id*8*6 + (src_id-1)*6 + dir];
    }

    inline void set(const NodeT& node,
                    MultiHolder<6,WaveT> child_level_waves[],
                    const NodeListT& parent_level_nodes,
                    const NodeListT& child_level_nodes)
    {
        OctreeIndexer idxer = node.get_idx();

        // figure out which target child nodes actually really exist-
        // don't bother doing interaction accumulations on non-existent
        // child nodes (or shadows)!
        for (int child_id=1; child_id <= 8; ++child_id)
        {
            targs[child_id-1] = false;
            try {
                const NodeT& child_node = dynamic_cast<NodeT&>(node.get_child(child_id));
                if (child_node.isShadow() == false && child_node.isDeleted() == false) {
                    targs[child_id-1] = true;
                }
            } catch (BadIndexer) {}
            
        }

        for (unsigned short nb=0; nb < 27; ++nb)
        {
            if (nb == 13) { continue; } // self node not necessary

            // clear pointers
            for (unsigned short src=1; src <= 8; ++src)
            {
                wave(nb, src, UP) = NULL;
                wave(nb, src, DOWN) = NULL;
                wave(nb, src, NORTH) = NULL;
                wave(nb, src, SOUTH) = NULL;
                wave(nb, src, EAST) = NULL;
                wave(nb, src, WEST) = NULL;
            }

            OctreeIndexer nb_idxer;
            try
            {
                // get the neighbour idxer
                nb_idxer = idxer.get_neighbour_idxer(nb);
            }
            catch (BadIndexer) { continue; }

            typename NodeListT::const_iterator find_it = parent_level_nodes.find(nb_idxer);
            if (find_it != parent_level_nodes.end())
            {
                const NodeT* nb_node = find_it->second;
                if (nb_node->empty()) { continue; }
                for (unsigned short src=1; src <= 8; ++src)
                {
                    try {
                        const NodeT& other_node = dynamic_cast<const NodeT&>(nb_node->get_child(src));
                        if (other_node.empty()) { continue; }
                        long pw_idx = other_node.get_pw_idx();
                        wave(nb, src, UP) = &(child_level_waves[pw_idx](UP));
                        wave(nb, src, DOWN) = &(child_level_waves[pw_idx](DOWN));
                        wave(nb, src, NORTH) = &(child_level_waves[pw_idx](NORTH));
                        wave(nb, src, SOUTH) = &(child_level_waves[pw_idx](SOUTH));
                        wave(nb, src, EAST) = &(child_level_waves[pw_idx](EAST));
                        wave(nb, src, WEST) = &(child_level_waves[pw_idx](WEST));
                    }
                    catch (BadIndexer) {}
                }
            }
        }

    }

    template<int begin, int end>
    size_t do_interaction_list(MultiHolder<8,WaveT> targets[6], const fmm::Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs) const
    {
        
        // return the actual number of plane-wave translations carried out -- is interesting for profiling
        // (and, lets face it, isn't exactly going to massively affect the number of clock cycles in here...)
        size_t num_xlations = 0;

        // get interaction list for this neighbour_id
        #pragma omp parallel
        {
            
            for (unsigned short ctr=0; ctr < ILIST_SIZE; ++ctr)
            {
                const INTERACTION& inter = interactions[ctr];

                const WaveT* source = wave(inter.neighbour, inter.src, inter.dir);
                if (source == NULL) { continue; } // if the source node doesn't exist, nothing to do
                if (targs[inter.targ-1] == false) { continue; } // if the target (child) node doesn't exist, nothing to do

                WaveT& target = targets[inter.dir](inter.targ-1);

                switch(inter.dir)
                {
                case(UP):
                    fmm::translate_up_wave<multiplicity,begin,end>(level_globs,
                                        *source,
                                        target,
                                        inter.dx, inter.dy, inter.dz);
                    break;
                case(NORTH):
                    fmm::translate_up_wave<multiplicity,begin,end>(level_globs,
                                        *source,
                                        target,
                                        inter.dz, inter.dx, inter.dy);
                    break;
                case(EAST):
                    fmm::translate_up_wave<multiplicity,begin,end>(level_globs,
                                        *source,
                                        target,
                                        -inter.dz, inter.dy, inter.dx);
                    break;

                case(DOWN):
                    fmm::translate_down_wave<multiplicity,begin,end>(level_globs,
                                            *source,
                                            target,
                                            inter.dx, inter.dy, inter.dz);

                    break;
                case(SOUTH):
                    fmm::translate_down_wave<multiplicity,begin,end>(level_globs,
                                            *source,
                                            target,
                                            inter.dz, inter.dx, inter.dy);

                    break;
                case(WEST):
                    fmm::translate_down_wave<multiplicity,begin,end>(level_globs,
                                            *source,
                                            target,
                                            -inter.dz, inter.dy, inter.dx);

                    break;
                }
                
                ++num_xlations;
            }
        }
        
        return num_xlations;
    }

    size_t do_interaction_list(MultiHolder<8,WaveT> targets[6], const fmm::Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs) const
    {
        return do_interaction_list<0,multiplicity>(targets, level_globs);
    }

    size_t do_interaction_list_A(MultiHolder<8,WaveT> targets[6], const fmm::Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs) const
    {
        return do_interaction_list<0,8>(targets, level_globs);
    }
    size_t do_interaction_list_B(MultiHolder<8,WaveT> targets[6], const fmm::Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs) const
    {
        return do_interaction_list<8,12>(targets, level_globs);
    }

private:

    WaveT* wave_ptrs[216*6];
    bool targs[8];

};

// Inherit from Octree so that we can add extra FMM gubbins to the top-level Octree class
template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n=1>
class FMM_Octree_T : public Octree< FMM_NodeT<ContentType>, ContentType >
{

public:

    typedef ContentType CType;
    typedef Octree<FMM_NodeT<ContentType>, ContentType> super;
    typedef FMM_NodeT<ContentType> FMM_Node;
    typedef typename super::NodeList NodeList;

    typedef MultiHolder<n,BaseMultipoleHolder<NTERMS> > MultipoleHolder;
    typedef MultiHolder<n,BasePlaneWaveHolder<NLAMBS, NWAVES> > PlaneWaveHolder;

    FMM_Octree_T() {}
    FMM_Octree_T(unsigned int max_objects_per_node,
            const Vector& centre,
            double edge_length,
                 OctreeIndexer root_idx=OctreeIndexer(0,0,0,0)) :
                super(max_objects_per_node, centre, edge_length, root_idx) {}
    
    virtual ~FMM_Octree_T() {}

    unsigned int get_num_charges() const
    {
        return super::size();
    }

    void solve(double);
    
    inline const TimeInfo& get_timing_info() const { return timing_info; }
    inline void zero_timing_info() const { 
        timing_info.zero();
    }
   
    template <typename EvalPtType>
    void evaluate(EvalPtType& eval_pt)
    {
        const FMM_Node& node = super::get_node(eval_pt);
        
        long start_fmm_evals = myclock();
        evaluate_FMM(node, eval_pt);
        timing_info.fmm_evaluations += myclock() - start_fmm_evals;
        
        long start_explicit_evals = myclock();
        
        // where to put neighbourhood list
#ifdef __DELETED__
        const typename boost::shared_ptr<std::vector< ContentType* > >& neighbourhood = get_neighbourhood_contents(node);
#else
        const typename boost::shared_ptr<std::vector< ContentType* > >& neighbourhood = super::get_neighbourhood_contents(node);
#endif
        evaluate_explicit_neighbours(node, eval_pt, *neighbourhood);
        
        timing_info.explicit_evaluations += myclock() - start_explicit_evals;

        return;

    }
    
    double calculate_potential(const Vector&);
    void calculate_potential_and_field(const Vector&, double& pot, Vector&);
    
    template<typename EvalPtType>
    void evaluate_many(std::vector<EvalPtType*>& eval_pts)
    {
        
        // sanity check the input -- if nothing to do, there's nothing to do!
        if (eval_pts.size() == 0) { return; }
        if (super::size() == 0) { return; }

        long shared_start = myclock();

        typedef std::map< const FMM_Node*, std::vector<EvalPtType*> > NodeMap;
        NodeMap node_map;

        typedef std::vector<EvalPtType*> EvalPtList;
        // create a list of eval points grouped by octree leaf node
        for (typename EvalPtList::iterator it=eval_pts.begin(), end=eval_pts.end(); it != end; ++it)
        {
            EvalPtType* eval = *it;

            // where in the octree should this be evaluated?
            // convert x,y,z locations
            const FMM_Node* node = &(super::get_node(*eval));

            // find node address in map
            typename NodeMap::iterator find_it = node_map.find(node);
            if (find_it == node_map.end())
            {
                //node_map.insert(std::pair< const FMM_Node*, std::vector<EvalPtType*> >(node, std::vector<EvalPtType*>()) );
                node_map[node] = EvalPtList();
                find_it = node_map.find(node);
            }
            find_it->second.push_back(eval);

        }

        timing_info.explicit_evaluations += (myclock() - shared_start)/2;
        timing_info.fmm_evaluations += (myclock() - shared_start)/2;

        // Now we have a map of nodes to evaluation points
        for (typename NodeMap::iterator it=node_map.begin(), end=node_map.end(); it != end; ++it)
        {
            const FMM_Node* node = it->first;
            std::vector<EvalPtType*>& eval_pts_in_node = it->second;
            
            // loop over eval points
            long start_fmm = myclock();
            if (node->get_level() >= 2)
            {
                #pragma omp parallel for
                for (int ii=0; ii < eval_pts_in_node.size(); ++ii)
                {
                    evaluate_FMM(*node, *(eval_pts_in_node[ii]));
                }
            }
            timing_info.fmm_evaluations += myclock() - start_fmm;

            // loop over eval points
            long start_explicit = myclock();
            typedef boost::shared_ptr< std::vector< ContentType* > > NBListPtr;
            NBListPtr neighbourhood = super::get_neighbourhood_contents(*node);
            
            // empty nodes --> nothing to do.
            if (neighbourhood->size() == 0) {
                continue;
            }
            
            #pragma omp parallel for
            for (int ii=0; ii < eval_pts_in_node.size(); ++ii)
            {
                evaluate_explicit_neighbours(*node, *(eval_pts_in_node[ii]), *neighbourhood);
            }
            timing_info.explicit_evaluations += myclock() - start_explicit;

        }

        return;

    }

#ifdef OPENCL

    template<typename EvalPtType>
    void evaluate_many(std::vector<EvalPtType*>& eval_pts, OpenCL_Handler& ocl_handler)
    {
        // sanity check the input -- if nothing to do, there's nothing to do!
        if (eval_pts.size() == 0) { return; }
        if (super::size() == 0) { return; }

        long shared_start = myclock();

        typedef std::map< FMM_Node*, std::vector<EvalPtType*> > NodeMap;
        NodeMap node_map;

        typedef std::vector<EvalPtType*> EvalPtList;
        // create a list of eval points grouped by octree leaf node
        for (typename EvalPtList::iterator it=eval_pts.begin(), end=eval_pts.end(); it != end; ++it)
        {
            EvalPtType* eval = *it;

            // where in the octree should this be evaluated?
            // convert x,y,z locations
            FMM_Node* node = &(super::get_node(*eval));

            // find node address in map
            typename NodeMap::iterator find_it = node_map.find(node);
            if (find_it == node_map.end())
            {
                //node_map.insert(std::pair< const FMM_Node*, std::vector<EvalPtType*> >(node, std::vector<EvalPtType*>()) );
                node_map[node] = EvalPtList();
                find_it = node_map.find(node);
            }
            find_it->second.push_back(eval);

        }
        
        timing_info.explicit_evaluations += (myclock() - shared_start)/2;
        timing_info.fmm_evaluations += (myclock() - shared_start)/2;

        long start_ocl_explicit = myclock();
        // OpenCL Max memory constraint -- to do the FMM explicit interactions in OpenCL
        // we need to allocate a chunk of global memory for thread-block level results.
        // For large FMM problems it is quite possible that we will exceed the maximum malloc
        // size (e.g. 256MB on GTX280) so need to break the problem into chunks if we find
        // that we are approaching the limit.
        size_t opencl_max_malloc = ocl_handler.get_max_malloc_size();

        // Now we have a map of nodes to evaluation points
        for (typename NodeMap::iterator it=node_map.begin(), end=node_map.end(); it != end; ++it)
        {
            const FMM_Node* node = it->first;
            std::vector<EvalPtType*>& eval_pts_in_node = it->second;

            typedef boost::shared_ptr< std::vector< ContentType* > > NBListPtr;
            NBListPtr neighbourhood = super::get_neighbourhood_contents(*node);

            // break problem into smaller chunks if necessary
            unsigned int num_x_chunks = FMM_Resources<EvalPtType>::calc_x_chunksize(eval_pts_in_node.size(), neighbourhood->size(), opencl_max_malloc);
            unsigned int num_y_chunks = FMM_Resources<EvalPtType>::calc_y_chunksize(eval_pts_in_node.size(), neighbourhood->size(), opencl_max_malloc);

            // if number of chunks is 1, then no need to split the problem up
            if (num_x_chunks == 1 && num_y_chunks == 1)
            {
                // IMPORTANT: this is heap allocated and then passed to the OpenCL handler which process the
                // work in a separate thread-- it is the job of the OpenCL handler to correctly delete this!
                FMM_Resources<EvalPtType>* fmm_res_ptr = new FMM_Resources<EvalPtType>(eval_pts_in_node, *neighbourhood, kappa_explicit);
                ocl_handler.add_work_to_queue(fmm_res_ptr);
            }
            else
            {
                // Annoyingly complicated: break the neighbourhood into equal sized lumps, and
                // create an FMM work unit for each
                size_t x_chunk_size = static_cast<size_t>(ceil(static_cast<double>(eval_pts_in_node.size()) / num_x_chunks));
                size_t y_chunk_size = static_cast<size_t>(ceil(static_cast<double>(neighbourhood->size()) / num_y_chunks));

                typename std::vector<EvalPtType*>::const_iterator eval_it=eval_pts_in_node.begin(), eval_end=eval_pts_in_node.end();
                while(eval_it != eval_end)
                {
                    std::vector<EvalPtType*> eval_pts_chunk;
                    eval_pts_chunk.reserve(x_chunk_size);
                    for (size_t x_chunk_ctr=0 ; eval_it != eval_end && x_chunk_ctr < x_chunk_size; ++eval_it)
                    {
                        eval_pts_chunk.push_back(*eval_it);
                        ++x_chunk_ctr;
                    }

                    typename std::vector< ContentType* >::const_iterator neigh_it=neighbourhood->begin(), neigh_end=neighbourhood->end();
                    while(neigh_it != neigh_end)
                    {
                        std::vector< ContentType* > neighbourhood_chunk;
                        neighbourhood_chunk.reserve(y_chunk_size);
                        for (size_t y_chunk_ctr=0 ; neigh_it != neigh_end && y_chunk_ctr < y_chunk_size; ++neigh_it)
                        {
                            neighbourhood_chunk.push_back(*neigh_it);
                            ++y_chunk_ctr;
                        }
                        FMM_Resources<EvalPtType>* fmm_res_ptr = new FMM_Resources<EvalPtType>(eval_pts_chunk, neighbourhood_chunk, kappa_explicit);
                        ocl_handler.add_work_to_queue(fmm_res_ptr);

                    }
                }
            }

        }
        timing_info.explicit_evaluations += (myclock() - start_ocl_explicit);

        long start_fmm = myclock();
        for (typename NodeMap::iterator it=node_map.begin(), end=node_map.end(); it != end; ++it)
        {
            const FMM_Node* node = it->first;
            if (node->get_level() < 2) { continue; }

            std::vector<EvalPtType*>& eval_pts_in_node = it->second;

            #pragma omp parallel for
            for (int ii=0; ii < eval_pts_in_node.size(); ++ii)
            {
                evaluate_FMM(*node, *(eval_pts_in_node[ii]));
            }
        }
        timing_info.fmm_evaluations += (myclock() - start_fmm);


        // wait for the GPU FMM evals to complete
        //std::cout << "Waiting for FMM local interactions on GPU" << std::endl;
        long waiting_explicits=myclock();
        while (ocl_handler.pending() > 0)
        {
            boost::this_thread::sleep(boost::posix_time::milliseconds(10));
        }
        timing_info.explicit_evaluations += myclock() - waiting_explicits;
        //std::cout << "FMM local interactions are complete" << std::endl;

        return;
    }
#endif

    inline void add(boost::shared_ptr<ContentType>& ptr) {
        super::add(ptr);
        return;
    }

    inline void add(const ContentType& ch)
    {
        super::add( ch );
        return;
    }

    inline void add(const std::vector<ContentType>& things) {
        super::add(things);
        return;
    }


private:

    // The all-important screening factor, kappa
    double kappa; // guaranteed non-zero- at least 1e-10 (otherwise special functions return Nans)
    double kappa_explicit; // may be zero - used for the explicit calcs
    
    template <typename EvalPtType>
    void evaluate_FMM(const FMM_Node& node,
                      EvalPtType& eval_pt) const
    {
        // If the node is above the second level then the multipole expansion will be zero
        if (node.get_level() < 2) {
            return;
        }

        double pot_local = 0;
        Vector field_local(0,0,0);

        const long local_idx = node.get_lc_idx();

		// There may be no levels if this is a simple sphere?
		if (lexp_arrays.find(node.get_level()) == lexp_arrays.end()) return;

		// Otherwise, evaluate the local expansions...
        boost::shared_array<MultipoleHolder> local_expansions = lexp_arrays.find(node.get_level())->second;
        const MultipoleHolder& locals = local_expansions[local_idx];

        double scale = FMM_Globals<NTERMS>::get_scale(kappa, super::universe_edge_length, node.get_level());
        evaluate_local_expansion_at_xyz(kappa,
                                        scale,
                                        eval_pt,
                                        node.get_centre(),
                                        locals);
        return;
    }

    void evaluate_explicit_neighbours(const FMM_Node& node,
                                      EvalPtBase& eval_pt,
                                      const std::vector< ContentType* >& charges_and_locations) const
    {
        double pot=0;
        Vector field(0,0,0);

        for (typename std::vector< ContentType* >::const_iterator it=charges_and_locations.begin(), end=charges_and_locations.end();
                it != end;
                ++it)
        {
            EvalPtBase::add_explicit_contrib(eval_pt, **it, pot, field, kappa_explicit);
        }

        eval_pt.add_field(field);
        eval_pt.add_potential(pot);
        return;
    }

    void evaluate_explicit_neighbours(const FMM_Node& node,
                                      EvalPt_2ndDerivs& eval_pt,
                                      const std::vector< ContentType* >& charges_and_locations) const
    {
        double pot=0;
        Vector field(0,0,0);
        GradField3x3 field2(0.0);

        for (typename std::vector< ContentType* >::const_iterator it=charges_and_locations.begin(), end=charges_and_locations.end();
                it != end;
                ++it)
        {
            EvalPt_2ndDerivs::add_explicit_contrib(eval_pt, **it, pot, field, field2, kappa_explicit);
        }

        eval_pt.add_potential(pot);
        eval_pt.add_field(field);
        eval_pt.add_field2(field2);
        return;
    }

    void calculate_multipole_expansions();
    void translate_multipole_expansions();
    void local_expansions(std::map<unsigned short, long> &);
    void allocate_memory(std::map<unsigned short, long>&);

    // a simple profiling class
    mutable TimeInfo timing_info;
    
    std::map<unsigned short, boost::shared_array<MultipoleHolder> > mexp_arrays;
    std::map<unsigned short, boost::shared_array<MultipoleHolder> > lexp_arrays;

    // FMM-related arrays (rotation matrices, scale factors...)
    FMM_Globals<NTERMS> fmm_globs;

    struct MpIdxSorter
    {
        inline bool operator()(const typename super::NodeT* ptr1, const typename super::NodeT* ptr2) const
        {
            return ptr1->get_mp_idx() < ptr2->get_mp_idx();
        }
    };
    
};

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::allocate_memory(std::map<unsigned short, long>& pw_ctrs)
{

    long start_allocate_memory = myclock();
    
    // clear any existing multipole/local expansion data
    mexp_arrays.clear();
    lexp_arrays.clear();

    // traverse the tree and figure out how much of each type of holder we will require.

    // TODO: just use arrays- std::vector is overkill
    // init. counters
    // also are unsigned longs really necessary??
    pw_ctrs.clear();
    std::map<unsigned short, long> mp_ctrs;
    std::map<unsigned short, long> lc_ctrs;

    for (int i=super::get_top_level(); i <= super::get_bottom_level(); ++i)
    {
        pw_ctrs[i] = 0;
        mp_ctrs[i] = 0;
        lc_ctrs[i] = 0;
    }

    //std::cout << "Traversing tree." << std::endl;
    for (int level=super::get_bottom_level(); level >= super::get_top_level(); --level)
    {
        long& pw_ctr = pw_ctrs[level];
        long& mp_ctr = mp_ctrs[level];
        long& lc_ctr = lc_ctrs[level];

        // loop over cubes in level
        for (typename NodeList::iterator it=super::get_node_list(level).begin(), end=super::get_node_list(level).end();
            it != end;
            ++it)
        {
            FMM_Node& node = *(it->second);

            // ALL nodes need multipole expansions (even apparently deleted ones for now...)
            if (!node.empty())
            {
                node.set_mp_idx(mp_ctr++);
            }
            
            // that's it for deleted nodes though
            if (node.isDeleted()) { continue; }

            // set plane-wave holder index -- this is to conserve memory during local expansion downward pass
            // note that we omit some nodes (e.g. ones below the leaf nodes in an adaptive tree) because we
            // don't need the local expansions that far down the tree (because it is cheap enough to evaluate
            // the potential at a higher level in the tree.)
            if (node.isShadow() == false)
            {
                node.set_lc_idx(lc_ctr++);
            }

            // we need multipole expansions in all non-empty nodes of the tree
            if (!node.empty())
            {
                node.set_pw_idx(pw_ctr++);
            }
        }
    }

    // Allocate some memory

    // create Multipole Expansion holders (and Local Expansion and Plane Wave holders).
    //std::cout << "Allocating memory." << std::endl;

    for (int level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        if (mp_ctrs[level] > 0)
        {
            boost::shared_array<MultipoleHolder> mpole(new MultipoleHolder[mp_ctrs[level]]);
            assert(mpole.get() != NULL);
            mexp_arrays[level] = mpole;
            memset(mpole.get(), 0, sizeof(MultipoleHolder)*mp_ctrs[level]);
        }

        if (lc_ctrs[level] > 0) 
        {
            boost::shared_array<MultipoleHolder> local(new MultipoleHolder[lc_ctrs[level]]);
            assert(local.get() != NULL);
            lexp_arrays[level] = local;
            memset(local.get(), 0, sizeof(MultipoleHolder)*lc_ctrs[level]);
        }
    }

    // add time taken to timing_info block
    timing_info.allocate_memory += (myclock() - start_allocate_memory);

    return;

}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::solve(double kappa_in)
{
    kappa_explicit = kappa_in;
    kappa = (kappa_in < 1e-10) ? 1e-10 : kappa_in;
    
    // note the time taken...
    long start_init = myclock();

    // ensure the neighbourhoods around cells are correct
    if (super::built_neighbourhoods == false)
    {
        super::build_neighbourhoods();
        //super::optimize_neighbourhoods();
        super::remove_empty_nodes();
    }

    // this holds the number of cubes which will require plane wave holders
    std::map<unsigned short, long> plane_wave_counters;
    allocate_memory(plane_wave_counters); // re-indexes the tree and allocates memory

    // .. in the profile counter
    timing_info.init += myclock() - start_init;
    
    if (get_num_charges() == 0) { return; }

    //std::cout << "Starting upward pass of tree." << std::endl;
    long start_upward_pass = myclock();

    // calculate multipole expansions at lowest level
    calculate_multipole_expansions();

    // upward pass of multipole expansions
    translate_multipole_expansions();
    
    // profile
    timing_info.total_upward_pass += (myclock() - start_upward_pass);

    //std::cout << "Starting downward pass of tree." << std::endl;

    // downward pass of local expansions
    // (this function is massive)
    long start_downward_pass = myclock();
    local_expansions(plane_wave_counters);
    timing_info.total_downward_pass += (myclock() - start_downward_pass);

    // clear the multipole array memory (not needed anymore)
    // NB: the local expansions are still allocated as they are
    // required to actually evaluate stuff.
    long start_dealloc = myclock();
    mexp_arrays.clear();
    timing_info.allocate_memory += myclock() - start_dealloc;

    // done
    //std::cout << "Done downward pass." << std::endl;

    return;
}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::calculate_potential_and_field(const Vector& xyz, double& pot, Vector& f)
{
    EvalPt eval_pt(xyz);
    evaluate(eval_pt);
    f = eval_pt.get_field();
    pot = eval_pt.get_potential();

    return;

}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
double FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::calculate_potential(const Vector& xyz)
{
    EvalPt eval_pt(xyz);
    evaluate(eval_pt);
    return eval_pt.get_potential();

}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::calculate_multipole_expansions()
{
    // start timer for profiling
    long start_calculate_mpoles = myclock();

    for (unsigned short level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        boost::shared_array< MultipoleHolder > multipoles = mexp_arrays[level];
        double scale = FMM_Globals<NTERMS>::get_scale(kappa, super::universe_edge_length, level);

        const NodeList& node_list = super::get_node_list(level);
        
        // workspace for legendre polynomials
        LegendreHolder<NTERMS> workspace_p;
        
        std::vector<typename super::NodeT*> ptrs;
        for (typename super::const_childless_iterator childless_it=super::childless_begin(level), childless_end=super::childless_end(level);
            childless_it != childless_end;
            ++childless_it)
        {
            ptrs.push_back(childless_it->second);
        }
        std::sort(ptrs.begin(), ptrs.end(), MpIdxSorter()); // default sort to arrange in order they appear in memory.  Why not.

        for (int ii=0; ii < ptrs.size(); ++ii)
        {
            const FMM_Node& node = *(ptrs[ii]);

            MultipoleHolder& mexp = multipoles[node.get_mp_idx()];
            yformmp(kappa,
                    node.get_centre(),
                    node.get_contents(),
                    mexp,
                    scale,
                    workspace_p,
                    fmm_globs.scale_factors);

        }
    }
    
    // profile
    timing_info.create_mpoles += myclock() - start_calculate_mpoles;

    return;
}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::translate_multipole_expansions()
{

    // start timer for profiling
    long start_xlating_mpoles = myclock();
    
    // upwardly translate the expansions to centres of parent cubes
    for (int from_level=super::get_bottom_level(); from_level > super::get_top_level(); --from_level)
    {

        // gonna need some precalcs to do the translations...
        //std::cout << "Pre-calculating MP shift coefficients for level " << from_level << std::endl;
        Level_Dependent_FMM_Globals<NTERMS,NLAMBS> fmm_level_globs(kappa, from_level, super::universe_edge_length);
        const DblMatrix3D& mp_translation_coefficients = fmm_level_globs.multipole_shift_coefficients;

        // do the translations
        //std::cout << "Translating multipoles from level " << from_level << " to level " << from_level-1 << std::endl;
        boost::shared_array<MultipoleHolder>& parent_multipoles = mexp_arrays[from_level-1];

        // loop over cubes in level
        boost::shared_array<MultipoleHolder>& child_multipoles = mexp_arrays[from_level];

        // heap allocate some workspaces
        BaseMultipoleHolder<NTERMS> workspace1;
        BaseMultipoleHolder<NTERMS> workspace2;

        // loop over cubes in level
        std::vector<FMM_Node*> ptrs;
        ptrs.reserve(super::get_node_list(from_level).size());
        for (typename NodeList::iterator it=super::get_node_list(from_level).begin(), end=super::get_node_list(from_level).end();
            it != end;
            ++it)
        {
            FMM_Node& node = *(it->second);

            // if this node has nothing in it then there's nothing to translate upwards except
            // a whole bunch of zeros.
            if (node.empty()) { continue; }
            ptrs.push_back(it->second);
            
        }
        
        #pragma omp parallel for private(workspace1,workspace2)
        for (int ii=0; ii < ptrs.size(); ++ii)
        {
            FMM_Node& node = *(ptrs[ii]);

            // get parent node
            FMM_Node& parent_node = dynamic_cast<FMM_Node&>(*node.parent_ptr);

            // generate slice of multipoles matrix
            MultipoleHolder& child_mexp  = child_multipoles[node.get_mp_idx()];
            MultipoleHolder& parent_mexp = parent_multipoles[parent_node.get_mp_idx()];

            int id = node.get_id_within_parent();
            const DblMatrix3D& rotation_matrix =  (id > 4) ? fmm_globs.rdmsq3 : fmm_globs.rdsq3;

            ympshift(   id,
                        child_mexp,
                        parent_mexp,
                        workspace1,
                        workspace2,
                        mp_translation_coefficients,
                        rotation_matrix);

        }
    }

    // profile
    timing_info.pass_mpoles_upwards += myclock() - start_xlating_mpoles;

    return;
}

template <typename ContentType, int NTERMS, int NLAMBS, int NWAVES, int n>
void FMM_Octree_T<ContentType, NTERMS, NLAMBS, NWAVES, n>::local_expansions(std::map<unsigned short, long>& plane_wave_counters)
{
    for (unsigned short level=super::get_top_level(); level < super::get_bottom_level(); ++level)
    {
        // some workspace
        BaseMultipoleHolder<NTERMS> workspace1;
        BaseMultipoleHolder<NTERMS> workspace2;

        long start_inheritance = myclock();
        
        //std::cout << "Pre-calculating local shift coefficients for level " << level << std::endl;
        Level_Dependent_FMM_Globals<NTERMS,NLAMBS> fmm_level_globs(kappa, level, super::universe_edge_length);

        //std::cout << "Translating locals downward from level " << level << " to level " << level+1 << std::endl;

        //CmplxMatrix3D& this_level_multipoles = mexp_arrays[level];
        boost::shared_array<MultipoleHolder> this_level_locals = lexp_arrays[level];
        boost::shared_array<MultipoleHolder> child_level_multipoles = mexp_arrays[level+1];
        boost::shared_array<MultipoleHolder> child_level_locals = lexp_arrays[level+1];

        // loop over cubes in level
        typename super::NodeList& this_level_node_list = super::get_node_list(level);
        std::vector<typename super::NodeT*> ptrs;
        ptrs.reserve(this_level_node_list.size());
        for (typename NodeList::iterator it=this_level_node_list.begin(), end=this_level_node_list.end();
            it != end;
            ++it)
        {
            FMM_Node& node = *(it->second);
            
            // if this is a leaf node, then this is the last point in the hierarchy at which we need local expansions,
            // so don't bother transmitting them to the next level down
            if (node.isDeleted() || node.isLeaf() || node.isShadow()) { continue; }
            
            ptrs.push_back(it->second);
        }
        
        #pragma omp parallel for private(workspace1, workspace2)
        for (int ii=0; ii < ptrs.size(); ++ii)
        {
            FMM_Node& node = *(ptrs[ii]);

            // generate slice of multipoles matrix
            MultipoleHolder& this_lexp = this_level_locals[node.get_lc_idx()];

            for (int child=1; child <= 8; ++child)
            {
                try {
                    const FMM_Node& child_node = dynamic_cast<const FMM_Node&>(node.get_child(child));
                    if (child_node.isShadow() || child_node.isDeleted()) { continue; }

                    // child local expansion holder
                    const long lc_idx = child_node.get_lc_idx();
                    MultipoleHolder& child_lexp = child_level_locals[lc_idx];

                    const DblMatrix3D& rotation_matrix = (child > 4) ? fmm_globs.rdsq3 : fmm_globs.rdmsq3;

                    ylcshift(   child,
                                this_lexp,
                                child_lexp,
                                workspace1,
                                workspace2,
                                fmm_level_globs.local_shift_coefficients,
                                rotation_matrix);
                } catch (BadIndexer) {}
            }
        }
        
        // profile
        timing_info.inherit_locals_downwards += myclock() - start_inheritance;

        long start_malloc = myclock();

        // planewave holders
        // big fat malloc for the plane wave expansion holders
        // could do this earlier since PlaneWaveHolder is made large enough
        // for maximum nexptotp
        long plane_waves_needed = plane_wave_counters[level+1];
        const double mb_size = plane_waves_needed * sizeof(MultiHolder<6,PlaneWaveHolder>)  / double(1024*1024);
        //std::cout << "Allocating memory for plane waves... (" << mb_size << " Mbytes)" << std::endl;
        boost::shared_array< MultiHolder<6,PlaneWaveHolder> > child_level_planewaves(new MultiHolder<6,PlaneWaveHolder>[plane_waves_needed]);
        memset(child_level_planewaves.get(), 0, sizeof(MultiHolder<6,PlaneWaveHolder>)*plane_waves_needed);

        // kindof counts as allocating memory I reckon, so add it to the profile
        timing_info.allocate_memory += myclock() - start_malloc;

        long start_convert_to_pw = myclock();
        
        //std::cout << "Converting multipoles to plane waves for children of level " << level << std::endl;
        NodeList& node_list = super::get_node_list(level+1);

        MultipoleHolder mpole_rotated;
        
        std::vector<FMM_Node*> child_ptrs;
        for (typename NodeList::iterator child_it=node_list.begin(), child_end=node_list.end(); child_it != child_end; ++child_it)
        {
            FMM_Node& child_node = *(child_it->second);

            // empty nodes will not have anything interesting to transmit to their neighbours
            if (child_node.isDeleted() || child_node.empty()) { continue; }

            child_ptrs.push_back(&child_node);
            
        }
        
        #pragma omp parallel for private(mpole_rotated)
        for (int ii=0; ii < child_ptrs.size(); ++ii)
        {
            FMM_Node& child_node = *(child_ptrs[ii]);
            const OctreeIndexer& idxer = child_node.get_idx();

            long pw_idx = child_node.get_pw_idx();
            assert(pw_idx < plane_waves_needed);
            PlaneWaveHolder& up = child_level_planewaves[pw_idx](UP);
            PlaneWaveHolder& down = child_level_planewaves[pw_idx](DOWN);
            PlaneWaveHolder& north = child_level_planewaves[pw_idx](NORTH);
            PlaneWaveHolder& south = child_level_planewaves[pw_idx](SOUTH);
            PlaneWaveHolder& east = child_level_planewaves[pw_idx](EAST);
            PlaneWaveHolder& west = child_level_planewaves[pw_idx](WEST);

            MultipoleHolder& child_mpole = child_level_multipoles[child_node.get_mp_idx()];

            // convert multipole expansion to planewave expansions
            convert_mp_to_exp(fmm_globs,
                            fmm_level_globs,
                            child_mpole,
                            up,
                            down);

            rotate_ztoy(child_mpole,
                        mpole_rotated,
                        fmm_globs.rdmpi2);

            // convert multipole expansion to planewave expansions
            convert_mp_to_exp( fmm_globs,
                                fmm_level_globs,
                                mpole_rotated,
                                north,
                                south);

            // rotate multipole expansion so that plane wave is translated in z-direction
            rotate_ztox(child_mpole,
                        mpole_rotated,
                        fmm_globs.rdpi2);

            // convert multipole expansion to planewave expansions
            convert_mp_to_exp( fmm_globs,
                                fmm_level_globs,
                                mpole_rotated,
                                east,
                                west);

        }
        
        // add time taken to rotate/convert to profile counter
        timing_info.convert_to_planewaves += myclock() - start_convert_to_pw;

        //std::cout << "Processing interaction lists for children of level " << level << " (this might be slow...)" << std::endl;
        typedef InteractionListManager<typename super::NodeT, typename super::NodeList, NTERMS, NLAMBS, NWAVES, n> InteractionList;

        const NodeList& parent_nodes = super::get_node_list(level);
        const NodeList& child_nodes = super::get_node_list(level+1);

        // workspace for converted plane-wave --> local expansion (heap allocated)
        //boost::scoped_ptr<MultipoleHolder> local_out_workspace_ptr(new MultipoleHolder);
        
        MultipoleHolder local_out;
        boost::scoped_ptr<InteractionList> ilm;
        
        // loop over cubes in level            
        
        for (int ii=0; ii < parent_nodes.size(); ++ii)
        {
            boost::shared_array< MultiHolder<8,PlaneWaveHolder> > accum(new MultiHolder<8,PlaneWaveHolder>[6]);
            
            if (!ilm) {
                ilm.reset(new InteractionList);
            }
            
            typename NodeList::const_iterator it=parent_nodes.begin();
            for (int jj=0; jj < ii; ++jj) { ++it; }
            
            const FMM_Node& node = *(it->second);

            // IMPORTANT: if is a leaf node, then can evaluate the FMM here without further subdivision;
            // do not carry out any translating on the lower level, as it is not necessary. If the node is
            // deleted then definitely don't bother; if it's a shadow then it is below a leaf node, so 
            // the higher level expansion is all that's needed, not those in the shadows.
            if (node.isDeleted() || node.isLeaf() || node.isShadow()) { continue; }

            // time taken for interaction list
            long start_ilist = myclock();
            
            memset(accum.get(), 0, sizeof(MultiHolder<8,PlaneWaveHolder>)*6);
            ilm->set(node, child_level_planewaves.get(), parent_nodes, child_nodes);
            size_t num_ilist_xlations = ilm->do_interaction_list(accum.get(), fmm_level_globs);
            
            // profile time taken for ilist
            timing_info.interaction_list_xlations += myclock() - start_ilist;
            
            // stick the number of ilist-xlations in too for good measure
            timing_info.num_ilist_xlations += num_ilist_xlations;

            // start conversion back again to local expansions
            long start_convert_from_pw = myclock();
            
            #pragma omp parallel for private(mpole_rotated,local_out)
            for (int child_id=1; child_id <=8; ++child_id)
            {
                try {
                    const FMM_Node& child_node = dynamic_cast<const FMM_Node&>(node.get_child(child_id));

                    if(child_node.isShadow() || child_node.isDeleted()) { continue; }

                    long lc_idx = child_node.get_lc_idx();
                    MultipoleHolder& child_lexp = child_level_locals[lc_idx];

                    // unravel plane wave expansion holders
                    const PlaneWaveHolder& up = accum[UP](child_id-1);
                    const PlaneWaveHolder& down = accum[DOWN](child_id-1);
                    const PlaneWaveHolder& north = accum[NORTH](child_id-1);
                    const PlaneWaveHolder& south = accum[SOUTH](child_id-1);
                    const PlaneWaveHolder& east = accum[EAST](child_id-1);
                    const PlaneWaveHolder& west = accum[WEST](child_id-1);

                    // UP/DOWN
                    pw_to_local(fmm_globs,
                                fmm_level_globs,
                                up,
                                down,
                                local_out);
                    child_lexp += local_out;

                    // NORTH/SOUTH
                    pw_to_local(fmm_globs,
                                fmm_level_globs,
                                north,
                                south,
                                local_out);

                    rotate_ytoz(local_out,
                                mpole_rotated,
                                fmm_globs.rdpi2);
                    child_lexp += mpole_rotated;

                    // EAST/WEST
                    pw_to_local(fmm_globs,
                                fmm_level_globs,
                                east,
                                west,
                                local_out);

                    rotate_ztox(local_out,
                                mpole_rotated,
                                fmm_globs.rdmpi2);
                    child_lexp += mpole_rotated;
                }
                catch (BadIndexer) {}
            }
            
            // add time taken from rotate/convert to profile counter
            timing_info.convert_from_planewaves += myclock() - start_convert_from_pw;

        }
        
    }

    return;
}

// Specialisation: the basic-issue FMM Octree takes Charge objects and calculates
// potentials, fields, and whatever as one would expect.
typedef FMM_Octree_T<Charge,9,9,67> FMM_Octree_3FIG_ACCURACY;
typedef FMM_Octree_T<Charge,18,18,300> FMM_Octree_6FIG_ACCURACY;

} // end fmm namespace


#endif /* FMM_OCTREE_H_ */
