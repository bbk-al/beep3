/*
* fmm_bem_hybrid.h
*
*  Created on: 16 Aug 2010
*      Author: david
*/

#ifndef FMM_BEM_HYBRID_H_
#define FMM_BEM_HYBRID_H_

#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <exception>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "../common/octree_indexer.h"
#include "../common/matrix_types.h"
#include "fmm_globals.h"
#include "fmm_octree.h"
#include "fmm_time_info.h"
#include "../bem/node_patch.h"

#include "interaction_list.h"

namespace fmm
{

template <int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
class FMM_BEM_Octree_Node : public FMM_NodeT<NodePatchT>
{

    typedef FMM_NodeT<NodePatchT> super;

public:

    typedef MultiHolder<12, BaseMultipoleHolder<NTERMS> > MultipoleHolder;
    typedef MultiHolder<12, BasePlaneWaveHolder<NLAMBS, NWAVES> > PlaneWaveHolder;

public:

    FMM_BEM_Octree_Node(const Vector& _centre, double _edge_length, const OctreeIndexer& root_idx) :
        super(_centre, _edge_length, root_idx) {}
    FMM_BEM_Octree_Node(FMM_NodeT<NodePatchT>& other, OctreeIndexer idxer, unsigned short child_id) :
        super(other, idxer, child_id) {}
    virtual ~FMM_BEM_Octree_Node() {}

};

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
class FMM_BEM_Octree : public Octree< FMM_BEM_Octree_Node<NTERMS, NLAMBS, NWAVES, NodePatchT>, NodePatchT >
{
public:

    typedef Octree< FMM_BEM_Octree_Node<NTERMS, NLAMBS, NWAVES, NodePatchT>, NodePatchT > super;
    typedef typename super::NodeT NodeT;

    std::map<unsigned short, size_t> plane_wave_counters;
    typedef typename NodeT::MultipoleHolder MultipoleHolder;
    typedef typename NodeT::PlaneWaveHolder PlaneWaveHolder;

private:

    TimeInfo timing_info; // simple timing class for (rough) profiling
    
    fmm::FMM_Globals<NTERMS> fmm_globs;
    typedef std::map<unsigned short, boost::shared_ptr<const Level_Dependent_FMM_Globals<NTERMS, NLAMBS> > > FMMLevelGlobalsMap;
    FMMLevelGlobalsMap fmm_level_globs_beta;
    FMMLevelGlobalsMap fmm_level_globs_beta0;

    std::map<unsigned short, boost::shared_array<MultipoleHolder> > mexp_arrays;
    std::map<unsigned short, boost::shared_array<MultipoleHolder> > lexp_arrays;

    // FMM functions
    void calculate_multipole_expansions(double beta, double beta0);
    void translate_multipole_upward(double beta, double beta0);
    void downward_pass(double beta, double beta0);

    void interaction_lists(double beta, double beta0);
    void inherit_downwards(double beta, double beta0);

    void pass_local_expansion_downwards(NodeT& octree_node,
                                        const fmm::FMM_Globals<NTERMS>& fmm_globs,
                                        const fmm::Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& level_globs_beta,
                                        const fmm::Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& level_globs_beta0)
    {

        BaseMultipoleHolder<NTERMS> workspace1;
        BaseMultipoleHolder<NTERMS> workspace2;

        for (int child_octant_id=1; child_octant_id <= 8; ++child_octant_id)
        {
            try {
                // this willthrow exception if the node doesn't exist
                const NodeT& child_node = super::get_node(octree_node.get_child_idx(child_octant_id));
                if (child_node.empty() || child_node.isShadow()) { continue; }
                MultipoleHolder& local_out = get_lpole(child_node);

                const DblMatrix3D& rotation_matrix = (child_octant_id > 4) ? fmm_globs.rdsq3 : fmm_globs.rdmsq3;

                fmm::ylcshift<12,0,8>(child_octant_id,
                                get_lpole(octree_node),
                                local_out,
                                workspace1,
                                workspace2,
                                level_globs_beta.local_shift_coefficients,
                                rotation_matrix);
                fmm::ylcshift<12,8,12>(child_octant_id,
                                get_lpole(octree_node),
                                local_out,
                                workspace1,
                                workspace2,
                                level_globs_beta0.local_shift_coefficients,
                                rotation_matrix);
            }
            catch (BadIndexer) {}
        }

        return;
    }

public:

    typedef std::map< OctreeIndexer, std::vector<NodeT*> > TransmitList;
    typedef std::pair< OctreeIndexer, std::vector<NodeT*> > TransmitListPair;

    // list of nodes which must be transmitted across chares
    TransmitList chare_transmissions;

    std::vector<OctreeIndexer> subdivides_across_chare_boundary;

    FMM_BEM_Octree() {}
    virtual ~FMM_BEM_Octree() {}

    FMM_BEM_Octree(unsigned int max_items_per_node, const Vector& centre, double edge_length, OctreeIndexer root_idx_=OctreeIndexer(0,0,0,0)) :
        super(max_items_per_node, centre, edge_length, root_idx_)
    {}

    void allocate_memory();

    void solve(double beta, double beta0, const double f_lhs[], const double h_lhs[], double f_rhs[], double h_rhs[]);
    void evaluate(double beta, double beta0, double f_rhs[], double h_rhs[]);

    inline const TimeInfo& get_timing_info() const { return timing_info; }
    inline TimeInfo& get_timing_info() { return timing_info; }
#ifdef __DELETED__
    inline void zero_timing_info() const { 
#else
    inline void zero_timing_info() { 
#endif
        timing_info.zero();
    }
    
    inline bool is_within_this_chare(const OctreeIndexer& idx) const
    {
        // if the passed idx when followed upwards to the chare level matches
        // the chare index, then the idx falls within the child-space of the chare,
        // and hence of this octree segment.
        return idx.translate_upwards_to_level(super::root_level) == super::root_idx;
    }

    inline const MultipoleHolder& get_mpole(const NodeT& node) const
    {
        return mexp_arrays.find(node.get_level())->second[node.get_mp_idx()];
    }

    inline MultipoleHolder& get_mpole(const NodeT& node)
    {
        return mexp_arrays.find(node.get_level())->second[node.get_mp_idx()];
    }

    inline const MultipoleHolder& get_lpole(const NodeT& node) const
    {
        return lexp_arrays.find(node.get_level())->second[node.get_lc_idx()];
    }

    inline MultipoleHolder& get_lpole(const NodeT& node)
    {
        return lexp_arrays.find(node.get_level())->second[node.get_lc_idx()];
    }
};

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::solve(double beta, double beta0, const double f_lhs[], const double h_lhs[], double f_rhs[], double h_rhs[])
{
    if (super::get_bottom_level() < 2) { return; }
    if (beta < beta0) { beta = beta0; }

    // note the time for initializations
    long start_init = myclock();

    // insert f_lhs / h_lhs values into the node patches
    for (typename super::content_iterator np_it=super::contents_begin(), np_end=super::contents_end(); np_it != np_end; ++np_it)
    {
        NodePatchT& np = **np_it;
        size_t idx = np.get_idx();
        np.f = f_lhs[idx];
        np.h = h_lhs[idx];
    }

    // this holds the number of cubes which will require plane wave holders
    allocate_memory(); // re-indexes the tree and allocates memory

    // create level dependent globals
    for (unsigned short level = super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        boost::shared_ptr<const Level_Dependent_FMM_Globals<NTERMS,NLAMBS> > ptr_beta(new Level_Dependent_FMM_Globals<NTERMS,NLAMBS>(beta, level, super::universe_edge_length));
        boost::shared_ptr<const Level_Dependent_FMM_Globals<NTERMS,NLAMBS> > ptr_beta0(new Level_Dependent_FMM_Globals<NTERMS,NLAMBS>(beta0, level, super::universe_edge_length));
        fmm_level_globs_beta[level] = ptr_beta;
        fmm_level_globs_beta0[level] = ptr_beta0;
    }
    
    // put init time into profile counter
    timing_info.init += myclock() - start_init;

    long start_upward = myclock();
    calculate_multipole_expansions(beta, beta0);
    translate_multipole_upward(beta, beta0);
    timing_info.total_upward_pass += myclock() - start_upward;
    
    long start_downward = myclock();
    downward_pass(beta, beta0);
    timing_info.total_downward_pass += myclock() - start_downward;
    
    long start_eval = myclock();
    evaluate(beta, beta0, f_rhs, h_rhs);
    timing_info.fmm_evaluations += myclock() - start_eval;
    
    long clear_mem = myclock();
    
    // free the multipole memory
    mexp_arrays.clear();
    
    timing_info.allocate_memory += myclock() - clear_mem;
    
//     size_t num_octree_nodes = super::get_total_num_nodes();
//     size_t num_explicit_calcs = super::calc_neighbourhood_interacts();
//     size_t max_octree_depth = super::get_bottom_level() - super::get_top_level();
//     size_t num_fmm_evals = (max_octree_depth > 1) ? super::get_total_items() : 0;
//     
//     std::cout << "Number of octree nodes: " << num_octree_nodes << std::endl;
//     std::cout << "Number of explicit calcs: " << num_explicit_calcs << std::endl;
//     std::cout << "Number of fmm evals: " << num_fmm_evals << std::endl;
//     
//     
//     // loop over all nodes, don't bother with leaf nodes or sub-leaf nodes
//     // add up the interactions between children of these nodes
//     size_t ilist=0;
//     
//     // in the FMM near-field.
//     for (int level=super::get_top_level(); level <= super::get_bottom_level() && level < super::max_depth; ++level)
//     {
//         for (typename super::NodeList::const_iterator node_it=super::get_node_list(level).begin(), node_end=super::get_node_list(level).end(); 
//                 node_it != node_end; ++node_it)
//         {
//             const typename super::NodeT& node = *(node_it->second);
//             if (node.isLeaf() == false && node.isShadow() == false) 
//             {
//                 // get interaction list for this neighbour_id
//                 for (unsigned short ctr=0; ctr < ILIST_SIZE; ++ctr)
//                 {
//                     const INTERACTION& inter = interactions[ctr];
//                     if (inter.neighbour == 13) { continue; }
//                     try {
//                     
//                         OctreeIndexer idxer1 = node.get_idx().get_neighbour_idxer(inter.neighbour).get_child_idx(inter.targ);
//                         const typename super::NodeT& ch1 = super::get_node(idxer1);
//                         OctreeIndexer idxer2 = node.get_idx().get_child_idx(inter.src);
//                         const typename super::NodeT& ch2 = super::get_node(idxer2);
//                         ++ilist;
//                     }
//                     catch (...) { }
//                 }
//             }
//         }
//     }
//     std::cout <<"Number of interaction list moves: " << ilist << std::endl;
//     
}

// note allocate_memory must have been called before calling this function!
template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::calculate_multipole_expansions(double beta, double beta0)
{

    if (super::size() == 0) { return; }

    long start_create_mpole = myclock();
    size_t ctr=0;

    // workspace for legendre polynomials
    LegendreHolder<NTERMS> workspace_p;
    for (unsigned short level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        boost::shared_array< MultipoleHolder > multipoles = mexp_arrays[level];
        double scale = FMM_Globals<NTERMS>::get_scale(beta, super::universe_edge_length, level);
        double scale0 = FMM_Globals<NTERMS>::get_scale(beta0, super::universe_edge_length, level);

        for (typename super::const_childless_iterator childless_it=super::childless_begin(level), childless_end=super::childless_end(level);
            childless_it != childless_end;
            ++childless_it)
        {
            const NodeT& node = *(childless_it->second);
            ++ctr;

            // skip empty nodes -- childless iterator should do this automatically
            assert(node.hasChildren()==false);

            MultipoleHolder& mexp = multipoles[node.get_mp_idx()];
            const typename super::NodeT::ContentList& contents = node.get_contents();

            Vector node_centre = node.get_centre();

            yformmp<NodePatchT,12,0,8>(beta,
                        node_centre,
                        contents,
                        mexp,
                        scale,
                        workspace_p,
                        fmm_globs.scale_factors);

            yformmp<NodePatchT,12,8,12>(beta0,
                        node_centre,
                        contents,
                        mexp,
                        scale0,
                        workspace_p,
                        fmm_globs.scale_factors);
 
        }
    }
    timing_info.create_mpoles += myclock() - start_create_mpole;

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::translate_multipole_upward(double beta, double beta0)
{

    long start_xlate_mpoles = myclock();
    
    // upwardly translate the expansions to centres of parent cubes
    for (int from_level=super::get_bottom_level(); from_level > super::get_top_level(); --from_level)
    {

        // gonna need some precalcs to do the translations...
        //std::cout << "Pre-calculating MP shift coefficients for level " << from_level << std::endl;

        double octant_edge_length = super::universe_edge_length / pow(2, from_level);
        double r0 = sqrt(3.0) * octant_edge_length / 2.0;

        const DblMatrix3D& mp_translation_coeffs_beta = fmm_level_globs_beta[from_level]->multipole_shift_coefficients;
        const DblMatrix3D& mp_translation_coeffs_beta0 = fmm_level_globs_beta0[from_level]->multipole_shift_coefficients;

        // do the translations
        //std::cout << "Translating multipoles from level " << from_level << " to level " << from_level-1 << std::endl;
        boost::shared_array<MultipoleHolder>& parent_multipoles = mexp_arrays[from_level-1];

        // loop over cubes in level
        boost::shared_array<MultipoleHolder>& child_multipoles = mexp_arrays[from_level];

        // heap allocate some workspaces
        boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace1(new BaseMultipoleHolder<NTERMS>);
        boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace2(new BaseMultipoleHolder<NTERMS>);

        // loop over cubes in level
        for (typename super::NodeList::iterator it=super::get_node_list(from_level).begin(), end=super::get_node_list(from_level).end();
            it != end;
            ++it)
        {
            NodeT& node = *(it->second);

            // if this node has nothing in it then there's nothing to translate upwards except
            // a whole bunch of zeros.
            if (node.empty()) { continue; }

            // get parent node
#ifdef __DELETED__
            NodeT& parent_node = get_node(node.get_idx().get_parent_idx());
#else
            NodeT& parent_node = super::get_node(node.get_idx().get_parent_idx());
#endif

            // generate slice of multipoles matrix
            MultipoleHolder& child_mexp  = child_multipoles[node.get_mp_idx()];
            MultipoleHolder& parent_mexp = parent_multipoles[parent_node.get_mp_idx()];

            int id = node.get_id_within_parent();
            const DblMatrix3D& rotation_matrix =  (id > 4) ? fmm_globs.rdmsq3 : fmm_globs.rdsq3;

            ympshift<12,0,8>(id,
                        child_mexp,
                        parent_mexp,
                        *workspace1,
                        *workspace2,
                        mp_translation_coeffs_beta,
                        rotation_matrix);

            ympshift<12,8,12>(id,
                        child_mexp,
                        parent_mexp,
                        *workspace1,
                        *workspace2,
                        mp_translation_coeffs_beta0,
                        rotation_matrix);

        }
    }
    
    timing_info.pass_mpoles_upwards += myclock() - start_xlate_mpoles;
    
    return;
}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::inherit_downwards(double beta, double beta0)
{
    // nothing to do if not at least 2 levels
    if (super::get_bottom_level() < 2) { return; }
    
    long start_inherit_down = myclock();

    // loop over each level, processing interaction lists
    for (unsigned short level=super::get_top_level(); level < super::get_bottom_level(); ++level)
    {
        for (typename super::NodeList::iterator node_it=super::get_node_list(level).begin(), node_end=super::get_node_list(level).end();
            node_it != node_end;
            ++node_it)
        {
            NodeT& node = *(node_it->second);
            if (node.hasChildren() && node.isLeaf() == false)
            {
                pass_local_expansion_downwards(node, fmm_globs, *(fmm_level_globs_beta[level]), *(fmm_level_globs_beta0[level]));
            }
        }

    }

    timing_info.inherit_locals_downwards += myclock() - start_inherit_down;

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::downward_pass(double beta, double beta0)
{
    // nothing to do if not at least 2 levels
    if ( (super::get_bottom_level() - super::get_top_level()) < 2 || (super::size() == 0) ) {
        return;
    }

    interaction_lists(beta, beta0);
    
    long start_inheritance = myclock();
    inherit_downwards(beta, beta0);
    timing_info.inherit_locals_downwards += myclock() - start_inheritance;

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::interaction_lists(double beta, double beta0)
{

    // workspace (heap allocated)
    boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace1(new BaseMultipoleHolder<NTERMS>);
    boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace2(new BaseMultipoleHolder<NTERMS>);

    for (unsigned short level=super::get_top_level(); level < super::get_bottom_level(); ++level)
    {
        boost::shared_array<MultipoleHolder> child_level_multipoles = mexp_arrays[level+1];
        boost::shared_array<MultipoleHolder> child_level_locals = lexp_arrays[level+1];

        // planewave holders

        // big fat malloc for the plane wave expansion holders
        // could do this earlier since PlaneWaveHolder is made large enough
        // for maximum nexptotp
        long start_malloc = myclock();
        size_t plane_waves_needed = plane_wave_counters[level+1];
        const double mb_size = plane_waves_needed * sizeof(MultiHolder<6,PlaneWaveHolder>)  / double(1024*1024);
        //std::cout << "Allocating memory for plane waves... (" << mb_size << " Mbytes)" << std::endl;
        boost::shared_array< MultiHolder<6,PlaneWaveHolder> > child_level_planewaves(new MultiHolder<6,PlaneWaveHolder>[plane_waves_needed]);
        memset(child_level_planewaves.get(), 0, sizeof(MultiHolder<6,PlaneWaveHolder>)*plane_waves_needed);
        timing_info.allocate_memory += myclock() - start_malloc;

        long start_convert_to_pw = myclock();
        //std::cout << "Converting multipoles to plane waves for children of level " << level << std::endl;
        for (typename super::NodeList::iterator child_it=super::get_node_list(level+1).begin(), child_end=super::get_node_list(level+1).end();
            child_it != child_end;
            ++child_it)
        {
            NodeT& child_node = *(child_it->second);
            const OctreeIndexer& idxer = child_node.get_idx();

            // empty nodes will not have anything interesting to transmit to their neighbours
            if (child_node.empty()) { continue; }

            // TODO:: tidy this further
            size_t pw_idx = child_node.get_pw_idx();
            assert(pw_idx < plane_waves_needed);
            PlaneWaveHolder& up = child_level_planewaves[pw_idx](UP);
            PlaneWaveHolder& down = child_level_planewaves[pw_idx](DOWN);
            PlaneWaveHolder& north = child_level_planewaves[pw_idx](NORTH);
            PlaneWaveHolder& south = child_level_planewaves[pw_idx](SOUTH);
            PlaneWaveHolder& east = child_level_planewaves[pw_idx](EAST);
            PlaneWaveHolder& west = child_level_planewaves[pw_idx](WEST);

            MultipoleHolder& child_mpole = child_level_multipoles[child_node.get_mp_idx()];

            convert_mpole_to_six_planewaves<12,0,8>(fmm_globs, *(fmm_level_globs_beta[level]), child_mpole, up, down, north, south, east, west);
            convert_mpole_to_six_planewaves<12,8,12>(fmm_globs, *(fmm_level_globs_beta0[level]), child_mpole, up, down, north, south, east, west);
        }
        timing_info.convert_to_planewaves += myclock() - start_convert_to_pw;

        //std::cout << "Processing interaction lists for children of level " << level << std::endl;
        typedef InteractionListManager<typename super::NodeT, typename super::NodeList, NTERMS, NLAMBS, NWAVES, 12> InteractionList;
        boost::scoped_ptr<InteractionList> ilm(new InteractionList);
        boost::shared_array< MultiHolder<8,PlaneWaveHolder> > accum( new MultiHolder<8,PlaneWaveHolder>[6]);

        const typename super::NodeList& parent_nodes = super::get_node_list(level);
        const typename super::NodeList& child_nodes = super::get_node_list(level+1);

        // loop over cubes in level
        for (typename super::NodeList::const_iterator it=parent_nodes.begin(), end=parent_nodes.end();
            it != end;
            ++it)
        {
            const NodeT& node = *(it->second);

            // IMPORTANT: if we can evaluate the FMM here without further subdivision, then do not
            // carry out any translating on the lower level, as it is not necessary
            if (node.isLeaf() || node.isShadow()) { continue; }

            // time taken for interaction list
            long start_ilist = myclock();

            memset(accum.get(), 0, sizeof(MultiHolder<8,PlaneWaveHolder>)*6);
            ilm->set(node, child_level_planewaves.get(), parent_nodes, child_nodes);

            // kludge!  Because cannot use default arguments in template functions (until C++0x anyway)
            // we have to call these two explicit functions which handle each half of the 16-fold BEM/FMM 
            // multiple FMM collection.
            size_t num_ilist_xlations=0;
            num_ilist_xlations += (*ilm).do_interaction_list_A(accum.get(), *(fmm_level_globs_beta[level]));
            num_ilist_xlations += (*ilm).do_interaction_list_B(accum.get(), *(fmm_level_globs_beta0[level]));
            
            // profile time taken for ilist
            timing_info.interaction_list_xlations += myclock() - start_ilist;
            
            // stick the number of ilist-xlations in too for good measure
            timing_info.num_ilist_xlations += num_ilist_xlations;

            long start_convert_from_pw = myclock();
            for (int child_id=1; child_id <=8; ++child_id)
            {
                const OctreeIndexer child_idxer = node.get_idx().get_child_idx(child_id);
                try {
                    NodeT& child_node = super::get_node(child_idxer);
                    if (child_node.isShadow()) { continue; } // don't bother converting waves on shadows

                    size_t lc_idx = child_node.get_lc_idx();
                    MultipoleHolder& child_lexp = child_level_locals[lc_idx];

                    // TODO:: tidy this
                    // unravel plane wave expansion holders
                    PlaneWaveHolder& up = accum[UP](child_id-1);
                    PlaneWaveHolder& down = accum[DOWN](child_id-1);
                    PlaneWaveHolder& north = accum[NORTH](child_id-1);
                    PlaneWaveHolder& south = accum[SOUTH](child_id-1);
                    PlaneWaveHolder& east = accum[EAST](child_id-1);
                    PlaneWaveHolder& west = accum[WEST](child_id-1);

                    convert_and_add_accumulated_planewaves<12,0,8>(fmm_globs, *(fmm_level_globs_beta[level]), up, down, north, south, east, west, child_lexp);
                    convert_and_add_accumulated_planewaves<12,8,12>(fmm_globs, *(fmm_level_globs_beta0[level]), up, down, north, south, east, west, child_lexp);
                }
                catch (BadIndexer) {}

            }
            timing_info.convert_from_planewaves += myclock() - start_convert_from_pw;
        }
    }
    return;
}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::evaluate(double beta, double beta0, double f_rhs[], double h_rhs[])
{
    // workspaces
    boost::scoped_ptr<DblMatrix1D> pots_ptr;
    boost::scoped_ptr<DblMatrix2D> fields_ptr;
    boost::scoped_ptr<DblMatrix3D> f2_ptr;

    for (unsigned short level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        typename super::NodeList& nodes = super::get_node_list(level);
        for (typename super::NodeList::const_iterator it=nodes.begin(), end=nodes.end(); it != end; ++it)
        {
            // slightly cumbersome syntax- it's because NodeList is a
            // std::map not a vector/list
            const NodeT& node = *(it->second);

            // only iterate over leaves
            if (node.isLeaf() == false) { continue; }

            // loop over NodePatchT objects
            const typename NodeT::ContentList& contents = node.get_contents();
            
            #pragma omp parallel for private(pots_ptr, fields_ptr, f2_ptr)
            for (int np_idx=0; np_idx < contents.size(); ++np_idx)
            {
                
                // init workspaces (this works for multithreading too)
                if (pots_ptr.get() == NULL) {
                    pots_ptr.reset(new DblMatrix1D(Range(0,12)));
                }
                if (fields_ptr.get() == NULL) {
                    fields_ptr.reset(new DblMatrix2D(Range(0,3), Range(0,12)));
                }
                if (f2_ptr.get() == NULL) {
                    f2_ptr.reset(new DblMatrix3D(Range(0,3), Range(0,3), Range(0,12)));
                }
                
                const NodePatchT& np = *(contents[np_idx]);
                double& f_accum = f_rhs[np.get_idx()];
                double& h_accum = h_rhs[np.get_idx()];

                fmm::evaluate_local_expansion_at_xyz(beta,
                                                    beta0,
                                                    FMM_Globals<NTERMS>::get_scale(beta, super::universe_edge_length, level),
                                                    FMM_Globals<NTERMS>::get_scale(beta0, super::universe_edge_length, level),
                                                    np,
                                                    node.get_centre(),
                                                    get_lpole(node),
                                                    *pots_ptr, *fields_ptr, *f2_ptr);

                ::FMM_results_to_BEM(f_accum, h_accum, np.get_normal(), *pots_ptr, *fields_ptr, *f2_ptr);
            }

        }
    }

    return;

}

template<int NTERMS, int NLAMBS, int NWAVES, typename NodePatchT>
void FMM_BEM_Octree<NTERMS, NLAMBS, NWAVES,NodePatchT>::allocate_memory()
{

    long start_allocate_memory = myclock();

    // clear any existing multipole/local expansion data
    mexp_arrays.clear();
    lexp_arrays.clear();

    // traverse the tree and figure out how much of each type of holder we will require.

    // init. counters
    plane_wave_counters.clear();
    std::map<unsigned short, size_t> mp_ctrs;
    std::map<unsigned short, size_t> lc_ctrs;

    for (int i=super::get_top_level(); i <= super::get_bottom_level(); ++i)
    {
        plane_wave_counters[i] = 0;
        mp_ctrs[i] = 0;
        lc_ctrs[i] = 0;
    }
    
    //std::cout << "Traversing tree." << std::endl;
    for (int level=super::get_bottom_level(); level >= super::get_top_level(); --level)
    {
        size_t& pw_ctr = plane_wave_counters[level];
        size_t& mp_ctr = mp_ctrs[level];
        size_t& lc_ctr = lc_ctrs[level];

        // loop over cubes in level
        for (typename super::NodeList::iterator it=super::get_node_list(level).begin(), end=super::get_node_list(level).end();
            it != end;
            ++it)
        {
            NodeT& node = *(it->second);

            // set plane-wave holder index -- this is to conserve memory during local expansion downward pass
            // note that we omit some nodes (e.g. ones below the leaf nodes in an adaptive tree) because we
            // don't need the local expansions that far down the tree (because it is cheap enough to evaluate
            // the potential at a higher level in the tree... which is the whole point of the *adaptive* tree)
            if (node.isShadow() == false) {
                node.set_lc_idx(lc_ctr++);
            }

            // we need multipole expansions in all non-empty nodes of the tree
            if (!node.empty())
            {
                node.set_mp_idx(mp_ctr++);
                node.set_pw_idx(pw_ctr++);
            }
        }
    }

    // Allocate some memory

    // create Multipole Expansion holders (and Local Expansion and Plane Wave holders).
    //std::cout << "Allocating memory." << std::endl;

    for (int level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
    {
        boost::shared_array<MultipoleHolder> mpole(new MultipoleHolder[mp_ctrs[level]]);
        assert(mpole.get() != NULL);
        mexp_arrays[level] = mpole;
        memset(mpole.get(), 0, sizeof(MultipoleHolder)*mp_ctrs[level]);

        boost::shared_array<MultipoleHolder> local(new MultipoleHolder[lc_ctrs[level]]);
        assert(local.get() != NULL);
        lexp_arrays[level] = local;
        memset(local.get(), 0, sizeof(MultipoleHolder)*lc_ctrs[level]);
    }

    timing_info.allocate_memory += myclock() - start_allocate_memory;

    return;
}

} // end namespace

#endif /* FMM_BEM_HYBRID_H_ */
