#ifndef __FMM_WORKER_H__
#define __FMM_WORKER_H__

#include "prerequisites.h"
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include "fh_values_nodegroup.decl.h"
#include "fh_values_nodegroup.h"
#include "fmm_globals_nodegroup.decl.h"
#include "fmm_globals_nodegroup.h"
#include "parallel_fmm_octree.decl.h"
#include "parallel_fmm_octree.h"
#include "fmm_worker.decl.h"
#include "../fmm/fmm.h"
#include "main.decl.h"

using fmm::FMM_Globals;
using fmm::Level_Dependent_FMM_Globals;
using beepp::MultipoleHolderT;
using beepp::PlaneWaveHolderT;

using namespace beepp;
namespace beepp
{

//typedef ParallelFMMOctree<CharmNodePatch, CProxy_FMMWorker> ParallelTree;
    
template<class HType>
class FMM_Msg : public CMessage_FMM_Msg<HType>, public HType
{

public:

    FMM_Msg() : CMessage_FMM_Msg<HType>(), HType(), origin(OctreeIndexer(0,0,0,0)), num_waves_included(0) {}
    FMM_Msg(const CkArrayIndexOctreeIndexer& idx) : CMessage_FMM_Msg<HType>(), HType(), origin(idx),num_waves_included(0) {}
    inline operator HType& () { return static_cast<HType&>(*this); }
    inline operator const HType& () { return static_cast<const HType&>(*this); }

    CkArrayIndexOctreeIndexer origin;
    int num_waves_included;
};

template<int NTERMS, int NLAMBS, int NWAVES>
class FMMWorkerT : public CBase_FMMWorkerT<NTERMS, NLAMBS, NWAVES>
{
    typedef CProxy_FMMWorkerT<NTERMS, NLAMBS, NWAVES> proxy;
    typedef CBase_FMMWorkerT<NTERMS, NLAMBS, NWAVES> super;
    typedef ParallelFMMOctree<CharmNodePatch,proxy> ParallelTree;
    typedef typename ParallelTree::NodeT NodeT;
    typedef CkIndex_FMMWorkerT<NTERMS, NLAMBS, NWAVES> IdxType;
    typedef CProxySection_FMMWorkerT<NTERMS, NLAMBS, NWAVES> SectionProxyType;
    typedef MultipoleHolderT<NTERMS> MultipoleHolder;
    typedef PlaneWaveHolderT<NLAMBS,NWAVES> PlaneWaveHolder;
    typedef FMM_Msg<MultipoleHolder> MultipoleMsg;
    typedef FMM_Msg< MultiHolder<6,PlaneWaveHolderT<NLAMBS,NWAVES> > > PlaneWaveMsg;
    
private:

    CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS> FMM_Globals_Proxy;
    CProxy_ParallelTree FMM_Tree_Proxy;
    CProxy_FH_Values_NodeGroup FH_Proxy;
    boost::shared_ptr< MultipoleHolder > mpole_ptr;
    boost::shared_ptr< MultipoleHolder > lpole_ptr;
    boost::shared_array< PlaneWaveHolder > incoming_waves;  // 1 holder per direction: up/down/north/south/east/west
    CkCallback stashed_cb; // this gets passed in at kickoff, and tells where to notify of completion

public:

    /// Constructors ///
    FMMWorkerT(CkMigrateMessage *msg) {}// std::cout << super::thisIndex << ": FMMWorkerT cannot migrate!" << std::endl; throw std::exception();}

    FMMWorkerT(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS> fmm_globals_proxy,
                CProxy_ParallelTree fmm_tree_proxy,
                CProxy_FH_Values_NodeGroup fh_proxy,
                double edge_length) :
        num_children(-1),
        expected_incoming_waves(-1)
    {
        // stash the proxies for the various things we're gonna need
        FMM_Globals_Proxy = fmm_globals_proxy;
        FMM_Tree_Proxy = fmm_tree_proxy;
        FH_Proxy = fh_proxy;

        universe_edge_length = edge_length;
        level = super::thisIndex.get_level();

        make_interaction_list();
    }

    void init_memory();

    void solve(double beta_,
               double beta0_,
               CkCallback cb)
    {
        // stash the callback -- send it when have completed work (see evaluate function)
        stashed_cb = cb;

        beta = beta_;
        beta0 = beta0_;
        
        init_memory();
        super::contribute(CkCallback(IdxType::form_multipoles(), super::thisProxy));
    }
    
    void form_multipoles()
    {
        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);

        // if it's a leaf node, and is either level 0 or 1, then
        // there is nothing to be done- go straight to 'evaluation'
        // which will notify completion
        if (level == 0 || level == 1) {
            if (node.isLeaf()) {
                super::thisProxy[super::thisIndex].evaluate();
            }
            return;
        }

        // only bottom level nodes need do anything to begin with...
        if (node.hasChildren() == true) { return; }

        // we should have delete empty nodes as an optimization before this point
        assert(node.empty() == false);

        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const double* fvals = FH_Proxy.ckLocalBranch()->fvals();
        const double* hvals = FH_Proxy.ckLocalBranch()->hvals();
        // insert f_lhs / h_lhs values into the node patches
        std::vector<CharmNodePatch*> contents = node.get_contents();

        for (std::vector<CharmNodePatch*>::const_iterator np_it=contents.begin(), np_end=contents.end(); np_it != np_end; ++np_it)
        {
            CharmNodePatch& np = **np_it;
            size_t idx = np.get_idx();
    #ifdef USE_JUFFER_PEAKS
            np.f += fvals[idx];
            np.h += hvals[idx];
    #else
            np.f = fvals[idx];
            np.h = hvals[idx];
    #endif

        }

        // workspace for legendre polynomials
        LegendreHolder<NTERMS> workspace_p;
        double scale = FMM_Globals<NTERMS>::get_scale(beta, universe_edge_length, level);
        double scale0 = FMM_Globals<NTERMS>::get_scale(beta0, universe_edge_length, level);

        fmm::yformmp<CharmNodePatch,12,0,8>(beta,
                    node.get_centre(),
                    node.get_contents(),
                    *mpole_ptr,
                    scale,
                    workspace_p,
                    fmm_globs.scale_factors);

        fmm::yformmp<CharmNodePatch,12,8,12>(beta0,
                    node.get_centre(),
                    node.get_contents(),
                    *mpole_ptr,
                    scale0,
                    workspace_p,
                    fmm_globs.scale_factors);

        //std::cout << super::thisIndex << "Created mpole: " << *mpole_ptr << std::endl;

    #ifdef USE_JUFFER_PEAKS
        for (std::vector<CharmNodePatch*>::const_iterator np_it=contents.begin(), np_end=contents.end(); np_it != np_end; ++np_it)
        {
            CharmNodePatch& np = **np_it;
            size_t idx = np.get_idx();
            np.f -= fvals[idx];
            np.h -= hvals[idx];
        }
    #endif

        super::thisProxy[super::thisIndex].pass_multipole_upwards();
        super::thisProxy[super::thisIndex].make_plane_waves(0);

        return;
    }

    void make_plane_waves(int=0);
    void pass_multipole_upwards();

    void process_interaction(const MultiHolder<6,PlaneWaveHolder>& waves, const OctreeIndexer& target, const INTERACTION& inter)
    {

        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level-1, universe_edge_length);
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta0 = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta0_level(beta0, level-1, universe_edge_length);

        ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        PlaneWaveMsg* msg;
        tree.get_waves_lock(target, msg);
        ++(msg->num_waves_included);
        MultiHolder<6,PlaneWaveHolder>& collected = *msg;
       
        // note that there's a coordinate rotation here such that the plane wave is always translated in the z direction
        // (see the make_plane_waves function -- you'll see a rotation of the multipole so that the direction of the plane wave
        // is always in the z direction, prior to the formation of the plane wave expansion from the multipole.)
        short shift_x=0, shift_y=0, shift_z=0;
        DIRECTION dir = inter.dir;
        switch (dir)
        {
        case UP:
        case DOWN:
            shift_x = inter.dx;
            shift_y = inter.dy;
            shift_z = inter.dz;
            break;
        case NORTH:
        case SOUTH:
            shift_x = inter.dz;
            shift_y = inter.dx;
            shift_z = inter.dy;
            break;
        case EAST:
        case WEST:
            shift_x = -inter.dz;
            shift_y = inter.dy;
            shift_z = inter.dx;
            break;
        default:
            // should not get here -- throw exception?
            CkPrintf("Unexpected direction passed to receive_plane_wave_from_interaction_list()\n");
            assert(false);
            break;
        }

        switch (dir)
        {
        case UP:
        case NORTH:
        case EAST:
            fmm::translate_up_wave<12,0,8>(level_globs_beta, waves[dir], collected[dir], shift_x, shift_y, shift_z);
            fmm::translate_up_wave<12,8,12>(level_globs_beta0, waves[dir], collected[dir], shift_x, shift_y, shift_z);
            break;
        case DOWN:
        case SOUTH:
        case WEST:
            fmm::translate_down_wave<12,0,8>(level_globs_beta, waves[dir], collected[dir], shift_x, shift_y, shift_z);
            fmm::translate_down_wave<12,8,12>(level_globs_beta0, waves[dir], collected[dir], shift_x, shift_y, shift_z);
            break;
        default:
            // should not get here -- throw exception?
            CkPrintf("Unexpected direction passed to receive_plane_wave_from_interaction_list()\n");
            assert(false);
            break;
        }

        
        // release the lock on the local plane waves collection
        tree.release_waves_lock(target);

    }
    
    void request_collected_waves()
    {
        // deal with the case of not requiring any incoming waves
        FMM_Tree_Proxy.request_data(super::thisProxy);            
    }
    
    void receive_incoming_wave(CkMessage* msg_in)
    {
        
        assert(expected_incoming_waves > 0 && done_gathered_lpoles == false);
        if (expected_incoming_waves > 0 && incoming_waves.get() == NULL)
        {
            incoming_waves = boost::shared_array< PlaneWaveHolder >(new PlaneWaveHolder[6]);
            for (int ii=0; ii < 6; ++ii) { incoming_waves[ii].reset(); }
        }

        // the incoming planewave holder contains all 6 directions from the source node
        PlaneWaveMsg* msg_recast = reinterpret_cast<PlaneWaveMsg*>(msg_in);
        const MultiHolder<6,PlaneWaveHolder>& multi_pw_holder = *msg_recast;
        for (int ii=0; ii < 6; ++ii) { 
            incoming_waves[ii] += multi_pw_holder[ii]; 
        }
        
        incoming_wave_ctr += msg_recast->num_waves_included;
        check_waves_complete(0);
    }    
    
    void inherit_lpole_from_parent(MultipoleMsg *lpole_msg);
    void receive_multipole_contribution_from_child(MultipoleMsg *mpole_msg);
    void check_waves_complete(int);

    void evaluate();
    void debug_chk();

    void pup(PUP::er &p);

private:

    bool done_evals;
    void make_interaction_list();
    
    void check_locals_complete();

    unsigned short level;
    double universe_edge_length;
    double beta;
    double beta0;

    int fmm_upward_pass_contribution_counter;
    int num_children;
    bool done_gathered_lpoles;
    bool done_inherited_lpole;
    bool passed_upwards;
    bool made_plane_waves;

    std::vector< std::pair<CkArrayIndexOctreeIndexer, INTERACTION> > my_outgoing_interaction_list;
    int expected_incoming_waves;
    int incoming_wave_ctr;

};

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::debug_chk() {
    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    const NodeT& node = tree.get_node(super::thisIndex);
    if (node.isLeaf() == true && done_evals == false)
    {
        std::cout << "DEBUG!! " << super::thisIndex << " " << node.hasChildren() << " " << node.isShadow() << " " << node.isLeaf() << std::endl;
        assert(false);
    }
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::init_memory()
{
    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    const NodeT& node = tree.get_node(super::thisIndex);

    // allocate memory
    mpole_ptr.reset(new MultipoleHolder);
    mpole_ptr->reset();
    if (node.isShadow() == false) {
        lpole_ptr.reset(new MultipoleHolder);
        lpole_ptr->reset();
    }
    incoming_waves.reset();

    level = super::thisIndex.get_level();
    incoming_wave_ctr = 0;

    fmm_upward_pass_contribution_counter = 0;
    done_gathered_lpoles = false;
    done_inherited_lpole = false;
    passed_upwards = false;
    made_plane_waves = false;

    // root level does not inherit from parent (it has no parent!)
    // also root level and level 1 do not have interaction lists, so
    // no need to gather anything either; since level 0 has no
    // local multipole, there is nothing for level 1 to inherit;
    // ditto level 2
    if (level == 0 || level == 1) {
        done_inherited_lpole = true;
        done_gathered_lpoles = true;
    }
    if (level == 2)
    {
        done_inherited_lpole = true;
    }

    return;

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::evaluate()
{
    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    typedef NodeT NodeT;
    const NodeT& node = tree.get_node(super::thisIndex);

    // only evaluation nodes should actually get here
    assert(node.isLeaf() == true);

    if (level == 0 || level == 1) {
        done_evals = true;
        stashed_cb.send();
        return;
    }

    //std::cout << "FMM Evaluating " << super::thisIndex << std::endl;

    unsigned short level = node.get_level();
    const double scale = FMM_Globals<NTERMS>::get_scale(beta, tree.get_edge_length(), level);
    const double scale0 = FMM_Globals<NTERMS>::get_scale(beta0, tree.get_edge_length(), level);
    const Vector& node_centre =	node.get_centre();

    //const MultipoleHolder& local_expansions = tree.get_lpole(node);
    const MultipoleHolder& local_expansions = *lpole_ptr;

    DblMatrix1D pots(Range(0,12));
    DblMatrix2D fields(Range(0,3), Range(0,12));
    DblMatrix3D grad_fields(Range(0,3), Range(0,3), Range(0,12));

    FH_Values_NodeGroup& fh_vals_group = *(FH_Proxy.ckLocalBranch());
    unsigned int total_num_patches = fh_vals_group.get_num_patches();
    boost::scoped_array<double> results(new double[total_num_patches*2]);
    memset(results.get(), 0, sizeof(double)*total_num_patches*2);

    // loop over CharmNodePatch objects
    const typename NodeT::ContentList& contents = node.get_contents();
    for (typename NodeT::ContentList::const_iterator np_it=contents.begin(), np_end=contents.end();
        np_it != np_end;
        ++np_it)
    {
        const CharmNodePatch& np = **np_it;
        const Vector& pt = np.get_node();
        double& f_result=results[np.get_idx()];
        double& h_result=results[np.get_idx() + total_num_patches];

        fmm::evaluate_local_expansion_at_xyz(beta,
                                            beta0,
                                            scale,
                                            scale0,
                                            np.get_node(),
                                            node_centre,
                                            local_expansions,
                                            pots, fields, grad_fields);

        ::FMM_results_to_BEM(f_result, h_result, np.get_normal(), pots, fields, grad_fields);

    }

    // add them into the node group holder
    fh_vals_group.add_results(results.get());

    // can now safely delete the memory allocated for local expansions
    lpole_ptr.reset();

    done_evals = true;
    stashed_cb.send();
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::pup(PUP::er &p) {
    super::pup(p);

    p | level;
    p | universe_edge_length;
    p | beta;
    p | beta0;

    p | fmm_upward_pass_contribution_counter;
    p | num_children;
    p | done_gathered_lpoles;
    p | done_inherited_lpole;
    p | passed_upwards;
    p | made_plane_waves;

    unsigned int num_interactions;
    if (p.isPacking())
    {
        num_interactions = static_cast<unsigned int>(my_outgoing_interaction_list.size());
    }
    p | num_interactions;

    p | my_outgoing_interaction_list;

    p | expected_incoming_waves;
    p | incoming_wave_ctr;

    bool mp_set,lp_set,pw_set;
    if (p.isPacking())
    {
        mp_set = mpole_ptr.get() != NULL;
        lp_set = lpole_ptr.get() != NULL;
        pw_set = incoming_waves.get() != NULL;
    }

    p | mp_set;
    p | lp_set;
    p | pw_set;

    if (p.isUnpacking())
    {
        if (mp_set) { mpole_ptr.reset(new MultipoleHolder); }
        if (lp_set) { lpole_ptr.reset(new MultipoleHolder); }
        if (pw_set) { incoming_waves.reset(new PlaneWaveHolder[6]); }
    }

    if (mpole_ptr.get()) { p | *mpole_ptr; }
    if (lpole_ptr.get()) { p | *lpole_ptr; }
    if (incoming_waves.get())
    {
        for (int ii=0; ii < 6; ++ii)
        {
            p | incoming_waves[ii];
        }
    }
    
    p | stashed_cb;

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::make_interaction_list()
{
    my_outgoing_interaction_list.clear();
    expected_incoming_waves = 0;

    // no interaction list for top level
    if (level > 0) { 
        
        const ParallelTree* tree_ptr = FMM_Tree_Proxy.ckLocalBranch();
        assert(tree_ptr != NULL);
        const ParallelTree& tree = *(tree_ptr);
        const NodeT& node = tree.get_node(super::thisIndex);

        unsigned short my_id = node.get_id_within_parent();
        const OctreeIndexer parent_idxer = node.get_idx().get_parent_idx();

        for (int ii=0; ii < ILIST_SIZE; ++ii)
        {
            const INTERACTION& inter = interactions[ii];
            if (static_cast<unsigned short>(inter.targ) != my_id) { continue; }
            try {
                OctreeIndexer nb_idxer = parent_idxer.get_neighbour_idxer(inter.neighbour);
                OctreeIndexer remote_idxer = nb_idxer.get_child_idx(inter.src);
                const NodeT& remote_node = tree.get_node(remote_idxer);

                // if we get here without raising exceptions then the remote node must
                // exist and be non-empty.
                // if it is also not under a leaf, then it is in the interaction list
                if (remote_node.isShadow() == false) {

                    std::pair<CkArrayIndexOctreeIndexer,INTERACTION> p=std::pair<CkArrayIndexOctreeIndexer,INTERACTION>(CkArrayIndexOctreeIndexer(remote_idxer), invert(inter));
                    my_outgoing_interaction_list.push_back(p);
                }

                if (node.isShadow() == false) {
                    ++expected_incoming_waves;
                }
            }
            catch (BadIndexer) {}

        }
    }
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::check_waves_complete(int)
{
    //std::cout << super::thisIndex << "expecting: " << expected_incoming_waves << " got: " << incoming_wave_ctr << std::endl;
    assert(incoming_wave_ctr <= expected_incoming_waves);
    if (expected_incoming_waves == incoming_wave_ctr)
    {
        // plane wave holders should have been allocated, or else we're in trouble
        assert(incoming_waves.get() != NULL);

        // Get FMM globs -- some are level dependent
        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level-1, universe_edge_length);
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta0 = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta0_level(beta0, level-1, universe_edge_length);

        // convert plane waves to local expansion
        const PlaneWaveHolder& up =    incoming_waves[UP];
        const PlaneWaveHolder& down =  incoming_waves[DOWN];
        const PlaneWaveHolder& north = incoming_waves[NORTH];
        const PlaneWaveHolder& south = incoming_waves[SOUTH];
        const PlaneWaveHolder& east =  incoming_waves[EAST];
        const PlaneWaveHolder& west =  incoming_waves[WEST];

        // Convert mpole expansion into PlaneWave expansions (up/down/north/south/east/west).
        assert(lpole_ptr.get() != NULL);
        fmm::convert_and_add_accumulated_planewaves<12,0,8>(fmm_globs, level_globs_beta, up, down, north, south, east, west, *lpole_ptr);
        fmm::convert_and_add_accumulated_planewaves<12,8,12>(fmm_globs, level_globs_beta0, up, down, north, south, east, west, *lpole_ptr);

        incoming_waves.reset(); // no need to retain this memory
        incoming_wave_ctr = 0; // reset the counter whilst we're here

        done_gathered_lpoles = true;
        check_locals_complete();

    }

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::check_locals_complete()
{
    //std::cout << super::thisIndex <<" check_locals_complete (gather) (inherit):" << done_gathered_lpoles << " " << done_inherited_lpole << std::endl;
    if (done_gathered_lpoles && done_inherited_lpole)
    {
        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);
        if (node.isShadow() == true) {
            return;
        }

        if (node.isLeaf())
        {
            // if this node is an evaluation node, then no need to
            // pass the local expansions any further.  Just get on
            // with evaluating the stuff in this node.
            super::thisProxy[super::thisIndex].evaluate();
        }
        else
        {
            const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level, universe_edge_length);
            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta0 = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta0_level(beta0, level, universe_edge_length);

            boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace1(new BaseMultipoleHolder<NTERMS>);
            boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace2(new BaseMultipoleHolder<NTERMS>);

            for (int child_octant_id=1; child_octant_id <= 8; ++child_octant_id)
            {
                OctreeIndexer child_idxer = node.get_idx().get_child_idx(child_octant_id);
                try {
                    const NodeT& test = tree.get_node(child_idxer);
                }
                catch (BadIndexer) { continue; }

                MultipoleMsg* local_out = new MultipoleMsg();
                workspace1->reset();
                workspace2->reset();
                const DblMatrix3D& rotation_matrix = (child_octant_id > 4) ? fmm_globs.rdsq3 : fmm_globs.rdmsq3;

                fmm::ylcshift<12,0,8>(child_octant_id,
                                *lpole_ptr,
                                static_cast<MultipoleHolder&>(*local_out),
                                *workspace1,
                                *workspace2,
                                level_globs_beta.local_shift_coefficients,
                                rotation_matrix);
                fmm::ylcshift<12,8,12>(child_octant_id,
                                *lpole_ptr,
                                static_cast<MultipoleHolder&>(*local_out),
                                *workspace1,
                                *workspace2,
                                level_globs_beta0.local_shift_coefficients,
                                rotation_matrix);
                super::thisProxy[child_idxer].inherit_lpole_from_parent(local_out);
            }
        }
    }

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::inherit_lpole_from_parent(MultipoleMsg *lpole_msg)
{
    *lpole_ptr += *lpole_msg;
    delete lpole_msg;

    done_inherited_lpole = true;

    check_locals_complete();
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::make_plane_waves(int)
{

    // if there is nowhere to send stuff, don't bother creating plane wave expansions
    if (my_outgoing_interaction_list.size() != 0) 
    { 

        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);
        OctreeIndexer parent_idxer = super::thisIndex.get_parent_idx();

        // Convert mpole expansion into PlaneWave expansions (up/down/north/south/east/west).
        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level-1, universe_edge_length);
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta0 = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta0_level(beta0, level-1, universe_edge_length);

        boost::scoped_ptr< MultiHolder<6, PlaneWaveHolder> > pw_ptr(new MultiHolder<6, PlaneWaveHolder>); // cast the message as the holder type
        MultiHolder<6, PlaneWaveHolder>& pw_out = *pw_ptr;
        PlaneWaveHolder& up =    pw_out[UP];
        PlaneWaveHolder& down =  pw_out[DOWN];
        PlaneWaveHolder& north = pw_out[NORTH];
        PlaneWaveHolder& south = pw_out[SOUTH];
        PlaneWaveHolder& east =  pw_out[EAST];
        PlaneWaveHolder& west =  pw_out[WEST];

        fmm::convert_mpole_to_six_planewaves<12,0,8>(fmm_globs, level_globs_beta, *mpole_ptr, up, down, north, south, east, west);
        fmm::convert_mpole_to_six_planewaves<12,8,12>(fmm_globs, level_globs_beta0, *mpole_ptr, up, down, north, south, east, west);
        
        unsigned int interactees = my_outgoing_interaction_list.size();
        for (unsigned int ii=0; ii < interactees; ++ii)
        {
            assert(ii < (6*6*6 - 27));
            process_interaction(pw_out, my_outgoing_interaction_list[ii].first, my_outgoing_interaction_list[ii].second);
        }
        made_plane_waves = true;
        if (passed_upwards && made_plane_waves) { mpole_ptr.reset(); }
        
    }
    
    // sync point
    CkCallback cb(this, IdxType::request_collected_waves());
    super::contribute(cb);
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::pass_multipole_upwards()
{
    //std::cout << "Passing upwards: " << super::thisIndex << (*mpole_ptr) << std::endl;

    if (level == 0) {
        // All done-- multipoles have reached the top of the tree.
        //std::cout << "Root mpole: " << (*mpole_ptr) << std::endl;
        return;
    }

    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    const NodeT& node = tree.get_node(super::thisIndex);

    // do the translations
    //std::cout << "Translating multipoles from level " << from_level << " to level " << from_level-1 << std::endl;

    // heap allocate some workspaces
    boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace1(new BaseMultipoleHolder<NTERMS>);
    boost::scoped_ptr<BaseMultipoleHolder<NTERMS> > workspace2(new BaseMultipoleHolder<NTERMS>);

    MultipoleMsg *msg = new MultipoleMsg;
    MultipoleHolder& parent_mexp = *msg;
    parent_mexp.reset();

    unsigned short id = node.get_id_within_parent();
    const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
    const DblMatrix3D& mp_translation_coeffs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level, universe_edge_length).multipole_shift_coefficients;
    const DblMatrix3D& mp_translation_coeffs_beta0 = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta0_level(beta0, level, universe_edge_length).multipole_shift_coefficients;
    const DblMatrix3D& rotation_matrix =  (id > 4) ? fmm_globs.rdmsq3 : fmm_globs.rdsq3;

    fmm::ympshift<12,0,8>(id,
                        *mpole_ptr,
                        parent_mexp,
                        *workspace1,
                        *workspace2,
                        mp_translation_coeffs_beta,
                        rotation_matrix);

    fmm::ympshift<12,8,12>(id,
                        *mpole_ptr,
                        parent_mexp,
                        *workspace1,
                        *workspace2,
                        mp_translation_coeffs_beta0,
                        rotation_matrix);

    // send up
    //std::cout << node.get_idx().get_parent_idx() << " Should get:" << *parent_mexp << std::endl;
    super::thisProxy[node.get_idx().get_parent_idx()].receive_multipole_contribution_from_child(msg);

    passed_upwards = true;
    if (passed_upwards && made_plane_waves) { mpole_ptr.reset(); }

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void FMMWorkerT<NTERMS,NLAMBS,NWAVES>::receive_multipole_contribution_from_child(MultipoleMsg *mpole_msg)
{

    // if we haven't counted how many children we have yet, do so now.
    if (num_children == -1)
    {
        num_children = 0;
        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);

        for (int child_id=1; child_id <= 8; ++child_id)
        {
            OctreeIndexer idxer = super::thisIndex.get_child_idx(child_id);
            try {
                const NodeT& test = tree.get_node(idxer);
                ++num_children;
            }
            catch (BadIndexer) {}
        }
    }

    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    const NodeT& node = tree.get_node(super::thisIndex);
    assert(node.hasChildren() == true);
    assert(num_children > 0);

    // add mpole_in to internal multipole holder.
    *mpole_ptr += *mpole_msg;
    delete mpole_msg;

    ++fmm_upward_pass_contribution_counter;
    assert(fmm_upward_pass_contribution_counter <= num_children);
    //std::cout << super::thisIndex << " upward pass counter: " << fmm_upward_pass_contribution_counter  << "(expecting: " << num_children << ")" << std::endl;
    if (fmm_upward_pass_contribution_counter == num_children)
    {
        pass_multipole_upwards();

        // whilst that is going on, we could do something useful such as convert mpoles to plane waves
        super::thisProxy[super::thisIndex].make_plane_waves(0);
    }

    return;
}

} // end namespace

#define CK_TEMPLATES_ONLY
#include "fmm_worker.def.h"
#undef CK_TEMPLATES_ONLY

#endif //__FMM_WORKER_H__
