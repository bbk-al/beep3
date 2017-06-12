#ifndef __VANILLA_FMM_WORKER_H__
#define __VANILLA_FMM_WORKER_H__

#include "vanilla_fmm_prerequisites.h"
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include "fmm_globals_nodegroup.decl.h"
#include "fmm_globals_nodegroup.h"
#include "parallel_fmm_octree.decl.h"
#include "parallel_fmm_octree.h"
#include "../fmm/fmm.h"
#include "main.decl.h"

// Charm++ header for multicast support on section proxies
//#include "ckmulticast.h"
//extern CProxy_CkMulticastMgr MulticastMgrProxy;

#include "vanilla_fmm_worker.decl.h"
    
using fmm::FMM_Globals;
using fmm::Level_Dependent_FMM_Globals;
using namespace vanilla_fmm;

namespace vanilla_fmm
{

template<class HType>
class FMM_Msg : public CMessage_FMM_Msg<HType>, public HType
{

public:

    FMM_Msg() : CMessage_FMM_Msg<HType>(), HType(), origin(OctreeIndexer(0,0,0,0)),num_waves_included(0) {}
    FMM_Msg(const CkArrayIndexOctreeIndexer& idx) : CMessage_FMM_Msg<HType>(), HType(), origin(idx), num_waves_included(0) {}
    inline operator HType& () { return static_cast<HType&>(*this); }
    inline operator const HType& () { return static_cast<const HType&>(*this); }

    CkArrayIndexOctreeIndexer origin;
    int num_waves_included;
    
};

template<int NTERMS, int NLAMBS, int NWAVES>
class VanillaFMMWorkerT : public CBase_VanillaFMMWorkerT<NTERMS, NLAMBS, NWAVES>
{
    typedef CProxy_VanillaFMMWorkerT<NTERMS, NLAMBS, NWAVES> proxy;
    typedef CBase_VanillaFMMWorkerT<NTERMS, NLAMBS, NWAVES> super;
    typedef ParallelFMMOctree<Charge,proxy> ParallelTree;
    typedef typename ParallelTree::NodeT NodeT;
    typedef CkIndex_VanillaFMMWorkerT<NTERMS, NLAMBS, NWAVES> IdxType;
    typedef CProxySection_VanillaFMMWorkerT<NTERMS, NLAMBS, NWAVES> SectionProxyType;
    typedef MultipoleHolderT<NTERMS> MultipoleHolder;
    typedef PlaneWaveHolderT<NLAMBS,NWAVES> PlaneWaveHolder;
    typedef FMM_Msg<MultipoleHolder> MultipoleMsg;
    typedef FMM_Msg< MultiHolder<6,PlaneWaveHolderT<NLAMBS,NWAVES> > > PlaneWaveMsg;

private:

    CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS> FMM_Globals_Proxy;
    CProxy_ParallelTree FMM_Tree_Proxy;
    boost::shared_ptr< MultipoleHolder > mpole_ptr;
    boost::shared_ptr< MultipoleHolder > lpole_ptr;
    boost::shared_array< PlaneWaveHolder > incoming_waves;  // 1 holder per direction: up/down/north/south/east/west
    CkCallback stashed_cb; // this gets passed in at kickoff, and tells where to notify of completion

public:

    /// Constructors ///
/*    VanillaFMMWorkerT() :
        ready_to_evaluate(false),
        num_children(-1),
        expected_incoming_waves(-1)
    {}*/
    
    VanillaFMMWorkerT(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS> fmm_globals_proxy,
                      CProxy_ParallelTree fmm_tree_proxy,
                      double edge_length) :
        ready_to_evaluate(false),
        num_children(-1),
        expected_incoming_waves(-1)
        
    {
        // stash the proxies for the various things we're gonna need
        FMM_Globals_Proxy = fmm_globals_proxy;
        FMM_Tree_Proxy = fmm_tree_proxy;

        universe_edge_length = edge_length;
        level = super::thisIndex.get_level();
        
        init_memory();

    }


    VanillaFMMWorkerT(CkMigrateMessage *msg) {} //{ std::cout << super::thisIndex << ": FMMWorkerT cannot migrate!" << std::endl; throw std::exception(); }

    void init_memory();

    void solve(double beta_, CkCallback cb)
    {
        
        beta = beta_;
        if (beta < 1e-10) { beta = 1e-10; } // beta cannot be exactly zero
        stashed_cb = cb;

        
        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);

        // if it's a leaf node, and is either level 0 or 1, then
        // there is nothing to be done- go straight to 'evaluation'
        // which will notify completion
        if (level == 0 || level == 1) {
            ready_to_evaluate = true;
            
            // contribute to the reduction on the stashed callback- that'll return
            // to the mainchare when we're totally done
            if (node.isLeaf()) {
                stashed_cb.send();
            }
            
            // sync point for completed plane waves
            CkCallback cb(this, IdxType::request_collected_waves());
            super::contribute(cb);
            
            return;
        }
        
        if (expected_incoming_waves == 0) {
            done_gathered_lpoles = true;
            check_locals_complete();
        }

        // only bottom level nodes need do anything to begin with...
        if (node.hasChildren() == true) { return; }

        // we should have delete empty nodes as an optimization before this point
        assert(node.empty() == false);

        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();

        // workspace for legendre polynomials
        LegendreHolder<NTERMS> workspace_p;
        double scale = FMM_Globals<NTERMS>::get_scale(beta, universe_edge_length, level);

        fmm::yformmp<Charge,1,0,1>(beta,
                    node.get_centre(),
                    node.get_contents(),
                    *mpole_ptr,
                    scale,
                    workspace_p,
                    fmm_globs.scale_factors);
                    
        //std::cout << super::thisIndex << "Created mpole: " << *mpole_ptr << std::endl;
        super::thisProxy[super::thisIndex].make_plane_waves(0);
        pass_multipole_upwards();

        return;
    }

    void make_plane_waves(int dummy);
    void pass_multipole_upwards();

    void process_interaction(const MultiHolder<6,PlaneWaveHolder>& waves, const OctreeIndexer& target, const INTERACTION& inter)
    {

        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level-1, universe_edge_length);

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
            fmm::translate_up_wave<1,0,1>(level_globs_beta, waves[dir], collected[dir], shift_x, shift_y, shift_z);
            break;
        case DOWN:
        case SOUTH:
        case WEST:
            fmm::translate_down_wave<1,0,1>(level_globs_beta, waves[dir], collected[dir], shift_x, shift_y, shift_z);
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
        
        assert(expected_incoming_waves > 0);
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

    void evaluate(Eval_Message* msg);
    void asynchronous_check(int dummy); // opencl only
    
    void pup(PUP::er &p);

    void check_waves_complete(int dummy);
private:

    bool ready_to_evaluate;
    void make_interaction_list();
    inline void check_locals_complete();

    unsigned short level;
    double universe_edge_length;
    double beta;

    int fmm_upward_pass_contribution_counter;
    int num_children;
    bool done_gathered_lpoles;
    bool done_inherited_lpole;
    bool passed_upwards;
    bool made_plane_waves;

    std::vector< std::pair<CkArrayIndexOctreeIndexer, INTERACTION> > my_outgoing_interaction_list;
    int expected_incoming_waves;
    int incoming_wave_ctr;
    
    // a list of the callbacks and OpenCL tracking ids which are pending
    // for the evaluate function (because the OpenCL handler is deliberately
    // asynchronous, don't want to block the thread waiting for it)
    typedef std::pair<Eval_Message*, std::vector<long> > Msg_Tracking_Pair;
    std::vector< Msg_Tracking_Pair > work_pending;

};

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::init_memory()
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
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::pup(PUP::er &p) 
{
    super::pup(p);

    if (p.isPacking())
    {
        if(work_pending.size() != 0)
        {
            std::cerr << "ERROR: cannot migrate a worker if it's waiting for OpenCL results." << std::endl;
            assert(false);
        }
    }
    
    p | level;
    p | universe_edge_length;
    p | beta;

    p | fmm_upward_pass_contribution_counter;
    p | num_children;
    p | ready_to_evaluate;
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

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::make_interaction_list()
{
    my_outgoing_interaction_list.clear();
    expected_incoming_waves = 0;

    // no interaction list for top level
    if (level > 0) { 
        
        const ParallelFMMOctree<Charge,proxy>* tree_ptr = FMM_Tree_Proxy.ckLocalBranch();
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
void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::check_waves_complete(int dummy)
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

        // convert plane waves to local expansion
        const PlaneWaveHolder& up =    incoming_waves[UP];
        const PlaneWaveHolder& down =  incoming_waves[DOWN];
        const PlaneWaveHolder& north = incoming_waves[NORTH];
        const PlaneWaveHolder& south = incoming_waves[SOUTH];
        const PlaneWaveHolder& east =  incoming_waves[EAST];
        const PlaneWaveHolder& west =  incoming_waves[WEST];

        // Convert mpole expansion into PlaneWave expansions (up/down/north/south/east/west).
        assert(lpole_ptr.get() != NULL);
        fmm::convert_and_add_accumulated_planewaves<1,0,1>(fmm_globs, level_globs_beta, up, down, north, south, east, west, *lpole_ptr);

        incoming_waves.reset(); // no need to retain this memory

        done_gathered_lpoles = true;
        check_locals_complete();

    }

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::check_locals_complete()
{
    //std::cout << super::thisIndex <<" check_locals_complete (gather) (inherit):" << done_gathered_lpoles << " " << done_inherited_lpole << std::endl;
    if (done_gathered_lpoles && done_inherited_lpole && ready_to_evaluate == false)
    {
        ready_to_evaluate = true;

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
            stashed_cb.send();
        }
        else
        {
            const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level, universe_edge_length);

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

                fmm::ylcshift<1,0,1>(child_octant_id,
                                *lpole_ptr,
                                static_cast<MultipoleHolder&>(*local_out),
                                *workspace1,
                                *workspace2,
                                level_globs_beta.local_shift_coefficients,
                                rotation_matrix);
                super::thisProxy[child_idxer].inherit_lpole_from_parent(local_out);
            }
        }
    }

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::inherit_lpole_from_parent(MultipoleMsg *lpole_msg)
{
    *lpole_ptr += *lpole_msg;
    delete lpole_msg;

    done_inherited_lpole = true;

    check_locals_complete();
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::make_plane_waves(int dummy)
{

    // if have not made interaction list yet, do so now.
    if (expected_incoming_waves == -1)
    {
        make_interaction_list();
    }
    
    // if there is nowhere to send stuff, don't bother creating plane wave expansions
    if (my_outgoing_interaction_list.size() != 0) 
    { 

        const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
        const NodeT& node = tree.get_node(super::thisIndex);
        OctreeIndexer parent_idxer = super::thisIndex.get_parent_idx();

        // Convert mpole expansion into PlaneWave expansions (up/down/north/south/east/west).
        const FMM_Globals<NTERMS>& fmm_globs = FMM_Globals_Proxy.ckLocalBranch()->get_fmm_globs();
        const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& level_globs_beta = FMM_Globals_Proxy.ckLocalBranch()->get_globs_for_beta_level(beta, level-1, universe_edge_length);

        boost::scoped_ptr< MultiHolder<6, PlaneWaveHolder> > pw_ptr(new MultiHolder<6, PlaneWaveHolder>); // cast the message as the holder type
        MultiHolder<6, PlaneWaveHolder>& pw_out = *pw_ptr;
        PlaneWaveHolder& up =    pw_out[UP];
        PlaneWaveHolder& down =  pw_out[DOWN];
        PlaneWaveHolder& north = pw_out[NORTH];
        PlaneWaveHolder& south = pw_out[SOUTH];
        PlaneWaveHolder& east =  pw_out[EAST];
        PlaneWaveHolder& west =  pw_out[WEST];

        fmm::convert_mpole_to_six_planewaves<1,0,1>(fmm_globs, level_globs_beta, *mpole_ptr, up, down, north, south, east, west);
        
        // set the pending counter and target array within the message, then send it
        // the first target location is the last one in the list, after this the 
        // receive_incoming_wave function is responsible for figuring out where to
        // pass the message onto next (it looks for entries in pw_msg->targets array
        // which are locally accessible.
        size_t interactees = my_outgoing_interaction_list.size();
        if (interactees)
        {
            
            CkVec<CkArrayIndexMax> elems;    // add array indices
            for (unsigned int ii=0; ii < interactees; ++ii)
            {
                assert(ii < (6*6*6 - 27));
                //elems.push_back(my_outgoing_interaction_list[ii]);
                process_interaction(pw_out, my_outgoing_interaction_list[ii].first, my_outgoing_interaction_list[ii].second);
            }
            //SectionProxyType sectProxy = SectionProxyType::ckNew(super::ckGetArrayID(), elems.getVec(), elems.size());

            //CkMulticastMgr *mCastGrp = MulticastMgrProxy.ckLocalBranch();
            //sectProxy.ckSectionDelegate(mCastGrp);  // initialize section proxy
            //sectProxy.receive_incoming_wave(pw_msg);
        }
        
        made_plane_waves = true;
        if (passed_upwards && made_plane_waves) { mpole_ptr.reset(); }
        
    }
    
    // sync point
    CkCallback cb(this, IdxType::request_collected_waves());
    super::contribute(cb);
    

    return;

}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::pass_multipole_upwards()
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
    const DblMatrix3D& rotation_matrix =  (id > 4) ? fmm_globs.rdmsq3 : fmm_globs.rdsq3;

    fmm::ympshift<1,0,1>(id,
                        *mpole_ptr,
                        parent_mexp,
                        *workspace1,
                        *workspace2,
                        mp_translation_coeffs_beta,
                        rotation_matrix);

    // send up
    //std::cout << node.get_idx().get_parent_idx() << " Should get:" << *parent_mexp << std::endl;
    super::thisProxy[node.get_idx().get_parent_idx()].receive_multipole_contribution_from_child(msg);

    passed_upwards = true;
    if (passed_upwards && made_plane_waves) { mpole_ptr.reset(); }

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::receive_multipole_contribution_from_child(MultipoleMsg *mpole_msg)
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
        super::thisProxy[super::thisIndex].make_plane_waves(0);
        pass_multipole_upwards();

/*        // whilst that is going on, we could do something useful such as convert mpoles to plane waves
        CkEntryOptions* opts = new CkEntryOptions;
        opts->setPriority(5);
        super::thisProxy[super::thisIndex].make_plane_waves(0, opts);*/
    }

    return;
}

} // end namespace

#define CK_TEMPLATES_ONLY
#include "vanilla_fmm_worker.def.h"
#undef CK_TEMPLATES_ONLY


#endif //__VANILLA_FMM_WORKER_H__
