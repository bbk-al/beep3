#include "prerequisites.h"

#include <boost/scoped_array.hpp>

#include "opencl_nodegroup.decl.h"
#include "rhs_handler.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "mesh_working_copy.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "vanilla_fmm_worker.decl.h"
#include "vanilla_fmm_evals.decl.h"
#include "fmm_worker.decl.h"
#include "main.decl.h"

#include "parallel_fmm_octree.h"
#include "rhs_handler.h"
#include "opencl_nodegroup.h"
#include "mesh_working_copy.h"


extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
extern /* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
extern /* readonly */ ComlibInstanceHandle streaming_strat;

typedef ParallelFMMOctree<Charge, vanilla_fmm::CProxy_VanillaFMMWorker> VanillaParallelTree;

RHS_Handler::RHS_Handler(unsigned int max_items_per_node, Vector universe_centre, double universe_edge_length) : pending(0)
{
    // 'Vanilla' FMM is what I'm calling the standard FMM, as opposed to the 12-fold FMM/BEM hybrid
    Vanilla_FMM_Globals_Proxy = vanilla_fmm::CProxy_FMM_Globals_NodeGroup::ckNew();
    VanillaFMMWorkerProxy = vanilla_fmm::CProxy_VanillaFMMWorker::ckNew();

    // note that this is the parallel_fmm_octree holding Charge objects NOT NodePatch objects -- there are two different
    // parallel fmm octrees knocking about in this program, one for the RHS vector calculation, and one for the BEM/FMM.
    VanillaTreeProxy = vanilla_fmm::CProxy_ParallelTree::ckNew(max_items_per_node, universe_centre, universe_edge_length,0, CkCallback(CkCallback::ignore));

    ComlibAssociateProxy(streaming_strat, VanillaTreeProxy);    
    
    return;
}

void RHS_Handler::add_charges(std::vector<Charge> charges)
{
    // the charges are put into parallel-global objects- i.e. the
    // tree structure holding the charges is replicated on all
    // nodes (this is so that the FMM workers can be distributed
    // yet lightweight as all they can access the tree structure on
    // whatever node they reside on, rather than also containing 
    // the charges themselves)
    VanillaTreeProxy.insert(charges);
}

void RHS_Handler::get_rhs(double Dsolv, unsigned int num_patches, CkCallback cb)
{
    start_clock = myclock();
	CkPrintf("Calculating RHS...\n");
    
    // stash the callback
    total_num_patches = num_patches;
    stashed_cb = cb;
    Dsolvent = Dsolv;
    results.reset(new double[total_num_patches*2]);
    memset(results.get(),0,sizeof(double)*total_num_patches*2);

    // Ensure that the parallel octree has all neighbours
    // This method is threaded, and will suspend until all neighbourhoods have done this
    VanillaTreeProxy.finalize(CkCallbackResumeThread());

    unsigned int num_workers=0;
    
    // Create worker chares
    VanillaParallelTree *parallel_octree = VanillaTreeProxy.ckLocalBranch();
    assert(parallel_octree != NULL);
    double universe_edge_length = parallel_octree->get_edge_length();
    for (unsigned short level=parallel_octree->get_top_level(); level <= parallel_octree->get_bottom_level(); ++level)
    {
        for (VanillaParallelTree::NodeList::const_iterator it=parallel_octree->get_node_list(level).begin(), end=parallel_octree->get_node_list(level).end();
                it != end;
                ++it)
        {
            const VanillaParallelTree::NodeT& node = *(it->second);
            assert(node.empty() == false);
            VanillaFMMWorkerProxy[node.get_idx()].insert(Vanilla_FMM_Globals_Proxy, VanillaTreeProxy, universe_edge_length);
            if (node.isLeaf()) {
                num_workers++;
            }
        }
    }
    VanillaFMMWorkerProxy.doneInserting();
    std::cout << "info: " << num_workers << " worker chares" << std::endl;

    // Now can solve the FMM...
    VanillaFMMWorkerProxy.solve(very_small_number, CkCallback(CkCallback::ignore));
    
    // Now trigger the individual meshes to evaluate their NodePatch locations for RHS vector
    MeshWorkingCopyProxy.calculate_rhs(VanillaTreeProxy, VanillaFMMWorkerProxy, CkCallback(CkIndex_RHS_Handler::process_rhs_results_from_working_mesh(NULL), thisProxy));

    // pending is the number of patches we need rhs values for
    // when the counter hits zero, we have fully built the rhs vector
    // and can return it in the stashed callback
    pending = total_num_patches;
}

void RHS_Handler::process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message *msg)
{
    // received rhs values from one of the working_mesh objects
    // put the values into the rhs vector
    
    size_t num_pts = msg->length;
    EvalPt* eval_pts = msg->data;
    for (size_t ii=0; ii < num_pts; ++ii)
    {
        const EvalPt& ep = eval_pts[ii];
        size_t idx = ep.get_idx();
        const Vector& hvec = ep.get_field();
        results[idx] = ep.get_potential() * ONE_OVER_4PI / Dsolvent;
        results[idx+total_num_patches] = (hvec.x + hvec.y + hvec.z) * ONE_OVER_4PI / Dsolvent;
    }
    delete msg; 

    assert(pending >= num_pts);
    pending -= num_pts;

    if (pending == 0)
    {
        // Completed RHS calculations
        // create a RHS message and send it via the callback
        int sizes[1];
        sizes[0] = 2*total_num_patches;
        RHS_Message* rhs_msg = new(sizes, 0) RHS_Message;
        rhs_msg->length = 2*total_num_patches;
        memcpy(rhs_msg->data, results.get(), sizeof(double)*rhs_msg->length);
        
        double time_taken = (myclock() - start_clock) / 1000.;
        
        std::ostringstream buf;
        buf << "Done calculating RHS: " << time_taken << " ms" << std::endl;
        CkPrintf(buf.str().c_str());

        stashed_cb.send(rhs_msg);

    }
    
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
RHS_Handler::RHS_Handler(CkMigrateMessage *msg) { assert(false); }

#define CK_TEMPLATES_ONLY
#include "parallel_fmm_octree.def.h"
#undef CK_TEMPLATES_ONLY

#include "rhs_handler.def.h"
