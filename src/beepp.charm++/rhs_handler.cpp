#include "prerequisites.h"

#include <boost/scoped_array.hpp>

#include "opencl_nodegroup.decl.h"
#include "rhs_handler.decl.h"
#include "mesh_working_copy.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "main.decl.h"

#include "rhs_handler.h"
#include "opencl_nodegroup.h"
#include "mesh_working_copy.h"

extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
extern /* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;

RHS_Handler::RHS_Handler(unsigned int max_items_per_node, Vector universe_centre, double universe_edge_length, unsigned int num_charges) : pending(0)
{
    // note that this is the parallel_fmm_octree holding Charge objects NOT NodePatch objects -- there are two different
    // parallel fmm octrees knocking about in this program, one for the RHS vector calculation, and one for the BEM/FMM.
    VanillaFMMProxy = CProxy_VanillaFMM<18,18,300>::ckNew(max_items_per_node, universe_centre, universe_edge_length, 0, CkCallback(CkCallback::ignore));
    total_num_charges = num_charges;
    return;
}

void RHS_Handler::add_charges(std::vector<Charge> charges)
{
    VanillaFMMProxy.insert(charges);
}

void RHS_Handler::get_rhs(double Dsolv, unsigned int num_patches, CkCallback cb)
{
    start_clock = myclock();
    CkPrintf("Calculating RHS...\n");

    VanillaFMMProxy.callback_when_filled(total_num_charges, CkCallbackResumeThread());
    
    // stash the callback
    total_num_patches = num_patches;
    stashed_cb = cb;
    Dsolvent = Dsolv;
    results.reset(new double[total_num_patches*2]);
    memset(results.get(),0,sizeof(double)*total_num_patches*2);

    // Now can solve the FMM...
    CkPrintf("Solving...\n");
    VanillaFMMProxy.solve(0, CkCallbackResumeThread());
    
    // Now trigger the individual meshes to evaluate their NodePatch locations for RHS vector
    CkPrintf("Evaluating...\n");
    MeshWorkingCopyProxy.calculate_rhs(VanillaFMMProxy, CkCallback(CkIndex_RHS_Handler::process_rhs_results_from_working_mesh(NULL), thisProxy));

    // pending is the number of patches we need rhs values for
    // when the counter hits zero, we have fully built the rhs vector
    // and can return it in the stashed callback
    pending = total_num_patches;
}

void RHS_Handler::process_rhs_results_from_working_mesh(EvalMessage *msg)
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
#include "vanilla_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#include "rhs_handler.def.h"
