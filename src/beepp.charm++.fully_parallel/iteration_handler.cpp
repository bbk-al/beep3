#include "prerequisites.h"

#include "opencl_nodegroup.decl.h"
#include "fmm_globals_nodegroup.decl.h"
#include "rhs_handler.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "iteration_handler.decl.h"
#include "main.decl.h"

#include "rhs_handler.h"
#include "parallel_fmm_octree.h"
#include "explicit_worker.decl.h"
#include "fh_values_nodegroup.h"
#include "fmm_worker.decl.h"
#include "opencl_nodegroup.h"
#include "iteration_handler.h"

#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>

using beepp::CProxy_FMM_Globals_NodeGroup;
using beepp::CProxy_FMMWorker;

extern /* readonly */ CProxy_Main MainProxy;
extern /* readonly */ beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
extern /* readonly */ beepp::CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
extern /* readonly */ CProxy_FMMWorker FMMWorkerProxy;
extern /* readonly */ CProxy_IterationHandler IterationHandlerProxy;
extern /* readonly */ CProxy_ExplicitWorker ExplicitWorkerProxy;
extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

IterationHandler::IterationHandler(double _universe_edge_length, double _beta) :
        universe_edge_length(_universe_edge_length), beta(_beta), beta0(1e-10), pending(0), num_workers(0)
{
    if (beta < beta0) { beta = beta0; }
    return;
}

void IterationHandler::do_bem_fmm_iteration(CkCallback cb, FH_Values fh_vals)
{
    start_time = myclock();
    pending = 0;
    stashed_callback = cb;

    // figure out how many node patches there are
    //unsigned int total_num_patches = static_cast<unsigned int>(fh_vals.size() / 2);

    //std::cout << "Setting FH values" << std::endl;
    FH_Values_NodeGroupProxy.set(fh_vals, CkCallback(CkIndex_IterationHandler::phase_two(), thisProxy));

}

void IterationHandler::phase_two()
{

    //std::cout << "Done setting FH values" << std::endl;

    //std::cout << "Solving FMM" << std::endl;
    FMMWorkerProxy.solve(beta, beta0, CkSelfCallback(CkIndex_IterationHandler::done_evals()));

    //std::cout << "Evaluating explicit neighbours" << std::endl;
    ExplicitWorkerProxy.evaluate(beta);

}

void IterationHandler::done_evals()
{
    ++pending;
    assert(pending <= num_workers);
    //std::cout << "pending: " << pending << " / " << num_workers << std::endl;
    if (pending == num_workers)
    {
        //std::cout << "Done evaluations!" << std::endl;
        OpenCL_NodeGroupProxy.collate_bem_results(CkCallback(CkIndex_IterationHandler::reduce_fh_results(), thisProxy));
    }

    return;
}

void IterationHandler::reduce_fh_results() {

    // reduce the results vector
    //std::cout << "Reducing solution across all compute nodes" << std::endl;
    FH_Values_NodeGroupProxy.reduce(CkCallback(CkIndex_IterationHandler::done_reduction(NULL), thisProxy));

}

void IterationHandler::set_num_workers(unsigned int num_workers_in)
{
    num_workers = num_workers_in;

}

void IterationHandler::done_reduction(CkReductionMsg* msg)
{
    //std::cout << "Done reduction" << std::endl;

    double* results = reinterpret_cast<double*>(msg->getData());
    const unsigned int total_num_patches = FH_Values_NodeGroupProxy.ckLocalBranch()->get_num_patches();
    assert(msg->getSize() / sizeof(double) == total_num_patches*2);

    // create the return message, and copy data into it
    RHS_Message* results_msg = new(2*total_num_patches, 0) RHS_Message;
    results_msg->length = 2*total_num_patches;
    memcpy(results_msg->data, results, sizeof(double)*2*total_num_patches);

    // send the results via the callback (Will unsuspend the calling thread)
    stashed_callback.send(results_msg);
    
    // delete the reduction message
    delete msg;
    
    std::ostringstream buf;
    buf << "Completed BEM/FMM evaluation in " << (myclock() - start_time) / 1000. << " ms\n";
    CkPrintf(buf.str().c_str());

    return;

}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
IterationHandler::IterationHandler(CkMigrateMessage *msg) {
    std::cerr << "AARGH! IterationHandler is migrating!" << std::endl;
}

#define CK_TEMPLATES_ONLY
#include "fmm_worker.def.h"
//#include "fmm_globals_nodegroup.def.h"
#undef CK_TEMPLATES_ONLY

#include "iteration_handler.def.h"
