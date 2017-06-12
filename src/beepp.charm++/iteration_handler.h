#ifndef __ITERATION_HANDLER_H__
#define __ITERATION_HANDLER_H__

#include "prerequisites.h"
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include "../fmm/fmm_octree.h"
#include "rhs_handler.decl.h"
#include "rhs_handler.h"

class IterationHandler : public CBase_IterationHandler
{

public:

    /// Constructors ///
    IterationHandler(double, double);
    IterationHandler(CkMigrateMessage *msg);
    void do_bem_fmm_iteration(CkCallback cb, FH_Values fh_vals);
    void phase_two();
    void done_evals();
    void reduce_fh_results();
    void done_reduction(CkReductionMsg* msg);
    void set_num_workers(unsigned int num_workers_in);
    void pup(PUP::er &p) {
    	CBase_IterationHandler::pup(p);
        p | universe_edge_length;
        p | beta;
        p | beta0;
        p | epsilon;
        p | num_workers;
        p | pending;
        p | stashed_callback;
    }
private:

    double universe_edge_length;
    double beta;
    double beta0;
    double epsilon;
    unsigned int num_workers;
    unsigned int pending;
    CkCallback stashed_callback;
    long start_time;
};

#endif // __ITERATION_HANDLER_H__
