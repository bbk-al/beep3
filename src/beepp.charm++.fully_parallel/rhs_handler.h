#ifndef __RHS_HANDLER_H__
#define __RHS_HANDLER_H__

#include "prerequisites.h"
#include <boost/scoped_array.hpp>
#include "../fmm/fmm_octree.h"

// Vanilla-FMM parallel implementation
#include "parallel_fmm_octree.decl.h"
#include "vanilla_fmm_prerequisites.h"
#include "vanilla_fmm_worker.decl.h"
#include "vanilla_fmm_worker.h"
#include "vanilla_fmm_evals.decl.h"
#include "vanilla_fmm_evals.h"

static const double very_small_number = 1e-10;

class RHS_Handler : public CBase_RHS_Handler
{

public:

    /// Constructors ///
    RHS_Handler(unsigned int max_items_per_node, Vector universe_centre, double universe_edge_length);
    RHS_Handler(CkMigrateMessage *msg);
    void add_charges(std::vector<Charge> charges);
    void get_rhs(double Dsolvent, unsigned int total_num_patches, CkCallback cb);
    void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message *msg);
    
    void pup(PUP::er &p) {
    	CBase_RHS_Handler::pup(p);
        //p | pts;
        //p | normals;
        p | quallocation_points;
        p | indexes;
    }
    
private:

//    std::vector<Vector> pts;
//    std::vector<Vector> normals;
//    std::vector<Vector> weights;
    std::vector<QuadPoint> quallocation_points;
    std::vector<unsigned int> indexes;

    // these are vanilla-fmm related parallel objects which are 
    // controlled by the RHS_Handler -- i.e. it instantiates the
    // parallel objects, puts charges into the vanilla-fmm tree
    // then triggers it to solve.  Then the potential/field
    // is evaluated at the node patch locations as required by the
    // BEM rhs vector.
    vanilla_fmm::CProxy_ParallelTree VanillaTreeProxy;
    vanilla_fmm::CProxy_FMM_Globals_NodeGroup Vanilla_FMM_Globals_Proxy;
    vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;
    vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;

    long start_clock;
    int pending;
    CkCallback stashed_cb;
    double Dsolvent;
    unsigned int total_num_patches;
    
    boost::scoped_array<double> results;
    
};

class RHS_Message : public CMessage_RHS_Message {

public:

    int length;
    double* data;

};

#endif //__RHS_HANDLER_H__
