#ifndef __MESH_WORKING_COPY_H__
#define __MESH_WORKING_COPY_H__

#include "prerequisites.h"
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <string>

// Vanilla-FMM parallel implementation
#include "vanilla_fmm_prerequisites.h"
#include "vanilla_fmm_worker.decl.h"
#include "vanilla_fmm_worker.h"
#include "vanilla_fmm_evals.decl.h"
#include "vanilla_fmm_evals.h"

class MeshWorkingCopy : public CBase_MeshWorkingCopy 
{

private:

    int mesh_library_idx;
    unsigned int numbering_offset;
    const Mesh& get_ref_mesh() const;
    
    Vector xyz;
    Quaternion rot;
    double Dprotein;
    
    vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;
    CkCallback stashed_cb;

public:

    /// Constructors ///
    MeshWorkingCopy();
    MeshWorkingCopy(CkMigrateMessage *msg);

    /// entry methods ///
    void init(unsigned int, const Quaternion, const Vector, unsigned int, double, double); // sets the index of the template mesh in the mesh library for this working copy
    void calculate_rhs(vanilla_fmm::CProxy_ParallelTree ParallelFMMOctreeProxy, 
                       vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy,
                       CkCallback cb);
    void process_returned_eval_pts(vanilla_fmm::Eval_Message *msg);
                       
    void calculate_energy(double kappa, double Dsolvent);
    void write_output(std::string output_filename);

    void pup(PUP::er &p) {
    	CBase_MeshWorkingCopy::pup(p);
		p | mesh_library_idx;
		p | numbering_offset;
		p | xyz;
		p | rot;
        p | Dprotein;
		p | stashed_cb;
	}
};

#endif // __MESH_WORKING_COPY_H__ 
