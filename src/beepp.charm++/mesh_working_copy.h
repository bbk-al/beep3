#ifndef __MESH_WORKING_COPY_H__
#define __MESH_WORKING_COPY_H__

#include "prerequisites.h"
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "vanilla_fmm.decl.h"
#include "vanilla_fmm.h"

class MeshWorkingCopy : public CBase_MeshWorkingCopy 
{

private:

    int mesh_library_idx;
    unsigned int numbering_offset;
    const Mesh& get_ref_mesh() const;
    
    Vector xyz;
    Quaternion rot;
    double Dprotein;
    unsigned int total_qual_pts;
    unsigned int qual_pts_per_triangle;
    unsigned int quad_pts_per_triangle;
    
    CProxy_VanillaFMM<18,18,300> VanillaFMMProxy;
    CkCallback stashed_cb;

public:

    /// Constructors ///
    MeshWorkingCopy();
    MeshWorkingCopy(CkMigrateMessage *msg);

    /// entry methods ///
    void init(unsigned int, const Quaternion, const Vector, unsigned int, double, double, unsigned int, unsigned int); // sets the index of the template mesh in the mesh library for this working copy
    void calculate_rhs(CProxy_VanillaFMM<18,18,300>, CkCallback);
    void process_returned_eval_pts(EvalMessage *msg);
                       
    void calculate_energy(double kappa, double Dsolvent);
    void write_output(std::string output_filename);

    void pup(PUP::er &p) {
    	CBase_MeshWorkingCopy::pup(p);
		p | mesh_library_idx;
		p | numbering_offset;
		p | xyz;
		p | rot;
        p | Dprotein;
        p | total_qual_pts;
        p | qual_pts_per_triangle;
        p | quad_pts_per_triangle;
        p | stashed_cb;
	}
};

#endif // __MESH_WORKING_COPY_H__ 
