#include "prerequisites.h"
#include "vanilla_fmm_prerequisites.h"

#include "fh_values_nodegroup.decl.h"
#include "mesh_library.decl.h"
#include "rhs_handler.decl.h"
#include "mesh_working_copy.decl.h"
#include "main.decl.h"
#include "main.h"
#include "vanilla_fmm.decl.h"
#include "vanilla_fmm.h"

#include "mesh_library.h"
#include "fh_values_nodegroup.h"
#include "mesh_working_copy.h"

#include <fstream>

extern /* readonly */ CProxy_BEM_FMM<NT,NL,NW,CharmNodePatch> BEM_FMM_Proxy;
extern /* readonly */ CProxy_MeshLibrary MeshLibraryProxy;
extern /* readonly */ CProxy_RHS_Handler RHS_HandlerProxy;
extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;

MeshWorkingCopy::MeshWorkingCopy() : mesh_library_idx(-1), numbering_offset(0), total_qual_pts(0), qual_pts_per_triangle(-1), quad_pts_per_triangle(-1)  {}

void MeshWorkingCopy::init(unsigned int mesh_lib_idx, const Quaternion rotation, const Vector translation, unsigned int offset, double protein_dielectric, double Dsolvent, unsigned int num_quad_points, unsigned int num_qual_points)
{
    xyz = translation;
    rot = rotation;
    Dprotein = protein_dielectric;
    quad_pts_per_triangle = num_quad_points;
    qual_pts_per_triangle = num_qual_points;
    
    //CkPrintf("Init'ing working mesh (on processor %d) : mesh_instance_id=%d; mesh_lib_id=%d\n", CkMyPe(), thisIndex, mesh_lib_idx);

    //int this_working_mesh_idx = thisIndex;
    mesh_library_idx = mesh_lib_idx;
    numbering_offset = offset;

    std::vector<CharmNodePatch> patches;
    std::vector<Charge> charges;
    
    // copy list of patches from reference set in the mesh library (which is globally available)
    const Mesh& ref_mesh = get_ref_mesh();
    const std::vector<BasicNodePatch>& ref_patches = ref_mesh.get_node_patches();
    patches.reserve(ref_patches.size());
    total_qual_pts = 0;
    size_t ctr=0;
    for (std::vector<BasicNodePatch>::const_iterator it=ref_patches.begin(), end=ref_patches.end(); 
         it != end;
         ++it)
     {
         patches.push_back(CharmNodePatch(MeshLibraryProxy, 
                                          *it, 
                                          mesh_library_idx, 
                                          ref_mesh.get_centre(),
                                          rot,
                                          xyz,
                                          quad_pts_per_triangle,
                                          qual_pts_per_triangle));
        total_qual_pts += patches.rbegin()->get_qualocation_points()->size();
        ++ctr;
     }

    // copy list of patches from reference set in the mesh library (which is globally available)
    const std::vector<Charge>& ref_charges = ref_mesh.get_charges();
    charges.reserve(ref_charges.size());
    charges.insert(charges.begin(), ref_charges.begin(), ref_charges.end());

    // set the actual coordinates of the mesh -- will be a rotation and translation relative to the
    // reference mesh.
    for (std::vector<CharmChargeHolder>::iterator it=charges.begin(), end=charges.end();
        it != end;
        ++it)
    {
        it->change_coordinate_frame(ref_mesh.get_centre(), rot, xyz);
    }
    //std::cout << "Inserting " << patches.size() << " patches..." << std::endl;

    // Set the numbering of the node patches
    for (unsigned int idx=0; idx < patches.size(); ++idx)
    {
        patches[idx].set_idx(idx + numbering_offset);
        patches[idx].set_dielectric_ratio(Dsolvent / Dprotein);
    }

    BEM_FMM_Proxy.insert(patches);
    RHS_HandlerProxy.add_charges(charges);

}

const Mesh& MeshWorkingCopy::get_ref_mesh() const
{

    const MeshLibrary *mesh_library = MeshLibraryProxy.ckLocalBranch();
    assert(mesh_library != NULL); // is nodegroup chare so should always exist
    const Mesh& ref_mesh = mesh_library->get_mesh(mesh_library_idx);
    return ref_mesh;
}

void MeshWorkingCopy::calculate_rhs(CProxy_VanillaFMM<18,18,300> _vanilla_proxy,
                                    CkCallback cb)
{
    //std::cout << "MeshWorkingCopy " << thisIndex  << " calculate_rhs..." << std::endl;
    VanillaFMMProxy = _vanilla_proxy;
    
    stashed_cb = cb;
    
    // iterate over patches from reference set in the mesh library (which is globally available)
    const Mesh& ref_mesh = get_ref_mesh();
    const std::vector<BasicNodePatch>& ref_patches = ref_mesh.get_node_patches();
    
    Quaternion rot_from_local_to_uni = rot;
    const Vector& centre_local_coords = ref_mesh.get_centre();
    const Vector& centre_universe_coords = xyz;

    // loop over the incoming eval points, and form the rhs vector contributions
    
    // create a message to hold the results -- init a set of EvalPts, one per BasicNodePatch
    // which contain the index of the BasicNodePatch, this will all get sent back to the rhs handler
    // as part of the final result.
    EvalMessage* vanilla_eval_work = new(total_qual_pts,0) EvalMessage;
    vanilla_eval_work->length = total_qual_pts;
    
    size_t ep_ctr=0;
    size_t np_ctr=0;
    for (std::vector<BasicNodePatch>::const_iterator it=ref_patches.begin(), end=ref_patches.end();
        it != end;
        ++it)
    {
        // make a mutable copy
        BasicNodePatch np(*it);
        np.set_qual_points_per_triangle(qual_pts_per_triangle);
        
        boost::shared_ptr<QuadList> qps = np.get_qualocation_points();
        for (std::vector<QuadPoint>::const_iterator qit=qps->begin(), qend=qps->end();
             qit != qend; ++qit)
        {
            QuadPoint qual_pt = *qit; // make a temporary copy!
            //std::cout << qual_pt << " --> ";
            qual_pt.change_coordinate_frame(centre_local_coords, rot_from_local_to_uni, centre_universe_coords);
            //std::cout << qual_pt << std::endl;
            EvalPt* ep = new (&(vanilla_eval_work->data[ep_ctr++])) EvalPt(qual_pt.pt());
        }
    }
    
    vanilla_eval_work->cb = CkSelfCallback(CkIndex_MeshWorkingCopy::process_returned_eval_pts(NULL));
    VanillaFMMProxy.ckLocalBranch()->evaluate(vanilla_eval_work);
    
}

void MeshWorkingCopy::process_returned_eval_pts(EvalMessage *msg)
{
    // iterate over patches from reference set in the mesh library (which is globally available)
    const Mesh& ref_mesh = get_ref_mesh();
    const std::vector<BasicNodePatch>& ref_patches = ref_mesh.get_node_patches();
    Quaternion rot_from_local_to_uni = rot;
    
    // create a message to hold the results -- init a set of EvalPts, one per BasicNodePatch
    // which contain the index of the BasicNodePatch, this will all get sent back to the rhs handler
    // as part of the final result.
    EvalMessage* results = new(ref_patches.size(),0) EvalMessage;
    results->length = ref_patches.size();
    
    size_t ep_ctr=0, ii=0;
    for (std::vector<BasicNodePatch>::const_iterator it=ref_patches.begin(), end=ref_patches.end();
        it != end;
        ++it)
    {
        EvalPt* ep = new (&(results->data[ii])) EvalPt();
        ep->set_idx(ii+numbering_offset);
        
        // make a mutable copy
        BasicNodePatch np(*it);
        np.set_qual_points_per_triangle(qual_pts_per_triangle);

        boost::shared_ptr<QuadList> qps = np.get_qualocation_points();
        for (std::vector<QuadPoint>::const_iterator qit=qps->begin(), qend=qps->end();
             qit != qend; ++qit)
        {
            const QuadPoint& qual_pt = *qit;
            double wt = qual_pt.weight();
            assert(ep_ctr < msg->length);
            EvalPt& ep_qp = msg->data[ep_ctr++];
            EvalPt& output = results->data[ii];
            //std::cout << ep_qp.get_potential() << "\n";
            output.add_potential(ep_qp.get_potential() * wt);
            Vector tmp = ep_qp.get_field() * wt;
            Vector rotated_normal = qual_pt.normal();
            
            // don't forget that the reference mesh is in some different coordinate frame
            // from the actual entities floating about in the universe frame (they have some
            // relative rotation and translation stored in the rot / xyz members).  So in general
            // normal vectors have to be (purely) rotated (independent of origin) whilst physical
            // locations (i.e. xyz coordinates of quad points) have to be rotated about the centre
            // of rotation and the centre of rotation placed at the corresponding place in 
            // universe coordinate frame.
            rotated_normal.apply_rotation(rot_from_local_to_uni); 
            tmp.x *= rotated_normal.x;
            tmp.y *= rotated_normal.y;
            tmp.z *= rotated_normal.z;
            output.add_field(tmp);
        }
        ++ii;
    }
        
    delete msg;
    
    //std::cout << "MeshWorkingCopy " << thisIndex  << " recieved eval pts, returning callback..." << std::endl;
    stashed_cb.send(results);
}

void MeshWorkingCopy::calculate_energy(double kappa, double Dsolvent)
{
    const FH_Values_NodeGroup *fh_vals_grp = FH_Values_NodeGroupProxy.ckLocalBranch();
    const double *fvals = fh_vals_grp->fvals() + numbering_offset;
    const double *hvals = fh_vals_grp->hvals() + numbering_offset;

    const Mesh& ref_mesh = get_ref_mesh();
    double E = ref_mesh.calculate_energy(kappa, Dprotein, Dsolvent, fvals, hvals);

    std::ostringstream buf;
    buf << "Energy for mesh " << thisIndex << " (lib_id=" << mesh_library_idx << ") [" << numbering_offset << "/" << ref_mesh.get_node_patches().size() << "] = " << E << std::endl;
    CkPrintf(buf.str().c_str());
    
    return;
};

void MeshWorkingCopy::write_output(std::string output_filename)
{

    const FH_Values_NodeGroup *fh_vals_grp = FH_Values_NodeGroupProxy.ckLocalBranch();
    const double *fvals = fh_vals_grp->fvals() + numbering_offset;
    const double *hvals = fh_vals_grp->hvals() + numbering_offset;

    const Mesh& ref_mesh = get_ref_mesh();
    const std::vector<BasicNodePatch>& ref_patches = ref_mesh.get_node_patches();
    size_t num_patches = ref_patches.size();
        
    if (output_filename != "")
    {
        std::ofstream fh_output;
            
        std::stringstream name_buf;
        name_buf << output_filename << "." << thisIndex;
        //std::cout << "Writing to: " << name_buf.str() << std::endl;
        fh_output.open(name_buf.str().c_str(), std::ios_base::out);    
        
        for (size_t ctr=0; ctr < num_patches; ++ctr)
        {
            fh_output << fvals[ctr] << " " << hvals[ctr] << "\n";
        }
        
        fh_output.close();
    }
        
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
MeshWorkingCopy::MeshWorkingCopy(CkMigrateMessage *msg) { }

#define CK_TEMPLATES_ONLY
#include "bem_fmm.def.h"
#include "vanilla_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#include "mesh_working_copy.def.h"
