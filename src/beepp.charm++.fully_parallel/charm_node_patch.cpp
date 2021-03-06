#include "prerequisites.h"
#include "charm_node_patch.h"
#include "mesh_library.decl.h"
#include "mesh_library.h"

CharmNodePatch::CharmNodePatch() : BasicNodePatch(), lib_idx(0) {}


CharmNodePatch::CharmNodePatch(const CProxy_MeshLibrary& mesh_library_proxy,
                               const BasicNodePatch& other, 
                               unsigned int _lib_idx,
                               const Vector& _centre_of_rotation, 
                               const Quaternion& _rotation, 
                               const Vector& _xyz_offset) :
                                   BasicNodePatch(other), 
                                   mesh_lib_proxy(mesh_library_proxy),
                                   lib_idx(_lib_idx),
                                   centre_of_rotation(_centre_of_rotation),
                                   rotation(_rotation),
                                   xyz_offset(_xyz_offset) 
{
    // convert reference (local) coordinates to universe coords
    Vector::change_coordinate_frame(centre_of_rotation, rotation, xyz_offset);
    BasicNodePatch::node.change_coordinate_frame(centre_of_rotation, rotation, xyz_offset);
    BasicNodePatch::centroid.change_coordinate_frame(centre_of_rotation, rotation, xyz_offset);
    BasicNodePatch::normal.apply_rotation(rotation); // normal vector just gets rotated
    BasicNodePatch::alt_normal.apply_rotation(rotation);

}

void CharmNodePatch::pup(PUP::er &p) {
    
    BasicNodePatch::pup(p);

    p | BasicNodePatch::dielectric_ratio;
    p | BasicNodePatch::f;
    p | BasicNodePatch::h;
    p | BasicNodePatch::energy_coefficient_f;
    p | BasicNodePatch::energy_coefficient_h;
    p | BasicNodePatch::force_coefficient_f;
    p | BasicNodePatch::force_coefficient_h;
    p | BasicNodePatch::gc;

    p | BasicNodePatch::node;
    p | BasicNodePatch::centroid;
    p | BasicNodePatch::normal;
    p | BasicNodePatch::alt_normal;
    p | BasicNodePatch::weighted_area;
    p | BasicNodePatch::bezier_area;
    p | BasicNodePatch::planar_area;
    p | BasicNodePatch::idx;
    p | BasicNodePatch::vertex_idx;
    
    p | mesh_lib_proxy;
    p | lib_idx;
    
    p | centre_of_rotation;
    p | rotation;
    p | xyz_offset;
    
    
    if (p.isUnpacking())
    {
        BasicNodePatch::single_qp = boost::shared_ptr<QuadList>(new QuadList);
        BasicNodePatch::single_qp->push_back(QuadPoint(*this, BasicNodePatch::normal, BasicNodePatch::weighted_area, BasicNodePatch::bezier_area));
    }
}

const Mesh& CharmNodePatch::get_ref_mesh() const
{
    const MeshLibrary *mesh_library = mesh_lib_proxy.ckLocalBranch();
    assert(mesh_library != NULL); // is nodegroup chare so should always exist
    const Mesh& ref_mesh = mesh_library->get_mesh(lib_idx);
    return ref_mesh;
}

