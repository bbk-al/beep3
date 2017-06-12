/*
* node_patch.cpp
*
*  Created on: 21 Jul 2010
*      Author: david
*/
#include "../common/math_vector.h"
#include "node_patch.h"
#include "mesh.h"
#include "triangle.h"
#include "vertex.h"
#include "bem_kernels.h"
#include "mesh_instance.h"

// BasicNodePatch::BasicNodePatch(const BasicTriangle& tri) : Vector(0,0,0), dielectric_ratio(0), f(0), h(0), vertex_idx(0), idx(0)
// {
//     PointNormal pn = tri.get_point_normal(1./3., 1./3.);
//     node = pn.pt();
//     static_cast<Vector&>(*this) = node;
//     normal = tri.get_planar_normal();
// 
//     bezier_area = tri.get_area();
//     planar_area = tri.get_planar_area();
//     weighted_area = bezier_area;
// 
//     tri.get_quad_points(quad_points, BasicTriangle::gauss_legendre_pts_1, 0);
//     tri.get_quad_points(quallocation_points, BasicTriangle::gauss_legendre_pts_1, 0);
//     //tri.get_quad_points(energy_quad_points, BasicTriangle::gauss_legendre_pts_4, 0);
//     tri.get_quad_points(galerkin_points, BasicTriangle::gauss_legendre_pts_2, 1);
// 
//     // loop over quadrature points to get the weighted area * normal vector
//     alt_normal = Vector(0,0,0);
//     centroid = Vector(0,0,0);
//     for (QuadList::const_iterator it=galerkin_points.begin(), end=galerkin_points.end();
// 	    it != end;
// 	    ++it)
//     {
// 	const QuadPoint& qp = *it;
// 
// 	// check for eccentric node patches -- should've been filtered out
// 	// by pre-processing stages.
// 	assert((node - qp).length() <= 10);
// 
// 	alt_normal += qp.normal() * qp.weight();
// 	weighted_area += qp.weight();
// 	centroid += qp * qp.weight();
//     }
//     
//     centroid /= weighted_area;
// 
//     gc = 0.5;
//     
// }

BasicNodePatch::BasicNodePatch(unsigned int _vertex_idx,
                               const Mesh& mesh) : 
                               Vector(0,0,0), 
                               mesh_ptr(&mesh), 
                               dielectric_ratio(0), 
                               f(0), h(0), 
                               vertex_idx(_vertex_idx)
{
    // mutex (thread-locking for generating quad points)
    mutex_ptr = new boost::mutex;
    
    force_coefficient_f = Vector(-1,-1,-1);
    force_coefficient_h = Vector(-1,-1,-1);

    
    const Vertex& v = mesh.vertices[vertex_idx];
    node = v;
    normal = v.get_normal();

    bezier_area = 0.0;
    planar_area = 0.0;

    for (std::vector<unsigned int>::const_iterator it=v.get_triangle_indices().begin(),
        end=v.get_triangle_indices().end();
        it != end;
        ++it)
    {
        const BasicTriangle& tri = *(mesh.triangle_ptrs[*it]);
        planar_area += tri.get_planar_area() / 3.0;
        //bezier_area += tri.get_area() / 3.0;
    }

    boost::shared_ptr<QuadList> temp_points = generate_points(BasicTriangle::gauss_legendre_pts_4());
    assert(temp_points->size() > 0);

    // loop over quadrature points to get the weighted area * normal vector
    alt_normal = Vector(0,0,0);
    centroid = Vector(0,0,0);
    for (QuadList::const_iterator it=temp_points->begin(), end=temp_points->end();
        it != end;
        ++it)
    {
        const QuadPoint& qp = *it;

        alt_normal += qp.normal() * qp.weight();
        bezier_area += qp.weight();
        centroid += qp.pt() * qp.weight();
    }
    normal = alt_normal.normalised();

    centroid /= bezier_area;
    static_cast<Vector&>(*this) = centroid;

    gc = 0.5;
    //calc_solid_angle();

    single_qual_pt = boost::shared_ptr<QuadList>(new QuadList);
    single_quad_pt = boost::shared_ptr<QuadList>(new QuadList);
    try {
        single_qual_pt->push_back(QuadPoint(*this, normal, 1.0));
        single_quad_pt->push_back(QuadPoint(*this, normal, bezier_area));
    }
    catch (BadQuadPoint)
    {
        std::cerr << "Bad single QP on node patch: " << *this << " vertex normal: " << v.get_normal() << std::endl;
        throw std::exception();
    }
    
    quad_points_per_triangle = mesh.get_quad_points_per_triangle();
    qual_points_per_triangle = mesh.get_qual_points_per_triangle();
    galerkin_points_per_triangle = 4;
    //galerkin_points_per_triangle = (quad_points_per_triangle > galerkin_points_per_triangle) ? quad_points_per_triangle : galerkin_points_per_triangle;
    //galerkin_points_per_triangle = (qual_points_per_triangle > galerkin_points_per_triangle) ? qual_points_per_triangle : galerkin_points_per_triangle;
    
    return;
}

    
BasicNodePatch::BasicNodePatch(const BasicNodePatch& other) : Vector(other), mesh_ptr(other.mesh_ptr)
{
    // init mutex
    mutex_ptr = new boost::mutex;
    
    dielectric_ratio = other.dielectric_ratio;
    f = other.f;
    h = other.h;
    energy_coefficient_f = other.energy_coefficient_f;
    energy_coefficient_h = other.energy_coefficient_h;
    force_coefficient_f = other.force_coefficient_f;
    force_coefficient_h = other.force_coefficient_h;
    gc = other.gc;

    node = other.node;
    centroid = other.centroid;
    normal = other.normal;
    alt_normal = other.alt_normal;

    planar_area = other.planar_area;
    bezier_area = other.bezier_area;

    idx = other.idx;
    vertex_idx = other.vertex_idx;
    
    single_qual_pt = other.single_qual_pt;
    single_qual_pt = boost::shared_ptr<QuadList>(new QuadList);
    single_qual_pt->push_back(*(other.single_qual_pt->begin()));

    single_quad_pt = other.single_quad_pt;
    single_quad_pt = boost::shared_ptr<QuadList>(new QuadList);
    single_quad_pt->push_back(*(other.single_quad_pt->begin()));

    quad_points_per_triangle = other.quad_points_per_triangle;
    qual_points_per_triangle = other.qual_points_per_triangle;
    galerkin_points_per_triangle = other.galerkin_points_per_triangle;

}

void BasicNodePatch::calc_solid_angle()
{
    gc = 0.0;
    boost::shared_ptr<QuadList> pts = get_galerkin_points();
    for (QuadList::const_iterator it=pts->begin(), end=pts->end(); it != end; ++it)
    {
        const QuadPoint& qp = *it;
        
        // solid angle is integral of ( n.dA / r2 )
        Vector rr(Vector(qp.pt()) - node);
        double r2 = rr.length2();
        rr /= sqrt(r2);
        gc += (rr.dot(qp.normal())*qp.weight()) / r2;
  
    }
    gc *= ONE_OVER_4PI;
    gc = 0.5 + gc;
    //std::cout << gc << std::endl;
    return;
    
}

boost::shared_ptr<QuadList> BasicNodePatch::generate_points(const QuadratureRule& rule, int subdivides) const
{
    boost::shared_ptr<QuadList> quads(new QuadList);

    const Mesh& mesh = get_ref_mesh();
    const Vertex& v = mesh.vertices[vertex_idx];
    const std::vector<BasicTriangle*>& triangles = mesh.triangle_ptrs;

    for (std::vector<unsigned int>::const_iterator it=v.get_triangle_indices().begin(),
            end=v.get_triangle_indices().end();
            it != end;
            ++it)
    {
        const BasicTriangle& tri = *(triangles[*it]);
 
        try {
            tri.get_constant_basis_quad_points(v, *quads, rule, subdivides);
        }
        catch (BadQuadPoint)
        {
            std::cerr << "Ignoring bad triangle: " << tri << " (this vertex: " << v << ")" << std::endl;
        }
    }

    // for derived classes this will make sure they're in the correct
    // coordinate frame
    change_coordinate_frame(*quads);

    return quads;
}

void BasicNodePatch::get_edge_points(std::vector<PointNormal>& edge_points) const 
{
    const Mesh& mesh = get_ref_mesh();
    const Vertex& v = mesh.vertices[vertex_idx];
    const std::vector<BasicTriangle*>& triangles = mesh.triangle_ptrs;

    // create connectivity in terms of indices
    std::vector<Vector> connectivity;

    // create list of vertices in order which define the node patch connectivity
    for (std::vector<unsigned int>::const_iterator it=v.get_triangle_indices().begin(),
            end=v.get_triangle_indices().end();
            it != end;
            ++it)
    {
        const Triangle& tri = dynamic_cast<Triangle&>(*(triangles[*it]));

        Vector a;
        Vector b;
        try {
            tri.get_anti_clockwise_vertices_from(v, a, b);
        }
        catch (Vertex::BadVertex &bv) {
            std::cout << "Vertex " << v << " not a member of Triangle " << tri << ".  Puzzling ... " << std::endl;
        }
        
        Vector mid_a = (v + a) / 2.0;
        Vector centre = (v + a + b) / 3.0;
        edge_points.push_back( PointNormal(mid_a,  tri.get_planar_normal()) );
        edge_points.push_back( PointNormal(centre, tri.get_planar_normal()) );
    }

    return;
}

Vector BasicNodePatch::calculate_force(double f,
                                       double h,
                                       double epsilon_int,
                                       double epsilon_ext,
                                       double kappa) const
{
    Vector force(0,0,0);
    double dielectric_ratio = epsilon_ext / epsilon_int;
    boost::shared_ptr<QuadList> galerkin_points = get_galerkin_points();
    for (QuadList::const_iterator it=galerkin_points->begin(), end=galerkin_points->end(); it != end; ++it)
    {
        const QuadPoint& qp = *it;
        assert(false);
    }
    return force;
}

NodePatch::NodePatch(const BasicNodePatch& other, const MeshInstance& mesh_instance, unsigned int _idx) : BasicNodePatch(other), minst(mesh_instance)
{
    BasicNodePatch::idx = _idx;
    
    // get transformation from local (reference Mesh) coordinates into Universe
    // coordinate frame (i.e. the coordinate frame in whcih the MeshInstance exists)
    const Vector& centre_of_rotation_old_frame = get_ref_mesh().get_centre();
    const Quaternion& rot_from_old_to_new = minst.get_rotation();
    const Vector& centre_of_rotation_new_frame = minst.get_xyz_offset();
    
    // convert reference (local) coordinates to universe coords
    Vector::change_coordinate_frame(centre_of_rotation_old_frame, rot_from_old_to_new, centre_of_rotation_new_frame);
    BasicNodePatch::node.change_coordinate_frame(centre_of_rotation_old_frame, rot_from_old_to_new, centre_of_rotation_new_frame);
    BasicNodePatch::centroid.change_coordinate_frame(centre_of_rotation_old_frame, rot_from_old_to_new, centre_of_rotation_new_frame);
    BasicNodePatch::normal.apply_rotation(rot_from_old_to_new); // normal vector just gets rotated
    BasicNodePatch::alt_normal.apply_rotation(rot_from_old_to_new);
    change_coordinate_frame(*BasicNodePatch::single_qual_pt);
    change_coordinate_frame(*BasicNodePatch::single_quad_pt);
    
    // set the number of qual/quad points according to the mesh instance
    BasicNodePatch::quad_points_per_triangle = mesh_instance.get_quad_points_per_triangle();
    BasicNodePatch::qual_points_per_triangle = mesh_instance.get_qual_points_per_triangle();
    
    BasicNodePatch::galerkin_points_per_triangle = 4;
    //BasicNodePatch::galerkin_points_per_triangle = (BasicNodePatch::quad_points_per_triangle > BasicNodePatch::galerkin_points_per_triangle) ? BasicNodePatch::quad_points_per_triangle : BasicNodePatch::galerkin_points_per_triangle;
    //BasicNodePatch::galerkin_points_per_triangle = (BasicNodePatch::qual_points_per_triangle > BasicNodePatch::galerkin_points_per_triangle) ? BasicNodePatch::qual_points_per_triangle : BasicNodePatch::galerkin_points_per_triangle;
    
}
    
// gets the reference mesh via the underlying MeshInstance pointer
const Mesh& NodePatch::get_ref_mesh() const
{
    return minst.get_ref_mesh();
}

void NodePatch::change_coordinate_frame(QuadList& qps) const
{
    const Vector& centre_of_rotation_old_frame = minst.get_ref_mesh().get_centre();
    const Quaternion& rot_from_old_to_new = minst.get_rotation();
    const Vector& centre_of_rotation_new_frame = minst.get_xyz_offset();
    for (QuadList::iterator it=qps.begin(), end=qps.end(); it != end; ++it)
    {
        it->change_coordinate_frame(centre_of_rotation_old_frame, rot_from_old_to_new, centre_of_rotation_new_frame);
    }
}



