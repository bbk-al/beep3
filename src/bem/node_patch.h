/*      Author: david fallaize   Created on: 21 Jul 2010 */
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file node_patch.h
 * \brief This module declares the node patch classes.
 *
 * The BasicNodePatch class describes the vectors associated with a vertex:
 * node, centroid, normal, alt_normal.  It provides services for numerical
 * integration using qualocation and quadrature.
 *
 *	This module contains the following public classes:
 *	- BasicNodePatch -- the base class
 *	- NodePatch -- subclass associated with a MeshInstance
 *
 *	It also defines the typdef for PatchList as a list of BasicNodePatch.
 */

#ifndef NODE_PATCH_H_
#define NODE_PATCH_H_

#include "../common/math_vector.h"
#include "../common/octree_indexer.h"
#include "spharm.h"
#include "quad_point.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "bezier.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/thread/mutex.hpp>

// See bem_kernels:
extern double fGeometricCorrection(double geometry_correction, double epsilon);
extern double hGeometricCorrection(double geometry_correction, double epsilon);

// fwd decl.
class Mesh;
class MeshInstance;
class Vertex;
class Triangle;

//#define CENTROID_COLLOCATION

class BasicNodePatch : public Vector
{

public:
	// static methods
    inline static const QuadratureRule& get_rule(int order);

	//! Default constructor
    BasicNodePatch()
	: Vector(0,0,0),
//TODO NB This list is incomplete...other constructors?
	  mesh_ptr(NULL),
	  idx(0),
	  vertex_idx(0),
	  quad_points_per_triangle(0),
	  qual_points_per_triangle(0),
	  galerkin_points_per_triangle(0)
	{
        mutex_ptr = new boost::mutex;
    }

    BasicNodePatch(const BasicNodePatch& other);		//! Copy constructor
	BasicNodePatch& operator=(const BasicNodePatch&);	//! Copy assignment
    inline BasicNodePatch(BasicNodePatch&& other);		//! Move constructor
//	BasicNodePatch& operator=(BasicNodePatch&&);		//! Move assignment

	//! Constructor from reference mesh vertex
    BasicNodePatch(unsigned int vertex_idx, const Mesh& mesh);

	//! Destructor
    virtual ~BasicNodePatch() { delete mutex_ptr; }
    
    // get_ methods
    unsigned int get_idx() const { return idx; }
    const Vector& get_node() const { return node; }
    const Vector& get_centroid() const { return centroid; }
    const Vector& get_normal() const { return normal; }
    const Vector& get_alt_normal() const { return alt_normal; }
    double get_dielectric_ratio() const { return dielectric_ratio; }
    double get_planar_area() const { return planar_area; }
    double get_bezier_area() const { return bezier_area; }
    boost::shared_ptr<QuadList> get_single_qual_pt() const {
		return single_qual_pt;
	}
    boost::shared_ptr<QuadList> get_single_quad_pt() const {
		return single_quad_pt;
	}
    
    const QuadratureRule& get_qualocation_rule() const {
		return get_rule(qual_points_per_triangle);
	}
    const QuadratureRule& get_quadrature_rule() const {
		return get_rule(quad_points_per_triangle);
	}
    const QuadratureRule& get_galerkin_rule() const {
		return get_rule(galerkin_points_per_triangle);
	}
    
	// get - other
    boost::shared_ptr<QuadList> get_qualocation_points() const;
    boost::shared_ptr<QuadList> get_quadrature_points() const;
    boost::shared_ptr<QuadList> get_galerkin_points() const;
    boost::shared_ptr<QuadList> get_energy_quad_points() const {
		return get_galerkin_points();
	}
    double get_charge(unsigned short charge_idx) const {
        return get_fmm_charge_equiv(charge_idx);
    }
    void get_edge_points(std::vector<PointNormal>& pts) const;
    
	// set_ methods
    inline void set_idx(unsigned int new_idx) { idx = new_idx; }
    void set_dielectric_ratio(double val) { dielectric_ratio = val; }

    void set_quad_points_per_triangle(unsigned int quad_points) {
		quad_points_per_triangle = quad_points;
	}
    void set_qual_points_per_triangle(unsigned int qual_points) { 
		qual_points_per_triangle = qual_points;
	}
    void set_galerkin_points_per_triangle(unsigned int galerkin_points) {
		galerkin_points_per_triangle = galerkin_points;
	}


    inline void obtain_shared_quad_ptrs
		(std::vector<boost::shared_ptr<QuadList> >& dump_here) const;

    boost::shared_ptr<QuadList> generate_points
		(const QuadratureRule& rule, int subdivides=0) const;

	// virtual methods
    virtual const Mesh& get_ref_mesh() const  { return *mesh_ptr; }

    // pure virtual as it's very important that the derived classes handle any
    // changes in coordinate frame correctly
    virtual void change_coordinate_frame(QuadList& qps) const { return; }
    

    // Python- not so great at automagically handling C++ style inheritance :-(
    //Vector& py_vector() { return static_cast<Vector&>(*this); }

    void calc_solid_angle();
    
    Vector calculate_force(double f,
                           double h,
                           double epsilon_int,
                           double epsilon_ext,
                           double kappa) const;

	// Attributes
	//NB These *public* attributes should be moved!
    
	// Ratio external:internal or solvent:protein dielectrics...
	//NB There are set_ and get_ functions!
    double dielectric_ratio;

	// Referenced mainly in Mesh or MeshInstance
	// Initialised via BEEP::gmres (solve)
    double f;  // Potential
    double h;  // Normal derivative
	// Referenced mainly in Mesh
	// Initialised on construction of Mesh via init_energy_precalcs
    double energy_coefficient_f;
    double energy_coefficient_h;
    Vector force_coefficient_f;
    Vector force_coefficient_h;

	// Referenced mainly in BEEP and Python
    double gc;	// refers to a geometry correction based on solid angles...

protected:

	void copy(const BasicNodePatch& other);
	void move(BasicNodePatch&& other);
    inline double get_fmm_charge_equiv(unsigned short equiv_charge_idx) const;

	// Attributes accessible to subclasses
    const Mesh* mesh_ptr;
    unsigned int idx; // unique number of this node patch
    unsigned int vertex_idx; // vertex index within a Mesh object -- need to know this to build quadrature points on the fly
    Vector node;
    Vector centroid;
    Vector normal;
    Vector alt_normal;
    double bezier_area;
    double planar_area;
    
    boost::shared_ptr<QuadList> single_qual_pt;
    boost::shared_ptr<QuadList> single_quad_pt;

    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;
    unsigned int galerkin_points_per_triangle;
    
    
private:
    
    mutable boost::weak_ptr<QuadList> weak_quad_ptr;
    mutable boost::weak_ptr<QuadList> weak_qual_ptr;
    mutable boost::weak_ptr<QuadList> weak_galerkin_ptr;
    
    // permanently store qual_pts (more efficient when using
    // OpenCL to compute local integrations, as these are
    // requested repeatedly rather than just once; less efficient
    // if not using OpenCL, but not that huge an overhead).
    mutable boost::shared_ptr<QuadList> qual_pts;

    mutable boost::mutex* mutex_ptr;
    
};

//#define __LOCAL_MOVES__

// The NodePatch class describes an actual node patch located in universe
// coordinates- it is linked conceptually to a MeshInstance which holds the
// rotation/xyz_offset for the node patch.
class NodePatch : public BasicNodePatch
{
public:
    
	// No default constructor - would not make sense here
	//! Constructor from node patch and mesh instance
    NodePatch(
		const BasicNodePatch& other,
		const MeshInstance& mesh_instance,
		unsigned int _idx=0);

    NodePatch(const NodePatch&) = default;				//! Copy constructor
    NodePatch& operator=(const NodePatch&) = default;	//! Copy assignment
    NodePatch(NodePatch&&) = default;					//! Move constructor
	// No move assignment - see base class
    virtual ~NodePatch() = default;						//! Destructor

	// virtuals
    // gets the reference mesh via the underlying MeshInstance pointer
    virtual const Mesh& get_ref_mesh() const;
#ifdef __LOCAL_MOVES__
	virtual void change_coordinate_frame
		(const Vector& centre_of_rotation_old_frame,
		 const Quaternion& rot_from_old_to_new,
		 const Vector& centre_of_rotation_new_frame);
	virtual void change_ql_coordinate_frame
		(QuadList& qps,
		 const Vector& centre_of_rotation_old_frame,
		 const Quaternion& rot_from_old_to_new,
		 const Vector& centre_of_rotation_new_frame) const;
#endif //  __LOCAL_MOVES__
    virtual void change_coordinate_frame(QuadList& qps) const;
    
private:
    
    const MeshInstance& minst;
};

// PatchList
typedef std::vector< boost::shared_ptr<BasicNodePatch> > PatchList;


// inlined methods
// Move operations need to protect mutex_ptr from double free
//NB This would be a lot neater with a Pimpl approach!
//which would also make move assignment worthwhile
//possibly compensating for the extra memory management involved
//-- left for now
inline BasicNodePatch::BasicNodePatch(BasicNodePatch&& other)
: Vector(std::move(other)),
	mesh_ptr{std::move(other.mesh_ptr)},
	dielectric_ratio{std::move(other.dielectric_ratio)},
	f{std::move(other.f)},
	h{std::move(other.h)},
	energy_coefficient_f{std::move(other.energy_coefficient_f)},
	energy_coefficient_h{std::move(other.energy_coefficient_h)},
	force_coefficient_f{std::move(other.force_coefficient_f)},
	force_coefficient_h{std::move(other.force_coefficient_h)},
	gc{std::move(other.gc)},
	node{std::move(other.node)},
	centroid{std::move(other.centroid)},
	normal{std::move(other.normal)},
	alt_normal{std::move(other.alt_normal)},
	planar_area{std::move(other.planar_area)},
	bezier_area{std::move(other.bezier_area)},
	idx{std::move(other.idx)},
	vertex_idx{std::move(other.vertex_idx)},
	quad_points_per_triangle{std::move(other.quad_points_per_triangle)},
	qual_points_per_triangle{std::move(other.qual_points_per_triangle)},
	galerkin_points_per_triangle
	{std::move(other.galerkin_points_per_triangle)}
{
	move(std::move(other));
}

inline void BasicNodePatch::obtain_shared_quad_ptrs(
	std::vector<boost::shared_ptr<QuadList> >& dump_here) const
{ 
	dump_here.push_back(get_quadrature_points());
	dump_here.push_back(get_qualocation_points());
}

inline const QuadratureRule& BasicNodePatch::get_rule(int order) {
	const QuadratureRule* rule;
	switch (order) {
	case(1):
		rule = &(BasicTriangle::gauss_legendre_pts_1());
		break;
	case(4):
		rule = &(BasicTriangle::gauss_legendre_pts_4());
		break;
	case(6):
		rule = &(BasicTriangle::triangle_order_6());
		break;
	case(7):
		rule = &(BasicTriangle::gauss_legendre_pts_7());
		break;
	case(16):
		rule = &(BasicTriangle::gauss_legendre_pts_16());
		break;
	default:
		std::cerr << "Bad rule specified: " << order << std::endl;
		throw std::exception();
	}
	return *rule;
}

// Candidate for constexpr except for throwing an exception
inline double BasicNodePatch::get_fmm_charge_equiv(
	unsigned short equiv_charge_idx) const
{

	switch(equiv_charge_idx)
	{
	case 0:
		return bezier_area*f*normal.x;
		break;
	case 1:
		return bezier_area*f*normal.y;
		break;
	case 2:
		return bezier_area*f*normal.z;
		break;
	case 3:
		return bezier_area*h;
		break;
	
	// epsilon weighted values
	case 4:
		return bezier_area*f*normal.x / dielectric_ratio;
		break;
	case 5:
		return bezier_area*f*normal.y / dielectric_ratio;
		break;
	case 6:
		return bezier_area*f*normal.z / dielectric_ratio;
		break;
	case 7:
		return bezier_area*h / dielectric_ratio;
		break;

	// epsilon weighted values for kappa=0 (Gpt kernels)
	// (and non-epsilon weighted h vals)
	case 8:
		return bezier_area*f*normal.x / dielectric_ratio;
		break;
	case 9:
		return bezier_area*f*normal.y / dielectric_ratio;
		break;
	case 10:
		return bezier_area*f*normal.z / dielectric_ratio;
		break;
	case 11:
		return bezier_area*h;
		break;

		
	default:
		// should not get here
		throw std::exception();
	}

	return 0.0;
}

#endif /* NODE_PATCH_H_ */
