/*      Author: david fallaize    Created on: 21 Jul 2010 */

/*! \file mesh.h
 * \brief This module declares the reference mesh and list classes.
 *
 * The Mesh class describes the model molecule surface in terms of
 * node patches, triangles and vertices, and also the partial charges
 * it contains.  This is used as a reference for the construction of
 * MeshInstance objects which also possess positions and orientations.
 *
 *	This module contains the following public classes:
 *	- Mesh -- the reference mesh class
 *	- ListedMesh -- a Mesh subclass with an id, stored in a MeshList
 *	- MeshList -- a list of ListedMesh, forming a reference Mesh library
 */

#ifndef MESH_H_
#define MESH_H_

#include <boost/shared_array.hpp>

#include <unordered_map>
#include "../common/math_vector.h"
#include "bem_kernels.h"
#include "node_patch.h"
#include "vertex.h"
#include "triangle.h"
#include "edge.h"
#include "../common/charge.h"
#include "spharm.h"
#include <exception>
#include <boost/shared_ptr.hpp>
#include "bezier.h"
#include "png1.h"
#include <fstream>
#include "mesh_tarball.h"

#ifdef __CHARMC__
#include <pup_stl.h>
#endif

namespace fs = boost::filesystem;	// Easier to swap to std::filesystem in '17

class Mesh {

public:

	//! Default constructor
	Mesh() :
		centre(Vector(0,0,0)),
		total_planar_area(0),
		total_bezier_area(0),
		net_charge(0),
		radius(0),
		done_energy_precalcs(false),
		quad_points_per_triangle(0),
		qual_points_per_triangle(0)
	{}

	//! Constructor from mesh-tar-zip file
	Mesh(const std::string& fname, bool force_planar=false) :
		done_energy_precalcs(false)
	{
		init(fname, force_planar);
	}

	//! Constructor from separate mesh and xyzq files
	Mesh(const std::string& mesh_filename,
	     const std::string& xyzq_filename,
	     bool force_planar=false);

	//! Copy constructor
	Mesh(const Mesh& other);

#ifdef PRETRIPTR
	//! Destructor
	virtual ~Mesh();
#endif // PRETRIPTR

#ifdef __CHARMC__
    virtual void pup(PUP::er &p) {

        p | node_patches;
        p | charges;
        p | centre;
        p | radius;

    }
#endif

    // exception that might be thrown if mesh cannot be init'd
    class MeshError : public std::exception {
    public:
        MeshError() : std::exception() {}
    };

	// get_ methods
    const std::vector<BasicNodePatch>& get_node_patches() const {
		return node_patches;
	}
    std::vector<BasicNodePatch>& get_node_patches() { return node_patches; }

#ifdef PREHYDROPHOBIC
    const std::vector<Charge>& get_charges() const { return charges; }
    std::vector<Charge>& get_charges() { return charges; }
#else
    const std::vector<Charge>& get_charges(bool all = false) const {
		return all ? allCharges : charges;
	}
    std::vector<Charge>& get_charges(bool all = false) {
		return all ? allCharges : charges;
	}
#endif // PREHYDROPHOBIC

    const std::vector<Triangle>& get_triangles() const { return triangles; }
    std::vector<Triangle>& get_triangles() { return triangles; }

#ifdef PRETRIPTR
    const std::vector<BasicTriangle*>& get_triangle_ptrs() const {
#else
    const std::vector<std::shared_ptr<Triangle>>& get_triangle_ptrs() const {
#endif // PRETRIPTR
		return triangle_ptrs;
	}
#ifdef PRETRIPTR
    std::vector<BasicTriangle*>& get_triangle_ptrs() { return triangle_ptrs; }
#else
    std::vector<std::shared_ptr<Triangle>>& get_triangle_ptrs() {
		return triangle_ptrs;
	}
#endif // PRETRIPTR

    const std::vector<Vertex>& get_vertices() const { return vertices; }
    std::vector<Vertex>& get_vertices() { return vertices; }

    inline double get_radius() const { return radius; }
    inline const Vector& get_centre() const { return centre; }
    unsigned int get_quad_points_per_triangle() const {
		return quad_points_per_triangle;
	}
    unsigned int get_qual_points_per_triangle() const {
		return qual_points_per_triangle;
	}
	
	// get_ index methods
    inline const BasicNodePatch& get_node_patch(unsigned int index) const;
#ifdef PREHYDROPHOBIC
    inline const Charge& get_charge(unsigned int index) const;
#else
    inline const Charge& get_charge(unsigned int index, bool all = false) const;
	inline unsigned int get_npcount(unsigned int ch_idx) const;
#endif // PREHYDROPHOBIC
    inline const Triangle& get_triangle(unsigned int index) const;
    inline const Vertex& get_vertex(unsigned int index) const;
#ifdef PRETRIPTR
    inline const BasicTriangle* get_triangle_ptr(unsigned int index) const;
#else
    inline const Triangle& get_triangle_ptr(unsigned int index) const;
#endif // PRETRIPTR

	// get counts
    unsigned int get_num_vertices() const { return vertices.size(); }
    unsigned int get_num_triangles() const { return triangles.size(); }
    unsigned int get_num_node_patches() const { return node_patches.size(); }
#ifdef PREHYDROPHOBIC
    unsigned int get_num_charges() const { return charges.size(); }
#else
    unsigned int get_num_charges(bool all = false) const {
		return all ? allCharges.size() : charges.size();
	}
#endif // PREHYDROPHOBIC
    unsigned int len() const { return get_num_node_patches(); }

	// get - other
    std::vector<Charge> get_ecm_charges(const std::string& ecm_filename,
									    const Quaternion&, const Vector&) const;

	// set_ methods
    inline void set_centre(const Vector& new_centre) {
        centre = new_centre;
    }
    inline void set_quad_points_per_triangle(unsigned int quad_points);
    inline void set_qual_points_per_triangle(unsigned int qual_points);

    // this method returns the centre and edge length of a cube such that all
    // vertices defined in the Mesh lie within the cube.
    // (note that there is no rotational optimization involved here
    // -- so this is probably not the smallest volume cube in which
    // the mesh can possibly fit.  This is meant to be used only to figure out
    // the scaling required to fit this Mesh into a unit Octree cube.)
    inline void get_bounding_cube(Vector &ccentre, double &edge_length) const;
    inline void get_bounding_cube_limits(Vector &max, Vector& min) const;

	// add a vertex to the mesh
    void add_vertex(const Vector& v, const Vector& vn) {
        vertices.push_back(Vertex(v,vn));
    }

	// add a triangle to the mesh
    inline void define_triangle
		(const unsigned int v1, const unsigned int v2, const unsigned int v3);

    void calculate_vertex_normals();

    std::string kinemage_node_patches() const;
    std::string kinemage_fh_vals
		(double fscale, double hscale, int num_colours) const;

    inline void reorder_vertex_triangles();

    // create Bezier curved patches for each triangle element
    // These will be used to create quasi-curved node patches
    // CurvedTriangleType should be Triangle (for planar); HybridBezierPatch
    // or PNG1_Triangle for curved
    template<typename CurvedTriangleType>
    inline void create_bezier_triangles();

    // Create BasicNodePatch objects from the triangulated surface mesh
    void create_node_patches();

	// This method sets the centre of the mesh by finding centroid of charges.
	// So obviously you should have defined the charges in the mesh first...!
	// (Better: use the init_centre() function to get the centre of the mesh
	// from an external file -- which lets you use e.g. hydropro to get a better
	// estimate of where the real centre of the molecule is, from a diffusey
	// point of view).
    Vector calculate_charges_centroid();

    Vector calculate_patches_centroid();

    double calculate_energy(
		double kappa,
        double Dprotein,
        double Dsolvent,
        const double fvals[],
        const double hvals[]) const;
    
	// Calculate electrostatic energy only for mesh in isolation
    double calculate_energy(
		double kappa,
		double Dprotein,
		double Dsolvent) const;

    KahanVector calculate_qe_force(
		double Dprotein,
        double Dsolvent,
        const double fvals[],
        const double hvals[]) const;

    void calculate_surface_integral_forces(
		double kappa,
        double Dprotein,
        double Dsolvent,
        const double fvals[],
        const double hvals[],
        KahanVector& MST_external,
        KahanVector& MST_internal,
        KahanVector& dbf,
        KahanVector& ionic) const;

    void calculate_forces(
		double kappa,
		double Dprotein,
		double Dsolvent,
		Vector& qE,
		Vector& MST_ext,
		Vector& MST_int,
		Vector& dbf,
		Vector& ionic) const;

    inline void set_dielectric_ratio(double Dsolvent, double Dprotein);

    double calculate_volume() const;
//TODO NB this should go back to being private, but add in reset_ if needed
    void init_energy_precalcs();

private:

    // init a Mesh object from tarball filename
    void init(const std::string&, bool force_planar=false);

    void init_mesh(const fs::path& mesh_filename);
    void init_charges(const fs::path& xyzq_filename);
    void init_fh_vals(const fs::path& fh_filename);
    void init_centre(const fs::path& centre_filename);

    void read_energy_precalcs(const fs::path& energies_filename);
    void renormalise_energy_precalcs();

    double calculate_radius();

    Vector centre;

    std::vector<Triangle> triangles;
#ifdef PRETRIPTR
    std::vector<BasicTriangle*> triangle_ptrs;
#else
    std::vector<std::shared_ptr<Triangle>> triangle_ptrs;
#endif // PRETRIPTR
    std::vector<Vertex> vertices;
    std::vector<BasicNodePatch> node_patches;
    std::vector<Charge> charges;
#ifndef PREHYDROPHOBIC
	// This does duplicate charges, but charges is needed for backwards 
	// compatibility.  Could have just kept neutrals here, but two lists
	// then have to be checked whenever hydrophobicity is being considered
	// Perhaps ideally would have neutrals then charges in one vector,
	// but how to do that without disturbing charges?
	// Using shared_ptr<Charge> messes up get_charges() and external users.
	// An unfortunate example of bad design:  List classes should be used.
	// On the other hand, most of the access needed is to an iterator, and it
	// should all be read-only...so is there another way around this?
	std::vector<Charge> allCharges;  // includes 0's
	// Map of allCharges - needs hash and equality functions...
	struct hashFunc{
		size_t operator()(const Charge& c) const{
			size_t s = 0;
			// Not expecting two charges in same place, so location is enough:
			for (int i = 0; i < 3; i++) boost::hash_combine(s, c(i));
			return s;
		}
	};
	struct equalsFunc{
		bool operator()(const Charge& lhs, const Charge& rhs) const{
			return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
		}
	};
	std::unordered_map<Charge, unsigned int, hashFunc, equalsFunc>
		allChargesMap;
	// Need to maintain number of node patches per charge for LJ calculation
	std::vector<unsigned int> nppc;
#endif  // PREHYDROPHOBIC

    double total_planar_area;
    double total_bezier_area;
    double net_charge;
    double radius;
    bool done_energy_precalcs;
    
    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;
};

// ListedMesh
#ifdef DELETED
typedef std::vector< boost::shared_ptr<Mesh> > MeshList;
#else
class ListedMesh : public Mesh
{
public:
	ListedMesh(unsigned int id, const std::string& filename, bool force_planar)
	: Mesh(filename, force_planar), lib_id(id) {};

	unsigned int get_id() const { return lib_id; };
private:
    const unsigned int lib_id;
};
	
// MeshList
class MeshList : public std::vector< boost::shared_ptr<ListedMesh> >
{
public:
	MeshList() {}

	boost::shared_ptr<ListedMesh>
		add(const std::string& filename, bool force_planar);

	void reset_energy_precalcs();

private:

};
#endif

// inlined Mesh methods
inline const BasicNodePatch& Mesh::get_node_patch(unsigned int index) const {
	if (index >= node_patches.size())
		throw std::out_of_range("index out of range for mesh");

	return node_patches[index];
}

inline const Triangle& Mesh::get_triangle(unsigned int index) const {
	if (index >= triangles.size())
		throw std::out_of_range("index out of range for mesh");

	return triangles[index];
}

#ifdef PREHYDROPHOBIC
inline const Charge& Mesh::get_charge(unsigned int index) const {
	if (index >= charges.size())
		throw std::out_of_range("index out of range for mesh");

	return charges[index];
}
#else
inline
const Charge& Mesh::get_charge(unsigned int index, bool all) const {
	if (index >= (all ? allCharges.size() : charges.size()))
		throw std::out_of_range("index out of range for mesh");

	return all ? allCharges[index] : charges[index];
}
inline unsigned int Mesh::get_npcount(unsigned int ch_idx) const {
	if (ch_idx > nppc.size())
		throw std::out_of_range("charge index out of range for mesh");
	return nppc[ch_idx];
}
#endif // PREHYDROPHOBIC

inline const Vertex& Mesh::get_vertex(unsigned int index) const {
	if (index >= vertices.size())
		throw std::out_of_range("index out of range for mesh");

	return vertices[index];
}

#ifdef PRETRIPTR
inline const BasicTriangle* Mesh::get_triangle_ptr(unsigned int index) const {
#else
inline const Triangle& Mesh::get_triangle_ptr(unsigned int index) const {
#endif // PRETRIPTR
	if (index >= triangle_ptrs.size())
		throw std::out_of_range("index out of range for mesh");

#ifdef PRETRIPTR
	return triangle_ptrs[index];
#else
	return *(triangle_ptrs[index]);
#endif // PRETRIPTR
}

inline void
Mesh::get_bounding_cube(Vector &ccentre, double &edge_length) const
{
	Vector max; // this will be the 'top right' corner
	Vector min; // and this will be 'bottom left'

	get_bounding_cube_limits(max, min);

	// figure out the maximum edge length in x/y/z dimension
	Vector diff = max - min;
	edge_length = diff.y > diff.x ? diff.y : diff.x;
	edge_length = diff.z > edge_length ? diff.z : edge_length;
	//edge_length = 1.0;

	// centre of the cube is the mid point of the two extremities
	ccentre = (max + min) / 2.0;
}

inline void Mesh::get_bounding_cube_limits(Vector &max, Vector& min) const
{
	// loop over node patch points in the mesh
	for (std::vector<BasicNodePatch>::const_iterator
			it=node_patches.cbegin(), end=node_patches.cend();
		 it != end; ++it)
	{
		const Vector& v = *it;
		if (it == node_patches.begin()) {
			max = v;
			min = v;
		}
		else {
			max.x = v.x > max.x ? v.x : max.x;
			max.y = v.y > max.y ? v.y : max.y;
			max.z = v.z > max.z ? v.z : max.z;

			min.x = v.x < min.x ? v.x : min.x;
			min.y = v.y < min.y ? v.y : min.y;
			min.z = v.z < min.z ? v.z : min.z;
		}
	}
	
	// loop over charges in the mesh
	// TODO PREHYDROPHOBIC - ok with mesh check above, but should be allCharges?
	for (std::vector<Charge>::const_iterator
			it=charges.cbegin(), end=charges.cend();
		 it != end; ++it)
	{
		const Vector& v = *it;
		if (it == charges.begin() && node_patches.size() == 0) {
			max = v;
			min = v;
		}
		else {
			max.x = v.x > max.x ? v.x : max.x;
			max.y = v.y > max.y ? v.y : max.y;
			max.z = v.z > max.z ? v.z : max.z;

			min.x = v.x < min.x ? v.x : min.x;
			min.y = v.y < min.y ? v.y : min.y;
			min.z = v.z < min.z ? v.z : min.z;
		}
	}
}

inline void Mesh::define_triangle(
	const unsigned int v1,
	const unsigned int v2,
	const unsigned int v3)
{
	unsigned int t_idx = triangles.size();
	triangles.push_back(Triangle(vertices, v1,v2,v3));

	vertices[v1].add_triangle(t_idx);
	vertices[v2].add_triangle(t_idx);
	vertices[v3].add_triangle(t_idx);

	//const Triangle& t = triangles[t_idx];
}

inline void Mesh::reorder_vertex_triangles() {
	// reorder the vertices for coherence
	for (std::vector<Vertex>::iterator it=vertices.begin(), end=vertices.end();
		 it != end; ++it)
	{
		it->reorder(triangles);
	}
}

template<typename CurvedTriangleType>
inline void Mesh::create_bezier_triangles()
{
	triangle_ptrs.reserve(triangles.size());
	for (std::vector<Triangle>::const_iterator tri_it=triangles.begin(), tri_end=triangles.end(); tri_it != tri_end; ++tri_it)
	{
#ifdef PRETRIPTR
		triangle_ptrs.push_back(new CurvedTriangleType(*tri_it));
#else
		triangle_ptrs.push_back(std::make_shared<CurvedTriangleType>(*tri_it));
#endif // PRETRIPTR
	}

}

inline void Mesh::set_dielectric_ratio(double Dsolvent, double Dprotein) {
	double epsilon_ratio = Dsolvent/Dprotein;
	for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
	{
		it->set_dielectric_ratio(epsilon_ratio);
	}
}

inline void Mesh::set_quad_points_per_triangle(unsigned int quad_points) { 
	quad_points_per_triangle = quad_points; 
	for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
	{
		it->set_quad_points_per_triangle(quad_points);
	}
}

inline void Mesh::set_qual_points_per_triangle(unsigned int qual_points) { 
	qual_points_per_triangle = qual_points; 
	for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
	{
		it->set_qual_points_per_triangle(qual_points);
	}
}

#endif /* MESH_H_ */
