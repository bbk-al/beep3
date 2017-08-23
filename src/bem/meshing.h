/*      Author: Adam Light        Created on: 22 Jul 2017 */

/*! \file meshing.h
 * \brief This module declares meshing class.
 *
 * The Meshing class is a kludge to enable both Mesh and MeshInstance
 * objects to use the same methods for examining and manipulating meshes.
 * The correct design would very likely have MeshInstance subclass Mesh
 * and BasicNodePatch subclass Vertex, with all this functionality in Mesh.
 * But that has far-reaching consequences, and hence this clunky approach
 * via template specialisations and lambdas.
 *
 * The methods enable:
 *	-	the construction of Mesh objects from other instances of either Mesh or
 *		MeshInstance;
 *	-	testing whether a point is internal to a Mesh or MeshInstance; and
 *	-	finding the nearest point on a Mesh or MeshInstance to a given point
 *		in a given direction.
 *
 *	Terminology:
 *	Mesh construction is from a primary mesh and an additional mesh, with
 *	the constructed mesh being the secondary mesh.
 *	Implemented: primary=SES, additional=extended SAS, secondary=HE mesh.
 *	Not implemented: primary=SES of mesh A, additional=HE of mesh B,
 *	secondary=recovered volume mesh (to calculate HE properly);
 *	primary=SES of mesh A, additional=SES of mesh B, secondary=low dielectric
 *	mesh (to reduce dielectric between nearby protein meshes).
 *
 *	This module contains the following public classes:
 *	- Meshing -- the meshing class
 */

#ifndef MESHING_H_
#define MESHING_H_

#include <unordered_map>
#include <unordered_set>
#include "../common/math_vector.h"
//#include "bem_kernels.h"
#include "node_patch.h"
#include "vertex.h"
#include "triangle.h"
#include "edge.h"
#include <exception>
#include <fstream>
#include <boost/functional/hash.hpp>

// Meshing about with types
namespace meshing {
	using Uint = unsigned int;

	template <typename T>
	using List = std::vector<T>;	// Convenience and in case vector --> ???

	template<Uint N>
	using Multiplet = std::array<Uint,N>;

	template <Uint N>
	using IdxListN = List<Multiplet<N>>;

	template <typename T> struct hashFn{
		size_t operator()(T v) const{ return boost::hash_value(v); }
	};
	template <typename T> struct eqFn {
		bool operator()(T lhs, T rhs) const{ return (lhs == rhs); }
	};
	template<Uint N, typename T = Uint>
	using IdxMapN = std::unordered_map<
				Multiplet<N>, T, hashFn<Multiplet<N>>, eqFn<Multiplet<N>> >;

	using IdxSet = std::unordered_set<Uint>;
	using IdxOrdSet = std::set<Uint>;

	template <typename T = Uint>
	using IdxMap = std::unordered_map<Uint, T>;
	using VecList = std::vector<Vector>;

	// Lazy structures - always use value initialisation, or explicitly set.
	// Description of a boundary node - midpoint of inner and outer vertices.
	struct BoundaryNode {
		Uint idx;		// Node vertex index - set in convert_boundary()
		Uint inner;		// Interior mesh vertex index
		Uint outer;		// Exterior mesh vertex index
		bool cavity;	// If true, inner of next in list defines extra triangle
		Vertex vert;	// Location and normal of this node (midpoint)
	};
	// Description of a boundary
	struct Boundary {
		List<BoundaryNode> blist;
		IdxOrdSet region;	// contained nodes by index
		Vector centroid;	// of boundary only
	};

	// Linkage from additional to primary mesh boundaries
	using APlink = std::unordered_map<Uint,Uint>;

	// An unordered set is used to detect repeat vertex hits in pt_triangle
	// Need two unary function objects:
	struct hashFunc{
		size_t operator()(const Vector& v) const{
			size_t s = 0;
			for (int i = 0; i < 3; i++) boost::hash_combine(s, v(i));
			return s;
		}
	};
	struct equalsFunc{
		bool operator()(const Vector& lhs, const Vector& rhs) const{
			return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
		}
	};
	using PtSet = std::unordered_set<Vector, hashFunc, equalsFunc>;
};

using namespace meshing;

// Class M is either Mesh or MeshInstance
template <class M>
class Meshing {

public:
	//! No default constructor
	Meshing() = delete;

	//! Constructor from primary and additional objects
	Meshing(const M& pm, const M& am);

	//! Copy constructor - currently assumed not required

	//! Destructor - default

    // exception that might be thrown if mesh cannot be init'd
    class Error : public std::exception {
    public:
        Error() : std::exception() {}
    };

	// get_ methods
	const auto& get_vertices(void) const { return vertex_list; };
	const auto& get_edges(void) const { return edge_list; };
	const auto& get_faces(void) const { return face_list; };

	//! Test if a node is the next one in the boundary
	bool test_next_node(
		const M& m,
		BoundaryNode& candidate,	//!\param in-out potential next node
		bool interior,				//!\param F=get exterior boundary
		Boundary& boundary			//!\param in-out current boundary
	);
	//! Find boundaries in this mesh within a set of node indices
	void find_boundaries(
		const M& m,
		const IdxSet& nodeSet,		//!\param Node indices set
		bool interior,				//!\param F=get exterior boundary
		List<Boundary>& boundaries	//!\param Returned boundaries
	);
	//! Match boundaries by proximity ensuring 1:1 _from_ aboundaries
	//! The additional mesh's boundaries must be no longer than this one's
	const APlink& match_boundaries(
		const List<Boundary>& aboundaries,	//!\param shorter list
		const List<Boundary>& pboundaries	//!\param longer list
	);
	//! Convert vertices, edges and faces from old to new lists and indices
	void convert_region(
		const M& m,
		const IdxOrdSet& region,	//!\param in, vertex indices
		IdxMap<>& vmap				//!\param in-out, v old to new
	);
	//! Convert vertices, edges and faces from old to new lists and indices
	void convert_boundary(
		const M& m,
		List<BoundaryNode>& boundary,	//!\param in, new nodes
		IdxMap<>& vmap					//!\param in-out, v old to new
	);
	//! Convert a single face to new lists and indices
	void convert_face(
		Uint v1, Uint v2, Uint v3	// New vertex indices
	);
	//! Link two sets of boundaries using the aplink
	void link_boundaries(List<Boundary>& aboundaries,
						 List<Boundary>& pboundaries);

	// Default map of hit vertices - must be reset by caller of pt_triangle
	static PtSet pt_hitVertices;

	// booleans
	static bool pt_triangle(const Vector& pt, const Vector* tv, int zdir,
						PtSet& hitVertices);
	static bool pt_triangle(const Vector& pt, const Vector* tv,
				const Vector& z, PtSet& hitVertices, Vector *ptr = nullptr);
	// This method has several modes of operation:
	// If dir is null, the question is only if pt is internal and a slightly
	// more efficient approach to this is used.
	// Otherwise, if *nearest == pt, then *nearest is returned with the nearest
	// crossing of m in the direction *dir;  else the first crossing found
	// that is nearer to pt than *nearest ends the search.
	// Normally, a true return means pt is internal; but in this last case of
	// return on first crossing, a true return means a nearer hit occurred.
	// Technically, dir == nullptr, nearest != nullptr is possible but
	// is not supported:  if it works it just uses the z-direction.
	// hitVertices should be explicitly set in a threaded context.
    static bool pt_is_internal(
		const M& m,								//!\param in the mesh/instance
		const Vector& pt,						//!\param in the point
		PtSet& hitVertices = pt_hitVertices,	//!\param in-out track hits
		const Vector* dir = nullptr,			//!\param in direction to use
		Vector* nearest = nullptr);				//!\param out if hit, pt on mesh

private:
	//void init(void);

	// Utility to obtain midpoint of vector and normal of two vertices
	Vertex midpoint(const Vertex& v1, const Vertex& v2) const;
	// Recursively locate best matches for match_boundaries()
	using APD = std::vector<std::vector<double>>;
	void setaplink(Uint pbi, Uint di, APD& apd);

	IdxMapN<2> emap;			//! new edge vert inds to new edge
	IdxMapN<3> tmap;			//! new face vert inds to new face
	List<Vertex> vertex_list;	//! vertices
	IdxListN<2> edge_list;		//! edge vertex indices
	IdxListN<3> face_list;		//! face vertex indices
	APlink aplink;				//! mesh pairings
	// Hacks: see note at end of link_boundaries()
	IdxMap<Multiplet<2>> efmap;	//! edge index to two face indices or max
	IdxMap<List<Uint>> vemap;	//! from vertex index to edge indices

	const M* primary;
	const M* additional;

	// vertex has to be by value:  there is no vertices member in MeshInstance
	static std::function<const Vertex(const M&, Uint)> vertex;
	static std::function<const Triangle&(const M&, Uint)> triangle;
	static std::function<const List<Uint>&(const M&, Uint)>
		vertex_triangles;
	static std::function<Uint(const M&, Multiplet<2>)> ecount;
	static std::function<Uint(const M&, Multiplet<3>)> fcount;
	static std::function<const std::vector<Triangle>&(const M&)> triangles;
	static std::function<const BasicNodePatch&(const M&, Uint idx)> patch;
	static std::function<const Vector&(const M&)> centre;
};

#endif /* MESHING_H_ */
