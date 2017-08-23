/*      Author: Adam Light        Created on: 22 Jul 2017 */

/*! \file meshing.cpp
 * \brief This module implements the meshing class.
 */

#include "../common/math_vector.h"
#include "meshing.h"
#include "mesh.h"
#include "mesh_instance.h"
#include "bem_kernels.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <cassert>
#include <exception>
#include <unordered_map>
#include <boost/scoped_ptr.hpp>
#include <boost/functional/hash.hpp>
#include <boost/range/combine.hpp>

using namespace meshing;

// Class attributes
// Specialisations of lambdas
// vertex -- consider a local noexcept version...
template <>
decltype(Meshing<Mesh>::vertex) Meshing<Mesh>::vertex
	= [](const Mesh& m, Uint idx) -> const Vertex
{ return m.get_vertex(idx); };

// There are no vertices in a MeshInstance... hence the return by value etc.
template<>
decltype(Meshing<MeshInstance>::vertex) Meshing<MeshInstance>::vertex
	= [](const MeshInstance& m, Uint idx) noexcept -> const Vertex
{
	BasicNodePatch& np = m.get_node_patch(idx);
	Vertex v(np.get_node(), np.get_normal());
	return v;
};

// triangle
template <>
decltype(Meshing<Mesh>::triangle) Meshing<Mesh>::triangle
	= [](const Mesh& m, Uint idx) -> const Triangle&
{ return m.get_triangle(idx); };

template <>
decltype(Meshing<MeshInstance>::triangle) Meshing<MeshInstance>::triangle
	= [](const MeshInstance& m, Uint idx) -> const Triangle&
{ return m.get_ref_mesh().get_triangle(idx); };

// triangles
template <>
decltype(Meshing<Mesh>::triangles) Meshing<Mesh>::triangles
	= [](const Mesh& m) noexcept -> const std::vector<Triangle>&
{ return m.get_triangles(); };

template <>
decltype(Meshing<MeshInstance>::triangles) Meshing<MeshInstance>::triangles
	= [](const MeshInstance& m) noexcept -> const std::vector<Triangle>&
{ return m.get_ref_mesh().get_triangles(); };

// vertex_triangles
template <>
decltype(Meshing<Mesh>::vertex_triangles) Meshing<Mesh>::vertex_triangles
	= [](const Mesh& m, Uint idx) -> const List<Uint>&
{ return m.get_vertex(idx).get_triangle_indices(); };

template <>
decltype(Meshing<MeshInstance>::vertex_triangles)
	Meshing<MeshInstance>::vertex_triangles
	= [](const MeshInstance& m, Uint idx) -> const List<Uint>&
{ return m.get_ref_mesh().get_vertex(idx).get_triangle_indices(); };

// patch
template<>
decltype(Meshing<Mesh>::patch) Meshing<Mesh>::patch
	= [](const Mesh& m, Uint idx) noexcept -> const BasicNodePatch&
{ return m.get_node_patches()[idx]; };

template<>
decltype(Meshing<MeshInstance>::patch) Meshing<MeshInstance>::patch
	= [](const MeshInstance& m, Uint idx) noexcept -> const BasicNodePatch&
{ return m.get_node_patch(idx); };

// ecount
template<>
decltype(Meshing<Mesh>::ecount) Meshing<Mesh>::ecount
	= [](const Mesh& m, Multiplet<2>ea) noexcept -> Uint
{ return m.edge_count(ea); };

template<>
decltype(Meshing<MeshInstance>::ecount) Meshing<MeshInstance>::ecount
	= [](const MeshInstance& m, Multiplet<2>ea) noexcept -> Uint
{ return m.get_ref_mesh().edge_count(ea); };

// fcount
template<>
decltype(Meshing<Mesh>::fcount) Meshing<Mesh>::fcount
	= [](const Mesh& m, Multiplet<3>fa) noexcept -> Uint
{ return m.face_count(fa); };

template<>
decltype(Meshing<MeshInstance>::fcount) Meshing<MeshInstance>::fcount
	= [](const MeshInstance& m, Multiplet<3>fa) noexcept -> Uint
{ return m.get_ref_mesh().face_count(fa); };

// centre
template<>
decltype(Meshing<Mesh>::centre) Meshing<Mesh>::centre
	= [](const Mesh& m) noexcept -> const Vector&
{ return m.get_centre(); };

template<>
decltype(Meshing<MeshInstance>::centre) Meshing<MeshInstance>::centre
	= [](const MeshInstance& m) noexcept -> const Vector&
{ return m.get_xyz_offset(); };


// Locals
static auto parr = [](std::pair<Uint, Uint>&& p) noexcept
				-> Multiplet<2>
{
	Multiplet<2> m = {std::get<0>(p), std::get<1>(p)};
	return m;
};

// Constructor
template <class M>
Meshing<M>::Meshing(const M& pm, const M& am) {
	primary = &pm;
	additional = &am;

	// reservations - can only offer a rough guide, 50% over max number of faces
	Uint fest = 3*(std::max(triangles(pm).size(), triangles(am).size()))/2;
	vertex_list.reserve((16+fest)/2);	// fewer triangles per vertex unlikely
	face_list.reserve(fest);
	edge_list.reserve(3*face_list.capacity()/2);
	// unordered maps and sets have efficient inserts, so no reservations needed
}

// Methods
template <class M>
Vertex Meshing<M>::midpoint(const Vertex& v1, const Vertex& v2) const {
	Vector v = (v1+v2)/2;
	Vector n = (v1.get_normal() + v2.get_normal())/2;
	n.normalise();
	Vertex vert(v, n);
	return vert;
}

// Supplied with only inner and outer indices for candidate
// If candidate passes tests, updates it with vertex and updates previous cavity
template <class M>
bool Meshing<M>::test_next_node(const M& m, BoundaryNode& candidate,
	bool interior, Boundary& boundary)
{
	auto& blist = boundary.blist;
	auto& centroid = boundary.centroid;
	auto& region = boundary.region;

	// Get the last node - caller guarantees there is one!
	BoundaryNode& bn = blist.back();

	// Check for neighbours in forward direction
	Multiplet<2> inside = parr(std::minmax(bn.inner, candidate.inner));
	Multiplet<2> outside = parr(std::minmax(bn.outer, candidate.outer));
	Multiplet<3> va = {0, 0, 0};
	// Same inner vertex but neighbouring outer vertices
	if (candidate.inner == bn.inner && ecount(m, outside) > 0)
		va[0] = bn.inner, va[1] = bn.outer, va[2] = candidate.outer;
	// Same outer vertex but neighbouring inner vertices
	else if (candidate.outer == bn.outer && ecount(m, inside) > 0)
		va[0] = bn.inner, va[1] = bn.outer, va[2] = candidate.inner;

	// If these form a triangle in the right sense
	if (va[0] != va[2] && fcount(m, va) > 0) {
		// Complete triangle specification
		bn.cavity = (candidate.outer == bn.outer);

std::cout << "Found a triangle! " << bn.cavity << " "
		<< va[0] << " " << va[1] << " " << va[2] << std::endl;
		// Add candidate onto boundary if it is not back at the start
		if (candidate.inner != blist.front().inner ||
			candidate.outer != blist.front().outer) {
			candidate.vert = midpoint(vertex(m, candidate.inner),
									  vertex(m, candidate.outer));
			blist.push_back(candidate);
			// Maintain centroid and region
			centroid = centroid + (candidate.vert-centroid)/blist.size();
			region.insert(interior ? candidate.inner : candidate.outer);
		}
		return true;
	}
	return false;
}

// No triangle list returned - worked out in create_triangles_from_vertices
template <class M>
void Meshing<M>::find_boundaries(const M& m, const IdxSet& nodeSet,
	bool interior, List<Boundary>& boundaries)
{
	// Create a list of (in,out) boundary indices and nodes from the nodeSet
	IdxMapN<2,Vertex> bidx;	// But vertex is not used
	for (const auto& npidx: nodeSet) {
		// Are all the neighbours of this node also in nodeSet?
		// Locate the neighbours via the triangles sharing this node
		const std::vector<unsigned int>& tri = vertex_triangles(m, npidx);
		Uint tvi[3];
		for (const auto& t: tri) {
			const Triangle& tr = triangle(m, t);
			tr.get_vertex_indices(tvi[0], tvi[1], tvi[2]);
			for (int i = 0; i < 3; i++) {
				if (tvi[i] == npidx) continue;
				if (nodeSet.count(tvi[i]) == 0) {
					Multiplet<2> k = {tvi[i], npidx};
					bidx[k] = vertex(m, npidx);	// One copy per npidx pair
				}
			}
		}
	}

	// Seek contiguous chains of nodes in the boundary set
	Multiplet<2> idx, idxo;	// [1] inner [0] outer indices, current and other
	while (bidx.size() > 0) {
		auto it = bidx.cbegin();
		idx = it->first;
		bidx.erase(it);

		// Start a new boundary with one incomplete node...
		boundaries.push_back(Boundary());
		auto& boundary = boundaries.back();
		auto& blist = boundary.blist;
		blist.push_back(BoundaryNode());
		auto& bnode = blist.back();
		bnode.inner = idx[interior];
		bnode.outer = idx[!interior];
		bnode.vert = midpoint(vertex(m, idx[0]), vertex(m, idx[1]));
		boundary.centroid = bnode.vert;
		boundary.region.insert(idx[1]);

		// Follow around the boundary
		bool found;
		do {
			found = false;
			BoundaryNode candidate;
			for (auto oit = bidx.cbegin(); oit != bidx.cend(); oit++) {
				candidate.inner = oit->first[interior];
				candidate.outer = oit->first[!interior];
				found = test_next_node(m, candidate, interior, boundary);
				if (found) {
					bidx.erase(oit);	// Invalidates oit..
					break;				// ... but only expecting one!
				}
			}
		} while (found);

		// Ensure the boundary ends meet up - will not update centroid
		if (!test_next_node(m, blist.front(), interior, boundary)) {
			// Oh dear
			std::cerr << "Error: boundary does not complete" << std::endl;
			throw std::exception();
		}
	}

	// Fill in the regions from the boundaries
	List<Uint> rlist;
	rlist.reserve(nodeSet.size());
	for (auto& boundary: boundaries) {
		auto& reg = boundary.region;
		int added = rlist.size();
		for (const auto& ridx: reg) rlist.push_back(ridx);
		unsigned int rsize;
		do {
			rsize = rlist.size();
			for (int r = added; r < rsize; r++) {
				const auto& ridx = rlist[r];
				const std::vector<unsigned int>& tri
					= vertex_triangles(m, ridx);
				Uint tvi[3];
				for (const auto& t: tri) {
					const Triangle& tr = triangle(m, t);
					tr.get_vertex_indices(tvi[0], tvi[1], tvi[2]);
					for (int i = 0; i < 3; i++) {
						if (tvi[i] == ridx) continue;
						if (reg.count(tvi[i])==0 && nodeSet.count(tvi[i])>0) {
							rlist.push_back(tvi[i]);
							reg.insert(tvi[i]);
						}
					}
				}
			}
			added = rsize;
		} while (rlist.size() > rsize);
	}
}

// Recursively locate best AP pairs in APD scores matrix (see match_boundaries)
// TODO this isn't quite right:  need to get min totals overall...
template <class M>
void Meshing<M>::setaplink(Uint pbi, Uint di, APD& apd) {
	// If this primary boundary is free, make a match
	if (aplink.count(pbi) == 0) {
		aplink[pbi] = di;
		return;
	}
	// This primary boundary has been taken, bump a boundary
	// [should really check seconds if equality, then thirds...]
	Uint odi = aplink[pbi];
	if (apd[odi][pbi] < apd[di][pbi])
		aplink[pbi] = di;
	else
		odi = di;
	apd[odi][pbi]
		= *std::min_element(apd[odi].cbegin(), apd[odi].cend())-1;
	auto it = std::max_element(apd[odi].cbegin(), apd[odi].cend());
	int npbi = it - apd[odi].cbegin();
	setaplink(npbi, odi, apd);
}

// To match boundaries a simple score is used that accounts for both the
// distance between boundaries and their overall orientation.  This is
// intended to eliminate buried regions on the primary mesh.
// However, the correct way to do this is to directly detect buried vertices
// and eliminate them.  This will reduce partially buried regions and shrink
// their boundaries, allowing what is left to pair with a valid secondary
// mesh boundary.  However, this is relatively expensive to do...
template <class M>
const APlink& Meshing<M>::match_boundaries(
	const List<Boundary>& aboundaries,
	const List<Boundary>& pboundaries)
{
	// Avoidance of division by zero
	// Will be a problem if ever the dividend is more than 100 - very unlikely
	auto NEAR0 = [](double l) -> double { return (l==0?(100*DBL_MIN):l); };

	// Match the two sets of boundaries
	// First set a value on each pair
	// Value is inverse to distance and stronger for normal alignment
	APD apd(aboundaries.size());
	Uint abi, pbi;
	for (abi = 0; abi < aboundaries.size(); abi++) {
		const auto& ab = aboundaries[abi].blist;
		const Vector& f = aboundaries[abi].centroid;
		auto& apdb = apd[abi];
		apdb.resize(pboundaries.size());
		for (pbi = 0; pbi < pboundaries.size(); pbi++) {
			const auto& pb = pboundaries[pbi].blist;
			Vector v = f - pboundaries[pbi].centroid;
			double lenv = NEAR0(v.length2());
			double ad = 0.0, pd = 0.0;
			for (const auto& bn: ab)
				ad += (v.dot(vertex(*additional, bn.inner).get_normal()))/lenv;

			for (const auto& bn: pb)
				pd += (v.dot(vertex(*primary, bn.inner).get_normal()))/lenv;

			apdb[pbi] = ad/ab.size() + pd/pb.size();
		}
	}
	// Then sort them out with uniqueness for additional mesh region pairings
	for (int di = 0; di < apd.size(); di++) {
		auto it = std::max_element(apd[di].cbegin(), apd[di].cend());
		pbi = it - apd[di].cbegin();
		setaplink(pbi, di, apd);
	}
	return aplink;
}

template <class M>
void Meshing<M>::convert_face(Uint v1, Uint v2, Uint v3)
{
	constexpr unsigned int perm[] = { 1, 2, 0 };
	Multiplet<3> ta = {v1, v2, v3};
	// Add the triangle to the interior list, if not seen before
	if (tmap.count(ta) == 0) {
		Uint ts = face_list.size();
		face_list.push_back(ta);

		// There are three triangle keys
		Uint ia[] = { 0, 1, 2 };
		for (int i = 0; i < 3; i++) {
			Multiplet<3> tak = {ta[ia[0]], ta[ia[1]], ta[ia[2]]};
			tmap.emplace(tak, ts);
			for (int j = 0; j < 3; j++) ia[j] = perm[ia[j]];
		}
	}
	
	// Add the edges, & map vertex indices to new edge indices
	Uint ti = tmap.at(ta);
	for (int i = 0; i < 3; i++) {
		Multiplet<2> ea = parr(std::minmax(ta[i],ta[perm[i]]));
		if (emap.count(ea) == 0) {
			emap.emplace(ea, edge_list.size());
			edge_list.push_back(ea);
		}
		Uint ei = emap.at(ea);
		for (int j = 0; j < 2; j++) {
			List<Uint>& ve = vemap[ea[j]];
			ve.push_back(ei);	// duplicates to be tolerated
		}
		if (efmap.count(ei) == 0) {
			Multiplet<2>& ef = efmap[ei];
			ef[0] = ti;
			ef[1] = std::numeric_limits<unsigned int>::max();
		}
		else {
			Multiplet<2>& ef = efmap.at(ei);
			if (ef[0] != ti) ef[1] = ti;
		}
	}
}

template <class M>
void Meshing<M>::convert_region(const M& m, const IdxOrdSet& region,
	IdxMap<>& vmap)			// Vertex indices old to new
{
	// Add the interior vertices, edges and triangles to the new lists
	// Add all the vertices first, as they are then referenced
	// Include boundary vertices, as these are known
	for (const auto& idx: region) {
		if (vmap.count(idx) == 0) {
			// Add the vertices, and map the old indices to new ones
			vmap.emplace(idx, vertex_list.size());
			const Vertex& vert = vertex(m, idx);
			vertex_list.push_back(Vertex(vert.get_vertex(), vert.get_normal()));
		}
	}
	// Now add the faces, and from them the edges
	for (const auto& idx: region) {
		// Add the faces, & map vertex indices to new triangle indices
		// Get each triangle associated with the vertex
		const List<Uint>& tlist = vertex_triangles(m, idx);
		for (const auto& t: tlist) {
			const Triangle& mt = triangle(m, t);

			// Check that the triangle is wholly interior and translate to new
			Multiplet<3> ta
				= {mt.get_v1_idx(), mt.get_v2_idx(), mt.get_v3_idx()};
			int i;
			for (i = 0; i < 3; i++) {
				if (vmap.count(ta[i]) == 0) break;
				ta[i] = vmap.at(ta[i]);
			}
			if (i < 3) continue;	// triangle crosses broundary
			convert_face(ta[0], ta[1], ta[2]);
		}
	}
}

template <class M>
void Meshing<M>::convert_boundary(const M& m, List<BoundaryNode>& boundary,
	IdxMap<>& vmap)			// Vertex indices old to new
{
	// This loop is unnecessary in current context, but necessary in general
	for (const auto& b: boundary) {
		if (vmap.count(b.inner) == 0) {
			// Add the vertices, and map the old indices to new ones
			vmap.emplace(b.inner, vertex_list.size());
			const Vertex& vert = vertex(m, b.inner);
			vertex_list.push_back(Vertex(vert.get_vertex(), vert.get_normal()));
		}
	}

	// Convert each boundary node into edges and faces
	Uint vbeg = vertex_list.size();
	Uint vend = vbeg + boundary.size() - 1;
	for (auto bnit = boundary.begin(); bnit != boundary.end(); bnit++) {
		bnit->idx = vertex_list.size();
		Uint next = (bnit->idx == vend ? vbeg : bnit->idx+1);
		vertex_list.push_back(bnit->vert);
		if (bnit->cavity) {
			auto bnit1 = bnit+1;
			if (bnit1 == boundary.end()) bnit1 = boundary.begin();
			convert_face(vmap.at(bnit->inner), vmap.at(bnit1->inner), next);
		}
		convert_face(bnit->idx, vmap.at(bnit->inner), next);
	}
}

// Utility class to facilitate boundary iterations
// Technically, this is a half reverse iterator as pvi goes backwards
template <typename T>
class HalfIterator {
	Uint ai, pi;			// half-iterators
	Uint alim, plim;		// half-iterator index limits

	List<T> *ab, *pb;		// the list of values being iterated
	Uint avi, pvi;			// value indexes
	T *abvio, *pbvio;		// old boundary values
	Uint astart, pstart;	// start value indexes
	Uint vin, apselect;		// next value index and whether a or p
	double step, d;			// to manage different a and p boundary lengths
	bool ended;				// used to enable a standard loop construct

public:
	// No use for default constructors
	HalfIterator() = delete;
	HalfIterator(const HalfIterator&) = delete;

	// Constructor from lists
	HalfIterator(List<T>& alist, Uint a_start, List<T>& plist, Uint p_start) {
		ab = &alist;
		pb = &plist;
		astart = a_start;
		pstart = p_start;

		avi = astart;
		pvi = pstart;

	    alim = ab->size();
		plim = pb->size();
		ai = 0;
		pi = plim-1;
		ended = (alim == 0 || plim == 0);
		if (ended) return;

		abvio = &ab->at(avi);
		pbvio = &pb->at(pvi);
		vin = avi;				// choose to start on a list
		apselect = 2;			// non-zero means a list, 2 indicates start
		step = double(ab->size())/pb->size();
		d = step;

		// Proposed change - needs to be before hit.update(), though...
		// ...so this would need to be built into the half iterator.
		// Initially: pamin, apmin DBL_MAX; amin ab.size(); pmin pb.size();
		// Adjust vertex positions to improve tube triangle shapes
		// Note that loop covers all vertices even after first iteration
		// if (apselect == 1) {  // first time around?
		//	Vertex& avert = vertex_list[ab[avi].idx];
		//	Vertex& pvert = vertex_list[pb[pvi].idx];
		//	// maintain min avi from pvi
		//	double aplen2 = (avert-pvert).length2();
		//	if (aplen2 < apmin) {
		//		apmin = aplen2;
		//		amin = avi;
		//	}
		//	// update avi from pvi and pvi from amin
		//	avert = (avert*2 + pvert)/3;
		//	if (amin < ab.size()) {
		//		pvert = (pvert*2 + vertex_list[ab[amin].idx])/3;
		//		amin = ab.size();
		//	}
		//	// reset pamin as a has changed
		//	pamin = std::numeric_limits<double>::max();
		// }
		// else if (apselect == 0) {
		//	-- converse of above => utility?
		// }

		// Need an initial update to complete this
		if (d > 1.0) {
			apselect = 1;
			vin = (++ai >= alim ? astart : ((ai + astart)%alim));
			d -= 1.0;
		}
		else {
			apselect = 0;
			vin = (pi >= plim ? pstart : ((pi-- + pstart)%plim)); //unsigned
			d += step;
		}
	};

	// Move on, favouring the more numerous nodes
	void update(void) {
		// Update indices
		if (apselect) {
			avi = vin;
			abvio = &ab->at(vin);
		}
		else {
			pvi = vin;
			pbvio = &pb->at(vin);
		}
		ended = (ai >= alim && pi >= plim);
		if (d > 1.0) {
			apselect = 1;
			vin = (++ai >= alim ? astart : ((ai + astart)%alim));
			d -= 1.0;
		}
		else {
			apselect = 0;
			vin = (pi >= plim ? pstart : ((pi-- + pstart)%plim)); //unsigned
			d += step;
		}
	};

	// Check for end condition, which is in terms of the driving iterators
	bool end(void) { return ended; };

	// Iterator value access
	T& a(void) { return ab->at(avi); };	// current a value
	T& ao(void) { return *abvio; };		// old a value
	T& p(void) { return pb->at(pvi); };	// current p value
	T& po(void) { return *pbvio; };		// old p value
	// next and current value, a or p - the value that changes
	T& n(void) {return (apselect == 0 ? pb->at(vin) : ab->at(vin)); };
	T& c(void) {return (apselect == 0 ? pb->at(pvi) : ab->at(avi)); };

};

// Utility for kludge at end of list_boundaries() --
static void edge2go(Uint e2g, const IdxMap<List<Uint>>& vemap,
	const IdxListN<2>& elist, IdxOrdSet& v2go)
{
	for (int i = 0; i < 2; i++) {
		const auto& vidx = elist[e2g][i];
		if (v2go.count(vidx) == 0) {
			v2go.emplace(vidx);
			for (Uint e2g2: vemap.at(vidx)) edge2go(e2g2, vemap, elist, v2go);
		}
	}
}

template <class M>
void Meshing<M>::link_boundaries(List<Boundary>& aboundaries,
	List<Boundary>& pboundaries)
{
	// Build node-pairings for each paired set of boundaries
	Uint abi, pbi;
	for (auto& ap: aplink) {
		std::tie(pbi, abi) = ap;
		auto& ab = aboundaries[abi].blist;
		auto& pb = pboundaries[pbi].blist;

		// Find the starting point as the pair that minimises the sum of
		// squared distances when working around the boundary
		Uint avi = 0;
		Uint pvi = 0;
		double dmin = std::numeric_limits<double>::max();
		for (unsigned int pvin = 0; pvin < pb.size(); pvin++) {
			double d = 0;
			for (HalfIterator<BoundaryNode> hit(ab, avi, pb, pvin);
				!hit.end(); hit.update()) {
				const Vector& av = vertex_list[hit.a().idx];
				const Vector& pv = vertex_list[hit.p().idx];
				d += (av - pv).length2();
			}
			if (d < dmin) {
				dmin = d;
				pvi = pvin;
			}
		}

		// avi pairs with pvi: loop through boundary adding edges and faces
		Uint fei = edge_list.size(), ti;
		for (HalfIterator<BoundaryNode> hit(ab, avi, pb, pvi);
			 !hit.end(); hit.update()) {

			// Update the check maps with the anticipated new edges and faces
			Uint ei = edge_list.size();
			vemap[hit.ao().idx].push_back(ei);
			vemap[hit.po().idx].push_back(ei);
			// Add the new triangles to the new edges
			ti = face_list.size();
			Multiplet<2>& ef = efmap[ei];
			ef[0] = ti-1;	// Wrong for first edge, but corrected later
			ef[1] = ti;
			// And add the triangle to the map for the existing edge
			ei = emap.at(parr(std::minmax(hit.n().idx, hit.c().idx)));
			Multiplet<2>& ef2 = efmap.at(ei);
			if (ef2[0] != ti) ef2[1] = ti;	// check should be redundant!

			// Add previous edge vertices in low-high order
			edge_list.push_back(parr(std::minmax(hit.ao().idx, hit.po().idx)));

			// Define a face in anti-clockwise order
			Multiplet<3> face = {hit.a().idx, hit.p().idx, hit.n().idx};
			face_list.push_back(face);
		}
		Multiplet<2>& ef = efmap[fei];	// Correct the first edge
		ef[0] = ti;
	}

	// The lists are now complete, but there could be isolated regions
	// resulting from the crude approach taken (see Mesh::create_mesh2).
	// Simply detect and eliminate these by building a sorted list of vertices
	// to go and then reducing the three lists.
	Uint ei;
	Multiplet<2> ef;
	IdxOrdSet vert2go;
	for (const auto& efp: efmap) {
		std::tie(ei, ef) = efp;
		// Detect absent face for this edge
		if (ef[1] >= face_list.size()) edge2go(ei, vemap, edge_list, vert2go);
	}
	IdxOrdSet e2go, f2go;
	for (Uint vidx: vert2go) {
		// erase the vertex
		vertex_list.erase(vertex_list.begin()+vidx);

		// Note the edges and faces to go
		for (Uint i = 0; i < edge_list.size(); i++) {
			const auto& edge = edge_list[i];
			if (edge[0] == vidx || edge[1] == vidx) e2go.emplace(i);
		}
		for (Uint i = 0; i < face_list.size(); i++) {
			const auto& face = face_list[i];
			if (face[0] == vidx || face[1] == vidx || face[2] == vidx)
				f2go.emplace(i);
		}
	}
	// Adjust the edge and face vertex indices
	for (auto rvit=vert2go.crbegin(); rvit!=vert2go.crend(); ++rvit) {
		const auto& vidx = *rvit;
		for (auto& edge: edge_list) {
			for (auto& evidx: edge) if (evidx > vidx) evidx--;
		}
		for (auto& face: face_list) {
			for (auto& fvidx: face) if (fvidx > vidx) fvidx--;
		}
	}
	// Now erase the noted elements - slightly more efficient in reverse
	for (auto reit=e2go.crbegin(); reit!=e2go.crend(); ++reit)
		edge_list.erase(edge_list.begin()+*reit);
	for (auto rfit=f2go.crbegin(); rfit!=f2go.crend(); ++rfit)
		face_list.erase(face_list.begin()+*rfit);
}

// Does a line // to z-axis in +ve direction from pt pass through triangle tv?
// This local function and object makes it easier to test simple cases

// Determine if z-line from pt crosses triangle
// NB z-line is here generalised to any semi-infinite line.
// Note ptr may not be changed if the z-line never reaches the triangle plane
// and if ptr is provided the caller must check for repeat hits.
static bool pt_triangle2(const Vector& pt, const Vector* tv, const Vector& z,
						PtSet& hitVertices, Vector *ptr = nullptr)
{
	constexpr unsigned int perm[6] = { 1, 2, 2, 0, 0, 1 };
	static const Vector origin = Vector(0,0,0);
	int i;		// general index

	// Find intersection of z-line with triangle plane
	for (i = 0; i < 3; i++) if (tv[i] != pt) break;
	if (i >= 3) {
		std::cerr << "MeshInstance~pt_triangle found null triangle"
				  << std::endl;
		throw std::exception();
	}
	const Vector& v = tv[i];
	// normal to node patch no use here - get a normal to triangle
	Vector n = (tv[2]-tv[1]).cross(tv[1]-tv[0]);
	Vector pz;
	double zdotn = z.dot(n);
	Vector q, r;
	// find a point in the plane that lies on the z-line from pt
	if (zdotn == 0.0) {				// z-dir is in the triangle's plane
		if ((v-pt).dot(n) == 0.0) {	// pt is also in the plane
			pz = pt;				// so use pt as the point in the plane
			// If the caller asked for the nearest point on the triangle but the
			// z-direction was in the plane, need to calculate the nearest point
			if (ptr) {
				// For each edge, find the nearest point in z-dir
				double d, dmin = (tv[3]-pt).length2();
				*ptr = tv[3];	// Just an initial point to be improved upon
				for (i = 0; i < 3; i++) {
					// Reduce the problem to 2 or 1 dimensions
					q = (tv[perm[2*i]] - tv[i]).normalised();  // not 0,0,0
					double zdotq = z.dot(q);
					r = (z-q*zdotq);		// orthogonal to q
					if (r == origin) {
						// Simple linear problem along q from tv[i]
						pz = tv[i];	// other vertices consided in loop
					}
					else {
						r = r.normalised();
						// 2-D problem, line1 = tv[i] + a q;  line2 = pt + b r
						// Intersect at:
						pz = pt + z*r.dot(tv[i]-pt)/sqrt(1-(zdotq)*(zdotq));
					}
					d = (pz-pt).length2();
					if (d < dmin) {
						dmin = d;
						*ptr = pz;
					}
				}
			}
		}
		else
			return false;			// cannot cross triangle if pt not in plane
	}
	else {							// z-dir is not in the plane
		double dz = ((v-pt).dot(n) / zdotn);
		if (dz < 0) return false;	// (z-line is _semi_-infinite)
		pz = dz * z + pt;
		if (ptr) *ptr = pz;
	}

	// Is pz in any of the half-planes that define the triangle?
	// In addition, for any triangle that is crossed at a vertex or edge,
	// there will be other triangles that are crossed at the same place.
	// This needs detecting to avoid multiple counts.  Although less than
	// half the triangles at this stage will be crossed, it is efficient
	// to carry out the vertex/edge checks here rather than later.
	double dot;
	double onedge = -1.0;
	for (i = 0; i < 3; i++) {
		// Get unit vector along an edge and vectors from opposite vertex and pz
		n = (tv[perm[2*i]]-tv[perm[2*i+1]]).normalised();
		q = tv[perm[2*i]]-tv[i];
		r = tv[perm[2*i]]-pz;

		// Check if pz is actually on the edge (or vertex)
		dot = q.dot(r);	// on edge if parallel, but can be antiparallel...
		if (dot >= 0.0 && dot*dot == q.length2()*r.length2() && onedge != 0.0)
			onedge = dot;

		// Obtain the normals to the selected edge from tv[i] and pt
		q = q - q.dot(n) * n;
		r = r - r.dot(n) * n;
#ifdef LOCAL_TEST // DEBUG
		std::cout << "MeshInstance~pt_triangle pt=" << pt << " pz=" << pz
				<< " T=" << tv[0] << "," << tv[1] << "," << tv[2]
				<< " n=" << n << " q=" << q << " r=" << r << " s=" << q.dot(r)
				<< std::endl;
#endif // DEBUG
		// If these normals are opposed, then pz is outside triangle
		if (q.dot(r) < 0) break;
	}
	// If loop did not complete then pt is outside triangle
	if (i < 3) return false;

	// We have a crossing, but if on an edge, then need to check for repeats.
	// However, that is up to the caller if they asked for the nearest hit.
	if (ptr || onedge < 0.0) return true;	// Strictly inside triangle

	// Check for repeat hits
	if (onedge == 0.0) {	// Hit a vertex!  Immensely unlikely!
		// This is actually quite annoying.  The only way to avoid repeats
		// is to maintain a list of hits;  but that list must be reset
		// for each caller call, so we have to provide a way to indicate that.
		// Since this function is local to the module for the one purpose,
		// (which makes testing it easier) it makes sense to let pt_is_internal
		// reset the list itself...
#if 0 // This check is now redundant, as edges are included here
		// First, though, check we did actually hit a vertex - precision?
		for (i = 0; i < 3; i++) if (tv[i] == pz) break;
		if (i >= 3) {
			std::cout << "MeshInstance::pt_is_internal confused! "
					  << pz << " " << tv[0] << " " << tv[1] << " " << tv[2]
					  << std::endl;
			std::cerr << "MeshInstance::pt_is_internal confused! Carrying on "
					  << std::endl;
			return true;
		}
#endif
std::cout << "MeshInstance~pt_triangle Way-hey! Hit a vertex!! " <<  pt << ": " << pz << "=" << tv[i] << std::endl;
#if 1 // edges included
	}
	else
std::cout << "MeshInstance~pt_triangle Hey! Hit an edge! " << pt << std::endl;
#endif // edges included
		if (hitVertices.count(pz) > 0) return false;  // Seen it already
std::cout << "MeshInstance~pt_triangle new hit" << std::endl;
		hitVertices.insert(pz);
		return true;
#if 0 // This old version is not thread-safe; edges handled with vertices
	}
std::cout << "MeshInstance~pt_triangle Hey! Hit an edge! " << pt << std::endl;
	// Not a vertex, so strictly an edge - this is easy!  Always hit twice...
	static bool flip = false;
	flip = !flip;	// This trick is only an issue if caller needs to identify
	return flip;	// each triangle that is hit;  this just ensures half count.
#endif
}

// This is the efficient version for just testing internality
template<class M>
bool Meshing<M>::pt_triangle(const Vector& pt, const Vector* tv,
							int zdir, PtSet& hitVertices)
{
	const Vector z = Vector(0,0,zdir);	// z-line used as indexing is neater
	int i;		// general index

	// Check straddles
	for (i = 0; i < 2; i++) {
		auto il = {tv[0](i), tv[1](i), tv[2](i)};
		if (pt(i) < std::min<double>(il) || pt(i) > std::max<double>(il)) break;
	}
	if (i < 2) return false;

	return pt_triangle2(pt, tv, z, hitVertices);
}

// General version to get the nearest node patch in a given direction
// See pt_triangle2 for parameters
template<class M>
bool Meshing<M>::pt_triangle(const Vector& pt, const Vector* tv,
				const Vector& z, PtSet& hitVertices, Vector *ptr)
{
	constexpr unsigned int perm[6] = { 1, 2, 2, 0, 0, 1 };
	int i;		// general index

	// Although more complicated than the version above, it is still much
	// more efficient to check for straddling than to skip it.
	// Find two normals to z in plane including origin
	for (i = 0; i < 3; i++) if (z(i) != 0.0) break;
	if (i >= 3) {
		// Trick question - null z
		std::cerr << "MeshInstance pt_triangle(pt,tv,z) null z" << std::endl;
		std::cout << "MeshInstance pt_triangle(pt,tv,z) null z" << std::endl;
		throw std::exception();
	}
	double znc[3];
	znc[i] = 0.0;
	znc[perm[2*i]] = z(perm[2*i+1]);
	znc[perm[2*i+1]] = -z(perm[2*i]);
	if (znc[0] == 0.0 && znc[1] == 0.0 && znc[2] == 0.0) {
		// Only one coordinate was non-zero
		znc[perm[2*i]] = 1.0;
		znc[perm[2*i+1]] = 0.0;
	}

	Vector zn[2];
	zn[0] = Vector(znc[0], znc[1], znc[2]);
	zn[1] = z.cross(zn[0]);

	// Check straddles - now very similar to the simple case above
	for (i = 0; i < 2; i++) {
		auto il = {zn[i].dot(tv[0]), zn[i].dot(tv[1]), zn[i].dot(tv[2])};
		double ptn = zn[i].dot(pt);
		if (ptn < std::min<double>(il) || ptn > std::max<double>(il)) break;
	}
	if (i < 2) return false;

	return pt_triangle2(pt, tv, z, hitVertices, ptr);
}

// Non-thread-safe default
template<class M>
PtSet Meshing<M>::pt_hitVertices;

template<class M>
bool Meshing<M>::pt_is_internal(const M& m, const Vector& pt,
					PtSet& hitVertices, const Vector *dir, Vector* nearest)
{
    // if the point is outside of the maximum radius then it cannot be
    // inside the mesh instance
    const Vector& xyz_offset = centre(m);
    if ((pt - xyz_offset).length() >= m.get_radius()) return false;

	// Find the triangles for this mesh and
	// count number of times a z-line from pt crosses a mesh triangle
	Vector tv[3];	// Triangle - copies because const pointers are a mess
	int cc = 0;	// crossings count
	const int zdir = 2*(pt.z - xyz_offset.z > 0.0)-1;	// z-direction to use
	hitVertices.clear();	// See interface to pt_triangle above
	auto& tlist = triangles(m);
	Vector retnearest;
	double dmin = std::numeric_limits<double>::max();
	bool hit, earlyExit = false;
	if (nearest && pt != *nearest) {
		earlyExit = true;
		dmin = (pt-*nearest).length2();
	}
    for (const auto t: tlist) {
		unsigned int vi[] = {t.get_v1_idx(), t.get_v2_idx(), t.get_v3_idx()};
		// This relies on vertices and patches having the same order...
		for (int i = 0; i < 3; i++) tv[i] = patch(m, vi[i]).get_node();

		// Does z-line does cross the mesh at this triangle?
		if (dir)
			hit = pt_triangle(pt, tv, *dir, hitVertices, nearest);
		else
			hit = pt_triangle(pt, tv, zdir, hitVertices);
		cc += hit;

		// If nearest was supplied, maintain the nearest value to pt
		if (nearest && hit) {
			double d =  (pt-*nearest).length2();
			if (d < dmin) {
				dmin = d;
				retnearest = *nearest;
				if (earlyExit) {
					cc = 1;
					break;
				}
			}
		}
	}  // each triangle in mesh
	if (nearest) *nearest = retnearest;
	return ((cc%2) == 1);  // pt is internal if odd number of crossings
}

// And now tell the compiler what to generate
template class Meshing<Mesh>;
template class Meshing<MeshInstance>;
