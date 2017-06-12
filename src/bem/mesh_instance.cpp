/*      Author: david fallaize    Created on: 27 Jul 2010 */
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file mesh_instance.cpp
 * \brief This module implements the mesh instance and list classes.
 *
 */

#include "mesh_instance.h"
#include <iostream>
#include <iomanip>
#include <unordered_map>


MeshInstance::MeshInstance(
	unsigned int mesh_lib_id,
	unsigned int mesh_instance_id,
	MeshList& mesh_library,
	const Vector& offset,
	const Quaternion& rot,
	double protein_dielectric,
	double solvent_dielectric,
	unsigned int num_quad_points,
	unsigned int num_qual_points,
	bool _silent
	) : instance_id(mesh_instance_id),
		silent(_silent),
		xyz_offset(offset),
		rotation(rot),
		quad_points_per_triangle(num_quad_points),
		qual_points_per_triangle(num_qual_points)
#ifdef USING_MAXD  // maxd no longer required
		,maxd(std::make_unique<double[]>(2))
#endif // maxd
{
    // rather essential- copy the shared_ptr to the underlying mesh type!
    assert(mesh_lib_id <= mesh_library.size());
	mesh_ptr = boost::shared_ptr<ListedMesh>(mesh_library[mesh_lib_id]);
    
    init();
    set_dielectrics(protein_dielectric, solvent_dielectric);
}

void MeshInstance::init()
{
    const Mesh& mesh = *mesh_ptr;

    // set BasicNodePatches
    patches.clear();
    const std::vector<BasicNodePatch>& library_node_patches
		= mesh.get_node_patches();
    for (std::vector<BasicNodePatch>::const_iterator
			it=library_node_patches.cbegin(), end=library_node_patches.cend();
		 it != end; ++it)
    {
        // copy the node patches from the reference mesh
        // create the NodePatch but immediately store it as a pointer to the
        // base class
        boost::shared_ptr<BasicNodePatch> ptr(new NodePatch(*it, *this)); 
        patches.push_back(ptr);
    }

    // set Charges
    charges.clear();
    const std::vector<Charge>& library_charges = mesh.get_charges();
    for (std::vector<Charge>::const_iterator
			it=library_charges.cbegin(), end=library_charges.cend();
		 it != end; ++it)
    {
        boost::shared_ptr<Charge> ptr(
			new Charge(*it, mesh.get_centre(), rotation, xyz_offset)
		);
        charges.push_back(ptr);
    }
    
    // set radius
    radius = mesh.get_radius();

    // create a simple octree of the mesh instance (for collision / grid checks)
	reset_mesh_tree();
}

void MeshInstance::reset_mesh_tree() {
    mesh_tree.reset(new Octree< Node<BasicNodePatch>, BasicNodePatch >
						(10, xyz_offset, radius*4.) );
    for (PatchList::iterator it=patches.begin(), end=patches.end();
		 it != end; ++it)
    {
        mesh_tree->add(*it);
    }
}

// Does a line // to z-axis in +ve direction from pt pass through triangle tv?
// This local function and object makes it easier to test simple cases
// An unordered map is OTT but why not?  See onedge handling in pt_triangle...
// Need two unary function objects:
struct hashFunc{
    size_t operator()(const Vector &v) const{
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
// Map of hit vertices - must be reset by caller of pt_triangle
static std::unordered_map<Vector, bool, hashFunc, equalsFunc> hitVertices;

// Determine if z-line from pt crosses triangle
static bool pt_triangle(const Vector& pt, const Vector* tv, int zdir) {
	constexpr unsigned int perm[6] = { 1, 2, 2, 0, 0, 1 };
	const Vector z = Vector(0,0,zdir);	// z-line used as indexing is neater
	int i;		// general index

	// Check straddles
	for (i = 0; i < 2; i++) {
		auto il = {tv[0](i),tv[1](i),tv[2](i)};
		if (pt(i) < std::min<double>(il) || pt(i) > std::max<double>(il)) break;
	}
	if (i < 2) return false;

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
	if (z.dot(n) == 0.0) {
		if ((v-pt).dot(n) == 0.0)
			pz = pt;
		else
			return false;
	}
	else {
		double dz = ((v-pt).dot(n) / z.dot(n));
		if (dz < 0) return false;  // z-line is _semi_-infinite
		pz = dz * z + pt;
	}

	// Is pz in any of the half-planes that define the triangle?
	// In addition, for any triangle that is crossed at a vertex or edge,
	// there will be other triangles that are crossed at the same place.
	// This needs detecting to avoid multiple counts.  Although less than
	// half the triangles at this stage will be crossed, it is efficient
	// to carry out the vertex/edge checks here rather than later.
	Vector q, r;
	double dot;
	double onedge = -1.0;
	for (i = 0; i < 3; i++) {
		// Get unit vector along an edge and vectors from opposite vertex and pz
		n = (tv[perm[2*i]]-tv[perm[2*i+1]]).normalised();
		q = tv[perm[2*i]]-tv[i];
		r = tv[perm[2*i]]-pz;

		// Check if pz is actually on the edge (or vertex)
		dot = q.dot(r);	// on edge if parallel, but can be antiparallel...
		if (dot >= 0.0 && dot*dot == q.length2() * r.length2()) onedge = dot;

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
	if (i < 3) return false;
	if (onedge < 0.0) return true;

	// We have a crossing, but it was on a shared edge or vertex...
	if (onedge == 0.0) {	// Hit a vertex!  Immensely unlikely!
		// This is actually quite annoying.  The only way to avoid repeats
		// is to maintain a list of hits;  but that list must be reset
		// for each caller call, so we have to provide a way to indicate that.
		// Since this function is local to the module for the one purpose,
		// (which makes testing it easier) it makes sense to let pt_is_internal
		// reset the list itself...
std::cout << "MeshInstance~pt_triangle Way-hey! Hit a vertex!!" << std::endl;
		if (hitVertices.count(v) > 0) return false;  // Seen it already
		hitVertices[v] = true;	// Dummy value, key is the key
		return true;
	}
std::cout << "MeshInstance~pt_triangle Hey! Hit an edge!" << std::endl;
	// Not a vertex, so strictly an edge - this is easy!  Always hit twice...
	static bool flip = false;
	flip = !flip;	// This trick is only an issue if caller needs to identify
	return flip;	// each triangle that is hit;  this just ensures half count.
}

bool MeshInstance::pt_is_internal(const Vector& pt) const {
// Alas, the simple algorithm deleted below does not work for meshes - the
// fix would require finding the nearest point on the triangle plane, then
// the nearest point to that on each triangle side.  The resulting effort
// is similar to that of the replacement algorithm but less easily made
// as efficient.
    // if the point is outside of the maximum radius then it cannot be
    // inside the mesh instance
    if ((pt - xyz_offset).length() >= radius) return false;

#if 1 // TESTING BOTH METHODS else delete, re-enable __DELETED__ & do returns
	bool retval;
	while (true) {
//#ifdef __DELETED__
    // ok so it's within the maximum radius, still not necessarily
    // internal -- find the nearest node patch and compare this 
    // point to the normal vector of the patch
    const BasicNodePatch& nearest_np = mesh_tree->get_nearest(pt);
    if ((pt - nearest_np).dot(nearest_np.get_normal()) > 0) {retval= false; break;}

    // if get here then the above test must indicate that the point
    // is on the internal side of the nearest node patch and is 
    // therefore internal to the mesh instance.
    retval= true; break;
//#else   // __DELETED __
	}

	Vector tv[3];	// Triangle - must be copies because frame changes
	int i, c;		// General and coordinate indices

#ifdef USING_MAXD  // maxd no longer required - sadly!
	// Lambda to test if triangle is too far from z-line
	auto tmaxd = [&]() mutable noexcept -> bool  // mutable to match smaxd
		{ return (tv[i](c) > pt(c)+maxd[c] || tv[i](c) < pt(c)-maxd[c]); };
	// Lambda to set maxd
	auto smaxd = [&]() mutable noexcept -> bool
		{ double d = abs(tv[i](c)-tv[(i+1)%3](c));
		  if (d > maxd[c]) maxd[c] = d; return false; };
	// Choice of lambda is made once
	std::function<bool()> amaxd[] = { tmaxd, smaxd };
	auto lmaxd = amaxd[maxd[0] <= 0.0];
#endif // maxd

	// Find the triangles for this mesh and
	// count number of times a z-line from pt crosses a mesh triangle
	int cc = 0;	// crossings count
	const int zdir = 2*(pt.z - xyz_offset.z > 0.0)-1;	// z-direction to use
	hitVertices.clear();	// See interface to pt_triangle above
    for (size_t np_ctr=0; np_ctr < patches.size(); ++np_ctr) {
        const NodePatch& np = dynamic_cast<NodePatch&>(*(patches[np_ctr]));
        
		tv[0] = np;
        Vector local_centre = mesh_ptr->get_centre();
        std::vector<PointNormal> edge_points;
        np.get_edge_points(edge_points);
        for (std::vector<PointNormal>::const_iterator
				it=edge_points.cbegin(), next_it, end=edge_points.cend();
			 it!=end; ++it)
        {
            next_it = it+1;
            if (next_it == end) next_it = edge_points.begin();
            
            tv[1] = it->pt();
            tv[1].change_coordinate_frame(local_centre, rotation, xyz_offset);
            tv[2] = next_it->pt();
            tv[2].change_coordinate_frame(local_centre, rotation, xyz_offset);

#ifdef USING_MAXD  // maxd no longer required
			// Depending on whether maxd is initialised, either:
			//	 Record the maximum x,y difference up to this triangle; or
			//   Skip any triangles too far from the line for it to cross them.
			for (c = 0; c < 2; c++) {	// coordinate number 0=x, 1=y
				for (i = 0; i < 3; i++)	// triangle vertex
					if (lmaxd()) break; // true => skip
				if (i < 3) break;
			}
			if (c < 2) continue;
#endif

			// Does z-line does cross the mesh at this triangle?
			cc += pt_triangle(pt, tv, zdir);
		}  // each triangle in node patch
	}  // each node patch
if (retval != ((cc%2) == 1)) { // Shout and record
std::cerr << "MeshInstance::pt_is_internal instance " << instance_id
		<< " new=" << cc << " " << ((cc%2) == 1)
		<< " old=" << retval << std::endl;
std::cout << "MeshInstance::pt_is_internal instance " << instance_id
		<< " new=" << cc << " " << ((cc%2) == 1)
		<< " old=" << retval << std::endl;
}
	return ((cc%2) == 1);  // pt is internal if odd number of crossings
#endif  // __DELETED__
}

double MeshInstance::get_potential_at_internal_pt(const Vector& pt) const
{
    assert(this->pt_is_internal(pt) == true);

    // Assuming the point is internal to the mesh, the potential
    // is the sum of the coulomb terms from the charges and the
    // induced charge/dipole contribs from the surface.
    double pot=0.0;

    if (charge_fmm_tree.get() == nullptr) {
        charge_fmm_tree.reset
			(new fmm::FMM_Octree_6FIG_ACCURACY(10000, xyz_offset, radius * 4.));
        for (std::vector<boost::shared_ptr<Charge> >::iterator
				charge_it=charges.begin(), charge_end=charges.end();
             charge_it != charge_end; ++charge_it)
        {
            charge_fmm_tree->add(*charge_it);
        }
		// solve for kappa=0, i.e. non-screened Coulomb potential
        charge_fmm_tree->solve(0.0);
    }
    pot = charge_fmm_tree->calculate_potential(pt);
    pot *= ONE_OVER_4PI / Dprotein;

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
			np_it=patches.begin(), np_end=patches.end(); 
         np_it != np_end; ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double Gpt_val = Gpt(pt, np);
        double dGpt_val = dGpt_dn(pt, np, n);
        accum += np.get_bezier_area()
					* (-Gpt_val*np.h*np.dielectric_ratio + dGpt_val*np.f);
    }

    return pot+accum;
}

double MeshInstance::get_potential_at_external_pt(
	const Vector& pt,
	double kappa) const
{
    assert(this->pt_is_internal(pt) == false);

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
			np_it=patches.cbegin(), np_end=patches.cend(); 
         np_it != np_end; ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double upt_val = upt(pt, np, kappa);
        double dupt_val = dupt_dn(pt, np, n, kappa);
        accum += np.get_bezier_area() * (dupt_val*np.f - upt_val*np.h);
    }
    
    return accum;
}

int MeshInstance::move(const Vector& translate, const Quaternion& rotate) {
//#define __LOCAL_MOVES__ 1 -- defined in node_patch.h and TODO remove #ifdef
#ifndef __LOCAL_MOVES__
	std::cerr << "MeshInstance::move() has not been implemented" << std::endl;
	// change the coordinate frame for the mesh and charges
	// and reset the octree values
	// Update current position and rotation

#else // __LOCAL_MOVES__
std::cout << "MeshInstance::move from " << xyz_offset << " to " << translate << std::endl;
	unsigned int npctr = 0;
    for (std::vector<boost::shared_ptr<BasicNodePatch>>::iterator
			it = patches.begin(), end = patches.end();
		 it != end; ++it)
    {
		(static_cast<NodePatch&>(**it)).change_coordinate_frame
										(xyz_offset, rotate, translate);
#if 1
		// Reset the potential and derivative
		// to the equivalent ref mesh node_patch values
		const BasicNodePatch& np = mesh_ptr->get_node_patch(npctr);
		(**it).f = np.f;
		(**it).h = np.h;
		(**it).energy_coefficient_f = np.energy_coefficient_f;
		(**it).energy_coefficient_h = np.energy_coefficient_h;
		(**it).force_coefficient_f = np.force_coefficient_f;
		(**it).force_coefficient_h = np.force_coefficient_h;
		(**it).gc = np.gc;
//TODO NB Other values to reset?  quads
		npctr++;
#endif // if 0
    }

#ifdef USING_MAXD  // maxd no longer required
	// Reset maximum dimensions of triangles - will be recalculated when needed
	for (int i = 0; i < 2; i++) maxd[i] = 0.0;
#endif // maxd

    // set Charges
    for (std::vector<boost::shared_ptr<Charge>>::iterator
			it = charges.begin(), end = charges.end();
		 it != end; ++it)
    {
		// Rotation and translation is being applied to a mesh instance
		// rather than a mesh reference.
		(**it).Vector::change_coordinate_frame(xyz_offset, rotate, translate);
    }
    
	// Update object copy of position and rotation - before resetting tree!!
    xyz_offset = translate;
    rotation *= rotate;

    // create a simple octree of the mesh instance (for collision / grid checks)
	reset_mesh_tree();
	
#endif // __LOCAL_MOVES__
}

double MeshInstance::calculate_energy(
	double kappa,
	double fvals[],
	double hvals[]) const
{
    // Dsolvent and Dprotein should have been set via set_dielectrics()
    double E = mesh_ptr->calculate_energy(kappa, Dprotein, Dsolvent,
	                                      fvals, hvals);
    std::cout << "Energy for mesh " << instance_id << " (lib_id="
	          << mesh_ptr->get_id() << ") = " << std::setprecision(10)
	          << E << std::endl;
    return E;
}

Vector MeshInstance::calculate_force(
	double kappa,
	double fvals[],
	double hvals[]) const
{
    // Dsolvent and Dprotein should have been set via set_dielectrics()
    KahanVector MST_ext,MST_int,dbf,ionic;
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent,
								 fvals, hvals, MST_ext, MST_int, dbf, ionic);
    Vector force(*MST_ext);
    force.apply_rotation(rotation);
    return force;
}

void MeshInstance::calculate_forces(
	double kappa,
	double fvals[],
	double hvals[],
	KahanVector& qE,
	KahanVector& MST_ext,
	KahanVector& MST_int,
	KahanVector& dbf,
	KahanVector& ionic) const
{
    // Dsolvent and Dprotein should have been via the set_dielectrics()
    qE = mesh_ptr->calculate_qe_force(Dprotein, Dsolvent, fvals, hvals);
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent,
								 fvals, hvals, MST_ext, MST_int, dbf, ionic);
    
    (qE).apply_rotation(rotation);
    (MST_int).apply_rotation(rotation);
    (MST_ext).apply_rotation(rotation);
    (ionic).apply_rotation(rotation);
    (dbf).apply_rotation(rotation);
}

void MeshInstance::kinemage_fh_vals(
	double fscale,
	double hscale,
	int num_colours,
	std::ostringstream& buf_f,
	std::ostringstream& buf_h) const
{
    for (size_t np_ctr=0; np_ctr < patches.size(); ++np_ctr) {
        const NodePatch& np = dynamic_cast<NodePatch&>(*(patches[np_ctr]));
        
        std::string f_name;
        std::string h_name;
		auto NEAR0 = [](double l) -> double { return (l==0?DBL_MIN:l); };

        // figure out the colour name corresponding to the fval
        {  // New scope
            std::stringstream s;
            int f_idx = static_cast<int>(round(
							static_cast<double>(num_colours)*np.f / NEAR0(fscale)
						) );
            if (f_idx < 1){
                f_idx = abs(f_idx);
                int fcol = (f_idx <= num_colours ? f_idx : num_colours+1);
                s << "red" << "_" << fcol;
            } else {
                int fcol = (f_idx <= num_colours ? f_idx : num_colours+1);
                s << "blue" << "_" << fcol;
            }
            f_name = s.str();
        }  // End scope

        // figure out the colour name corresponding to the hval
        {  // New scope
            std::stringstream s;
            int h_idx = static_cast<int>(round(
							static_cast<double>(num_colours)*np.h / NEAR0(hscale)
						) );
            if (h_idx < 1){
                h_idx = abs(h_idx);
                int hcol = h_idx <= num_colours ? h_idx : num_colours+1;
                s << "red" << "_" << hcol;
            } else {
                int hcol = h_idx <= num_colours ? h_idx : num_colours+1;
                s << "blue" << "_" << hcol;
            }
            h_name = s.str();
        }  // End scope

        Vector local_centre = mesh_ptr->get_centre();
        std::vector<PointNormal> edge_points;
        np.get_edge_points(edge_points);
        for (std::vector<PointNormal>::const_iterator
				it=edge_points.cbegin(), next_it, end=edge_points.cend();
			 it!=end; ++it)
        {
            next_it = it+1;
            if (next_it == end) next_it = edge_points.begin();
            
            Vector here = it->pt();
            here.change_coordinate_frame(local_centre, rotation, xyz_offset);
            Vector next = next_it->pt();
            next.change_coordinate_frame(local_centre, rotation, xyz_offset);

            // kinemage a quadrilateral with the fh values
            buf_f << "X " << f_name << " "
					    << np.x << " " << np.y << " " << np.z << " "
                        << f_name << " "
                        << here.x << " " << here.y << " " << here.z << " "
                        << f_name << " "
                        << next.x << " " << next.y << " " << next.z << "\n";

            // kinemage a quadrilateral with the fh values
            buf_h << "X "  << h_name << " "
					    << np.x << " " << np.y << " " << np.z << " "
                        << h_name << " "
                        << here.x << " " << here.y << " " << here.z << " "
                        << h_name << " "
                        << next.x << " " << next.y << " " << next.z << "\n";
        }
    }
}

void MeshInstance::set_unique_patch_id(unsigned int patch_ctr)
{
	for (PatchList::iterator nit=patches.begin(), nend=patches.end();
		 nit != nend; ++nit)
	{
		(**nit).set_idx(patch_ctr++);
	}
}

bool MeshInstance::init_fh_vals(
	double *x,
	unsigned int& xctr,
	unsigned int offset,
	bool preconditioned) const
{
	const Mesh& mesh = *mesh_ptr;
	const std::vector<BasicNodePatch>& patches = mesh.get_node_patches();
	for (unsigned int ctr=0; ctr < patches.size(); ++ctr) {
		double f_preset = patches[ctr].f;
		double h_preset = patches[ctr].h;
		preconditioned = preconditioned || f_preset != 0.0 || h_preset != 0.0;
		x[xctr] = f_preset;
		x[xctr+offset] = h_preset;
		xctr++;
	}
	return preconditioned;
}

//NB where do f_lhs etc really belong?
unsigned int MeshInstance::reset_fh_vals(
	const boost::shared_array<double>& f_lhs, 
	const boost::shared_array<double>& h_lhs,
	unsigned int start_ctr)
{
	unsigned int inc_ctr = start_ctr;
	Mesh& mesh = *mesh_ptr;
	std::vector<BasicNodePatch>& patches = mesh.get_node_patches();
	for (unsigned int ctr=0; ctr < patches.size(); ++ctr) {
		BasicNodePatch& np = patches[ctr];
		np.f = f_lhs[inc_ctr];
		np.h = h_lhs[inc_ctr];
		++inc_ctr;
	}
	return inc_ctr;
}


// MeshInstanceList

// Add a new mesh instance to the list
boost::shared_ptr<MeshInstance>
MeshInstanceList::add(
	unsigned int mesh_lib_id,
	unsigned int mesh_instance_id,
	const Vector& offset,
	const Quaternion& rotation,
	double Dprotein, double Dsolvent,
	unsigned int num_quad_points,
	unsigned int num_qual_points,
	bool _silent)
{
	boost::shared_ptr<MeshInstance> ptr(
		new MeshInstance(mesh_lib_id, mesh_instance_id, library,
							offset, rotation, Dprotein, Dsolvent,
							num_quad_points, num_qual_points, _silent)
	);

//NB Comments inherited from BEEP::insert_mesh_instance
// minimum separation 5A between node patches
//     const double minimum_separation = 5.;
//     
//     // loop over the node patch points and ensure that they are not within   an angstrom
//     // of any other node patch already in the system
//     for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator      minst_it = begin(), minst_end=end(); minst_it != minst_end; ++    minst_it)
//     {
//         const MeshInstance& minst = **minst_it;
//         bool further_check_required = false;
//         if ((minst.xyz_offset - ptr->xyz_offset).length() <    (ptr->radius + minst.radius + minimum_separation)) {
//             further_check_required = true;
//         }
// 
//         if (further_check_required)
//         {
// 
//             throw std::exception();
//             
// /*            for (PatchList::const_iterator
//                 
//                 
//             for (PatchList::const_iterator np_it = minst.patches.begin(),    np_end=minst.patches.end(); np_it != np_end; np_it != np_end)
//             {
//                 const BasicNodePatch& np = **np_it;
//                 np
//             }*/
//             
//         }
//         
//     }

	// set unique id on each node patch
	unsigned int patch_ctr = get_total_patches();
	ptr->set_unique_patch_id(patch_ctr);

	push_back(ptr);
	return ptr;
}

// Move an instance in the list
// This exists because the simpler MeshInstance::move() is not working:
// the energy calculation for the moved instance is too high and it is not
// apparent why this should be.  This approach mimics clear and insert
// which appear to work individually.
// Note: cannot just copy construct the new MeshInstance for the same reason.
boost::shared_ptr<MeshInstance>
MeshInstanceList::move(
	unsigned int mesh_instance_id,
	const Vector& offset,
	const Quaternion& rotation,
	double Dprotein, double Dsolvent,	//TODO NB Passed because not accessible!
	bool _silent)						// Usage?? Can they be removed?
{
	// Validate the id
    if (mesh_instance_id >= size()) {
        std::cerr << "Bad value for mesh instance_id: " << mesh_instance_id
			<< " (" << size() << " mesh instances defined)" << std::endl;
        throw std::exception();
    }
	
	// Obtain require values from existing MeshInstance
	MeshInstance& m = *(this->at(mesh_instance_id));
	unsigned int mesh_lib_id
		= static_cast<const ListedMesh&>(m.get_ref_mesh()).get_id();
	unsigned int num_quad_points = m.get_quad_points_per_triangle();
	unsigned int num_qual_points = m.get_qual_points_per_triangle();
	unsigned int patch_ctr = (**(m.get_node_patches().begin())).get_idx();
	
	// Create a new instance
	boost::shared_ptr<MeshInstance> ptr(
		new MeshInstance(mesh_lib_id, mesh_instance_id, library,
							offset, rotation, Dprotein, Dsolvent,
							num_quad_points, num_qual_points, _silent)
	);

	// set unique id on each node patch
	ptr->set_unique_patch_id(patch_ctr);

	// Insert the new instance at the required point
	erase(begin()+mesh_instance_id);
	insert(begin()+mesh_instance_id, ptr);
	return ptr;
}


// INIT INITIAL VALUES to what was set in input fh files
bool MeshInstanceList::init_library_fh_vals(double *x, unsigned int offset)
{
    bool preconditioned = false;
    unsigned int xctr=0;
    for (MeshInstanceList::const_iterator mit=begin(), mend=end();
		 mit!=mend; ++mit)
    {
        const MeshInstance& minst = **mit;
		preconditioned = preconditioned
		               || minst.init_fh_vals(x, xctr, offset, preconditioned);
    }
	return preconditioned;
}

// Reset the f/h node patch values of the mesh definitions
void MeshInstanceList::reset_library_fh_vals(
	const boost::shared_array<double>& f_lhs,
	const boost::shared_array<double>& h_lhs)
{
    unsigned int total_ctr=0;
    for (std::vector< boost::shared_ptr<MeshInstance> >::iterator
			mit=begin(), mend=end();
		 mit != mend; ++mit)
    {
        MeshInstance& minst = **mit;
		total_ctr += minst.reset_fh_vals(f_lhs, h_lhs, total_ctr);
	}
}

size_t MeshInstanceList::get_total_patches() const
{
    unsigned int total_np=0;
    for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
             it=cbegin(), end=cend();
         it != end; ++it)
    {
        const MeshInstance& m = **it;
        total_np += m.get_num_node_patches();
    }
    return total_np;
}

