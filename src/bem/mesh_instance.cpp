/*      Author: david fallaize    Created on: 27 Jul 2010 */
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file mesh_instance.cpp
 * \brief This module implements the mesh instance and list classes.
 *
 */

#include "mesh_instance.h"
#include <iostream>
#include <iomanip>


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
{
    // rather essential- copy the shared_ptr to the underlying mesh type!
    assert(mesh_lib_id <= mesh_library.size());
#ifdef __DELETED__
	mesh_ptr = mesh_library[mesh_lib_id];
#else
	mesh_ptr = boost::shared_ptr<ListedMesh>(mesh_library[mesh_lib_id]);
#endif
    
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

bool MeshInstance::pt_is_internal(const Vector& pt) const
{
    // if the point is outside of the maximum radius then cannot be
    // inside the mesh instance
    if ((pt - xyz_offset).length() >= radius) return false;

    // ok so it's within the maximum radius, still not necessarily
    // internal -- find the nearest node patch and compare this 
    // point to the normal vector of the patch
    const BasicNodePatch& nearest_np = mesh_tree->get_nearest(pt);
    if ((pt - nearest_np).dot(nearest_np.get_normal()) > 0) return false;

    // if get here then the above test must indicate that the point
    // is on the internal side of the nearest node patch and is 
    // therefore internal to the mesh instance.
    return true;
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
	// and reset the octree
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
//TODO NB Other values to reset?
		npctr++;
#endif // if 0
    }

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
#ifndef __DELETED__
    std::cout << "Energy for mesh " << instance_id << " (lib_id="
	          << mesh_ptr->get_id() << ") = " << std::setprecision(10)
	          << E << std::endl;
#else
    std::cout << "Energy for mesh " << instance_id << " = "
			  << std::setprecision(10) << E << std::endl;
    return E;
#endif // __DELETED__
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

        // figure out the colour name corresponding to the fval
        {  // New scope
            std::stringstream s;
            int f_idx = static_cast<int>(round(
							static_cast<double>(num_colours)*np.f / fscale
						) );
            if (f_idx < 1){
                f_idx = abs(f_idx);
                int fcol = f_idx <= num_colours ? f_idx : num_colours+1;
                s << "red" << "_" << fcol;
            } else {
                int fcol = f_idx <= num_colours ? f_idx : num_colours+1;
                s << "blue" << "_" << fcol;
            }
            f_name = s.str();
        }  // End scope

        // figure out the colour name corresponding to the fval
        {  // New scope
            std::stringstream s;
            int h_idx = static_cast<int>(round(
							static_cast<double>(num_colours)*np.h / hscale
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
            buf_f << "X " << np.x << " " << np.y << " " << np.z << " "
                        << here.x << " " << here.y << " " << here.z << " "
                        << f_name << " "
                        << next.x << " " << next.y << " " << next.z << "\n";

            // kinemage a quadrilateral with the fh values
            buf_h << "X "  << np.x << " " << np.y << " " << np.z << " "
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

#ifndef __DELETED__
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
//     const double minimum_separation = 5.; //     
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
#endif // ! __DELETED__
