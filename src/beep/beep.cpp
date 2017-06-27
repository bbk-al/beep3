/* Author: david fallaize   Created on: 23 Jul 2010 */

/*!\file beep.cpp
 * \brief This module implements the BEEP class for boundary element
 *	electrostatic potential calculations.
 */
#include "beep.h"
#include "../bem/mesh.h"
#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_cblas.h>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include "../fmm/fmm_octree.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include "../bem/node_patch.h"
#include "../common/useful_clock.h"
#include "../bem/constants.h"
#include <boost/format.hpp>

BEEP::BEEP(const ConfigFile& config, 
           int cmdline_quad_points, 
           int cmdline_qual_points,
           int _cmdline_bem_nbsize,
           double _cmdline_kappa,
           bool planar) : 
	manually_set_bounding_cube(false), 
	beta0 (1e-10), 
	cmdline_bem_nbsize(_cmdline_bem_nbsize), 
	cmdline_quad_points_per_triangle(cmdline_quad_points),
	cmdline_qual_points_per_triangle(cmdline_qual_points), 
	cmdline_kappa(_cmdline_kappa),
	force_planar(planar),
#ifdef PREHYDROPHOBIC
	resolve(true),
	previd(1),	// Will be invalid on BEEP creation (as required)
	prevloc(Vector(0,0,0)),
	unrot(Quaternion(1,0,0,0)),
#endif // PREHYDROPHOBIC
	skipping_precalcs(false)
{
    init(config);
}

// This constructor is for python post-processing scripts
//NB Automatically calculates energies and forces, so single scenario
BEEP::BEEP(const std::string& config_filename, bool planar, bool read_fh) : 
    manually_set_bounding_cube(false), 
    beta0 (1e-10), 
    cmdline_bem_nbsize(-1), 
    cmdline_quad_points_per_triangle(-1), 
    cmdline_qual_points_per_triangle(-1), 
    cmdline_kappa(-1),
    force_planar(planar),
#ifdef PREHYDROPHOBIC
	resolve(true),
	previd(1),	// Will be invalid on BEEP creation (as required)
	prevloc(Vector(0,0,0)),
	unrot(Quaternion(1,0,0,0)),
#endif // PREHYDROPHOBIC
    skipping_precalcs(false)
{
    // parse xml config
    ConfigFile config(config_filename);
    
    // load the meshes and instantiate them as objects in space
    init(config);
    
    // read in fh files from results files and do post-processy stuff
    // (assuming read_fh == true)
    if (read_fh)
    {
        unsigned int num_patches = get_total_patches();
        const std::string& fh_filename = config.output_file; // stem of output files, which we will read
        
        f_lhs = boost::shared_array<double>(new double[num_patches]);
        h_lhs = boost::shared_array<double>(new double[num_patches]);
        
        for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
				mit=meshes.cbegin(), mend=meshes.cend();
			 mit != mend; ++mit)
        {
            MeshInstance& minst = **mit;
            
            // silent meshes don't write their fh values
            if (minst.isSilent()) continue;
            
            std::stringstream name_buf;
            name_buf << fh_filename << "." << minst.get_id();
            std::cout << "Reading from: " << name_buf.str() << std::endl;
            
            std::ifstream fh_input;
            fh_input.open(name_buf.str().c_str(), std::ios_base::in);    
            
            for (PatchList::iterator
					nit = minst.get_node_patches().begin(),
					nend = minst.get_node_patches().end();
				 nit != nend; ++nit)
            {
                size_t patch_idx = (**nit).get_idx();
                fh_input >> f_lhs[patch_idx] >> h_lhs[patch_idx];
                fh_input.ignore(std::numeric_limits<std::streamsize>::max(),
				                '\n' );
            }
            
            fh_input.close();
        }
        
        calculate_energies();
        calculate_forces();
    }
}

void BEEP::init(const ConfigFile& config) {
#ifndef PREHYDROPHOBIC
	resolve = true;	// BEEP needs to solve
#endif
    for (ConfigFile::MeshLibrary::const_iterator
			it=config.mesh_library.cbegin(), end=config.mesh_library.cend();
		 it != end; ++it)
    {
        const std::string& mesh_fname = it->mesh_filename;
        unsigned int mesh_id = it->mesh_id;

        // Enforce that meshes are inserted into the mesh library in the
        // same order as they appear in the config file
        //NB assume this is just a way of enforcing config consistency
        assert(mesh_id == meshes.librarySize());

        // add mesh to library
        meshes.addMesh(mesh_fname, force_planar);
   }

    // Get these from the config.xml (can be overridden by command-line parms)
    Dsolvent = config.solvent_dielectric;
    kappa = (cmdline_kappa == -1) ? config.solvent_kappa : cmdline_kappa;
    
    // loop over the mesh instances in the config layout
    for (ConfigFile::Layout::const_iterator
			it=config.layout.cbegin(), end=config.layout.cend();
		 it != end; ++it)
    {
        // extract the mesh id
        unsigned int mesh_lib_id = it->mesh_id;
        const Vector& offset = it->offset;
        const Quaternion& rotation = it->rotation;

        unsigned int mesh_instance_id = meshes.size();
        assert(mesh_instance_id == it->instance_id);

        unsigned int quad_points_per_triangle=DEFAULT_QUADS;
        if (cmdline_quad_points_per_triangle != -1) {
            quad_points_per_triangle =
				static_cast<unsigned int>(cmdline_quad_points_per_triangle);
        }
        else if (it->quad_points != -1) {
            quad_points_per_triangle =
				static_cast<unsigned int>(it->quad_points);
        }

        unsigned int qual_points_per_triangle=DEFAULT_QUALS;
        if (cmdline_qual_points_per_triangle != -1) {
            qual_points_per_triangle =
				static_cast<unsigned int>(cmdline_qual_points_per_triangle);
        }
        else if (it->qual_points != -1) {
            qual_points_per_triangle =
		 		static_cast<unsigned int>(it->qual_points);
        }
        
        std::cout << "Mesh (lib_id=" << mesh_lib_id << ") @ " << offset
			<< " rot=" << rotation.str() << " quads="
			<< quad_points_per_triangle << " quals="
			<< qual_points_per_triangle << std::endl;

		meshes.add(mesh_lib_id, mesh_instance_id, 
			       offset, rotation, it->dielectric, Dsolvent,
			       quad_points_per_triangle, qual_points_per_triangle);
    }

    assert(meshes.size() > 0);

    // set unique id on each node patch
    unsigned int uid_ctr=0;
    for (MeshInstanceList::iterator mit=meshes.begin(), mend=meshes.end();
		 mit != mend; ++mit)
    {
        MeshInstance& minst = **mit;
        for (PatchList::iterator
				nit = minst.get_node_patches().begin(),
				nend = minst.get_node_patches().end();
			 nit != nend; ++nit)
        {
            (**nit).set_idx(uid_ctr++);
        }
    }
}

Mesh& BEEP::load_library_mesh(const std::string& mtz_filename) {
    // add mesh to library
    boost::shared_ptr<ListedMesh> mesh_ptr
		= meshes.addMesh(mtz_filename, force_planar);
    
    unsigned int num_quad_points = (cmdline_quad_points_per_triangle != -1) ? static_cast<unsigned int>(cmdline_quad_points_per_triangle) : DEFAULT_QUADS;
    unsigned int num_qual_points = (cmdline_qual_points_per_triangle != -1) ? static_cast<unsigned int>(cmdline_qual_points_per_triangle) : DEFAULT_QUALS;
    
    mesh_ptr->set_quad_points_per_triangle(num_quad_points);
    mesh_ptr->set_qual_points_per_triangle(num_qual_points);
    
    return *mesh_ptr;
}

void BEEP::clear_mesh_instances(unsigned int start, int end) {
#ifndef PREHYDROPHOBIC
	resolve = true;	// BEEP needs to solve
#endif
	// This allows python-like end-point specification with +1 increment
	if (start >= meshes.size()) return;  // Silently ignore (as if already done)
	if (end < 0) end = meshes.size() + 1 + end;
	if (end < start) return;  // no-op
	// Silently translate beyond-the-end to just end
    MeshInstanceList::iterator endit = meshes.end();
	if (end < meshes.size()) endit = meshes.begin()+end;
    meshes.erase(meshes.begin()+start, endit);
}

MeshInstance& BEEP::insert_mesh_instance(
	unsigned int mesh_lib_id,
	const Vector& offset, const Quaternion& rotation,
	const double protein_dielectric)
{
    // check that the mesh_lib_id is valid
    if (mesh_lib_id >= meshes.librarySize()) {
        std::cerr << "Bad value for mesh_lib_id: " << mesh_lib_id
			<< " (" << meshes.librarySize() << " meshes defined)" << std::endl;
        throw std::exception();
    }

    // load mesh instance
    unsigned int num_quad_points = (cmdline_quad_points_per_triangle != -1)
		? static_cast<unsigned int>(cmdline_quad_points_per_triangle)
		: DEFAULT_QUADS;
    unsigned int num_qual_points = (cmdline_qual_points_per_triangle != -1)
		? static_cast<unsigned int>(cmdline_qual_points_per_triangle)
		: DEFAULT_QUALS;

	unsigned int id = meshes.size();
	MeshInstance& m = *meshes.add(mesh_lib_id, id, offset, rotation, 
								  protein_dielectric, Dsolvent,
								  num_quad_points, num_qual_points);
#ifndef PREHYDROPHOBIC
	resolve = true;				// BEEP needs to solve
	m.update_energy(meshes);	// Revise non-electrostatic energies

	// Store the reversion data
	previd = meshes.size();		// No reversion is possible from an insert
#endif
	return m;
}

MeshInstance& BEEP::move_mesh_instance(
	unsigned int id,
	const Vector& offset, const Quaternion& rotation,
	const double protein_dielectric)
{
	static const Quaternion no_rotation{Quaternion(1,0,0,0)};

#ifdef MI_MOVE
    // check that the mesh id is valid
    if (id >= meshes.size()) {
        std::cerr << "Bad value for mesh instance_id: " << id
			<< " (" << meshes.size() << " mesh instances defined)" << std::endl;
        throw std::exception();
    }
	// Ignore null moves
	MeshInstance& m = *(meshes[id]);
	Vector loc = m.get_xyz_offset();	// Copy
	if (offset == loc && rotation == no_rotation) return m;

	// Make the move
	m.move(offset, rotation);
#else  // MI_MOVE
	MeshInstance& m = *meshes.move(id, offset, rotation,
									protein_dielectric, Dsolvent);
#endif // MI_MOVE
#ifndef PREHYDROPHOBIC
	// Energy updates are required
	resolve = true;		// BEEP needs to solve
	if (id == previd && offset == prevloc && rotation == unrot)
		m.revert_energy(meshes); // This is a move reversal, revert quickly
	else
		m.update_energy(meshes); // Revise non-electrostatic energies now

	// Store the reversion data
	previd = id;
	prevloc = loc;				// offset is absolute location
	unrot = rotation.inverse();	// but rotation is relative
#endif
	return m;
}

int BEEP::get_instance_id(const Vector& pt, int skip_id) const {
    for (MeshInstanceList::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        const MeshInstance& minst = **mit;
		if (minst.get_id() == skip_id) continue;
		if (minst.pt_is_internal(pt)) return minst.get_id();
	}
	return -1;  // Invalid id = no id TODO (should be a constexpr)
}

#ifdef FUTURE
void BEEP::create_kinemage(
	const std::string& filename,
	int num_colours,
	const std::string& preamble) const
{
	create_kinemage(filename, 0.0, 0.0, num_colours, preamble);
}
#endif

// Future preference is discard this version, merge into above without scales
// pybeep.py will then need updating
void BEEP::create_kinemage(
	const std::string& filename,
	double fscale, double hscale, int num_colours,
	const std::string& preamble) const
{
    std::ostringstream buf_f, buf_sf, buf_h, buf_sh;
    buf_f << "@trianglelist {fvals} off\n";
    buf_h << "@trianglelist {hvals} \n";
    buf_sf << "\n@trianglelist {silent-fvals} off\n";
    buf_sh << "\n@trianglelist {silent-hvals} \n";
#ifndef PREHYDROPHOBIC
    std::ostringstream buf_hy, buf_s, buf_e, buf_he, buf_lj;
    buf_hy << "\n@trianglelist {hydrophobicities} off\n";
    buf_s << "\n@trianglelist {Lennard-Jones-sigma} off\n";
    buf_e << "\n@trianglelist {Lennard-Jones-epsilon} off\n";
    buf_he << "\n@trianglelist {Hydrophobic Effect} off\n";
    buf_lj << "\n@trianglelist {Lennard-Jones} off\n";
	double scales[7] = {0.0, 0.0, 0.0, 0.0, 0.0, fscale, hscale};
#endif // PREHYDROPHOBIC
    for (const auto& spmi: meshes) {
        if (spmi->isSilent() == false)
#ifdef PREHYDROPHOBIC
            spmi->kinemage_fh_vals(fscale, hscale, num_colours, buf_f, buf_h);
        else
            spmi->kinemage_fh_vals(fscale, hscale, num_colours, buf_sf, buf_sh);
#else
            spmi->kinemage_vals(num_colours, buf_hy, buf_s, buf_e,
								buf_he, buf_lj, buf_f, buf_h, scales);
		else
            spmi->kinemage_vals(num_colours, buf_hy, buf_s, buf_e,
								buf_he, buf_lj, buf_sf, buf_sh, scales);
#endif // PREHYDROPHOBIC
    }

    std::ostringstream charge_buf;
    charge_buf << "@spherelist {charges}\n";
    for (const auto& spmi: meshes) {
//TODO allCharges now?
        for (const auto& spch: spmi->get_charges()) {
            const Charge& ch = *spch;
            charge_buf << "r=" << ch.get_radius() << " {} "
			           << ch.x << " " << ch.y << " " << ch.z << "\n";
        }
    }    
    
    std::ofstream kin_out(filename.c_str());
    kin_out << "@kinemage\n";
    kin_out << preamble << std::endl;
    kin_out << buf_f.str() << std::endl;
    kin_out << buf_h.str() << std::endl;
    kin_out << buf_sf.str() << std::endl;
    kin_out << buf_sh.str() << std::endl;
    kin_out << charge_buf.str() << std::endl;
#ifndef PREHYDROPHOBIC
    kin_out << buf_hy.str() << std::endl;
    kin_out << buf_s.str() << std::endl;
    kin_out << buf_e.str() << std::endl;
    kin_out << buf_he.str() << std::endl;
    kin_out << buf_lj.str() << std::endl;
#endif // PREHYDROPHOBIC
    kin_out.close();
}

// void BEEP::evaluate_ecm()
// {
//     // for each mesh instance which has effective charges defined (in the .mtz file) evaluate
//     // the energy according to the potential at the ecm position multiplied by the effective 
//     // charge
//     
// }

std::string BEEP::benchmark() {
    
    //std::cout << "sizeof a node patch and QP is: " << sizeof(NodePatch) << " " << sizeof(QuadPoint) << std::endl;

	unsigned int num_patches = reset_fh_vals();

    boost::scoped_ptr<double> x(new double[2*num_patches]);
    double* xx = x.get();
 
    // init FMM octree
    calculate_bounding_cube();
    long init_fmm_time = init_fmm_octree();
    std::cout << "Init FMM tree took: "<< init_fmm_time / 1000. << " ms"
	          << std::endl;
    std::cout << "total num patches = " << num_patches << std::endl;
    
    long time_start = myclock();
    long precalc = precalc_neighbour_interactions();
    long bem_time = do_bem_fmm_iteration(num_patches, xx, xx);
        
//     // Create the return value RunInfo
//     RunInfo stats;
//     
//     stats.bem_fmm_per_iter = fmm->get_timing_info() / 1000.;
//     stats.num_iterations = 1;
//     stats.time_total = (myclock() - time_start) / 1000.;
//     stats.num_explicits = count_explicit_integrations();
//     stats.num_patches = num_patches;
//     stats.neighbourhood_size = average_bem_explicit_neighbs;
//     
    // free the memory associated with the FMM
    fmm.reset();
    
    std::stringstream buf;
    buf << num_patches << "," << num_patches*num_patches << ","
	    << bem_time / 1000.;
    
    return buf.str();
}

RunInfo BEEP::solve(double gmres_stop_criteria, int max_iterations) {
#ifndef PREHYDROPHOBIC
	resolve = false;	// BEEP is now calculating electrostatics...
#endif
    vanilla_fmm_timer.zero();
    bem_fmm_timer.zero();
    
    long time_start_solve = myclock();
    
    // num patches
	unsigned int num_patches = reset_fh_vals();

    // Create the return value RunInfo
    RunInfo stats;

    // calculate the bounding cube (centre and edge length of universe)
    // NB: this can be set manually using the set_bounding_cube function
    // in which case the calculate bounding cube function will just return
    calculate_bounding_cube();

    // set the RHS vector
    long time_rhs = set_rhs();   // updates rhs_octree attribute
    stats.rhs_fmm = rhs_octree->get_timing_info();
    rhs_octree.reset();
        
    // init FMM octree  
    long time_init_fmm = init_fmm_octree();  // updates fmm attribute

    // pre-calc the local interactions
    size_t num_explicit_integrations = count_explicit_integrations();
    long time_precalcs = precalc_neighbour_interactions();

    // Run iterative GMRES solver to get BEM solution for electrostatics
    long time_gmres = myclock();
    unsigned int num_iterations = gmres(gmres_stop_criteria, max_iterations);
    time_gmres = myclock() - time_gmres;

    // set the f/h values within the node patches to their new converged values
    for (std::vector< boost::shared_ptr<MeshInstance> >::iterator
			mit=meshes.begin(), mend=meshes.end();
		 mit != mend; ++mit)
    {
        MeshInstance& minst = **mit;
        for (PatchList::iterator
				nit = minst.get_node_patches().begin(),
				nend = minst.get_node_patches().end();
			 nit != nend; ++nit)
        {
            BasicNodePatch& np = **nit;
            unsigned int patch_idx = np.get_idx();
            np.f = f_lhs[patch_idx];
            np.h = h_lhs[patch_idx];
        }
    }
    
    // calculate energies/forces on the molecules
    //double total_energy = calculate_energies();
    //Vector total_force = calculate_forces();

    //std::cout << "Total energy: " << total_energy << std::endl;

//     std::cout << "f/h vals:\n";
//     for (unsigned int ii=0; ii < 10; ++ii)
//     {
//         const BasicNodePatch& np = *(get_patch(ii));
//         std::cout << np.get_normal() << " " << np.get_alt_normal() << " " << f_lhs[ii] << " " << h_lhs[ii] << std::endl;
//     }

    long total_time = myclock() - time_start_solve;
    
/*
    vanilla_fmm_timer /= 1000.; // convert to ms
    bem_fmm_timer /= 1000.; // convert to ms
    
    std::cout << "Times:\n";
    std::cout << "(rhs, init_fmm, precalcs, gmres, total)\n";
    std::cout << time_rhs << ", " << time_init_fmm << ", " << time_precalcs << ", " << total_gmres << ", " << total_time << std::endl;
    std::cout << "FMM data (vanilla, bem/fmm):\n";
    std::cout << vanilla_fmm_timer.describe();
    std::cout << vanilla_fmm_timer << bem_fmm_timer << std::endl;*/


    stats.bem_fmm_per_iter = fmm->get_timing_info();
    stats.num_iterations = num_iterations;
    stats.time_total = (myclock() - time_start_solve);
    stats.time_rhs = time_rhs;
    stats.time_gmres = time_gmres;
    stats.time_precalcs = time_precalcs;
    stats.num_explicits = num_explicit_integrations;
    stats.num_patches = num_patches;
    stats.neighbourhood_size = (cmdline_bem_nbsize == -1)
		? DEFAULT_BEM_NEIGHBOURHOOD_SIZE : cmdline_bem_nbsize;
    
    // free the memory associated with the FMM
    fmm.reset();
    
    return stats;
}

size_t BEEP::reset_fh_vals() {
    // resets the f/h node patch values of the mesh definitions
	unsigned int num_patches = get_total_patches();

    std::cout << "total num patches = " << num_patches << std::endl;
    
    // allocate memory - initialised in set_rhs
    f_rhs = boost::shared_array<double>(new double[num_patches]);
    h_rhs = boost::shared_array<double>(new double[num_patches]);
    f_lhs = boost::shared_array<double>(new double[num_patches]);
    h_lhs = boost::shared_array<double>(new double[num_patches]);

    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        f_rhs[ctr] = 0.0;
        h_rhs[ctr] = 0.0;
        f_lhs[ctr] = 0.0;
        h_lhs[ctr] = 0.0;
    }
	return num_patches;
}

void BEEP::reset_library_fh_vals() {
    meshes.reset_library_fh_vals(f_lhs, h_lhs);
}

void BEEP::write_fh(const std::string& output_filename) {
    if (output_filename == "") { return; }

    std::ofstream fh_output;
    
    for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        MeshInstance& minst = **mit;
        
        // silent meshes don't write their fh values
        if (minst.isSilent()) continue;
        
        std::stringstream name_buf;
        name_buf << output_filename << "." << minst.get_id();
        std::cout << "Writing to: " << name_buf.str() << std::endl;
        fh_output.open(name_buf.str().c_str(), std::ios_base::out);
        
        fh_output << std::setprecision(12);
        for (PatchList::iterator
				nit = minst.get_node_patches().begin(),
				nend = minst.get_node_patches().end();
			 nit != nend; ++nit)
        {
            size_t patch_idx = (**nit).get_idx();
            fh_output << f_lhs[patch_idx] << " " << h_lhs[patch_idx] << "\n";
        }
        
        fh_output.close();
    }
}

double BEEP::calculate_energies() {
    double total_E = 0.0;

	// TODO would be nice to avoid any calculations if not required:
	// i.e. add double base_energy to reset with resolve bool set true?
	// Then resolve && base_energy != unset => return base_energy

    unsigned int offset=0;
    for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        const MeshInstance& m = **mit;
        if (m.isSilent()) continue;
#ifdef PREHYDROPHOBIC
        total_E += m.calculate_energy(kappa, &f_lhs[offset], &h_lhs[offset]);
#else
		// If resolve required, then do not calculate electrostatics...
        total_E += m.calculate_energy(!resolve, kappa,
									  &f_lhs[offset], &h_lhs[offset]);

#endif // PREHYDROPHOBIC
        offset += m.get_num_node_patches();
    }
    return total_E;
}

//TODO hydrophobic etc forces?
void BEEP::calculate_forces() {
    unsigned int offset=0;
    std::cout << std::setprecision(6);
    for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        const MeshInstance& m = **mit;
        
        Vector h_squared = m.get_h_squared();
        std::cout << "Total h-squared: " << h_squared << std::endl;
        
        KahanVector qE;
        KahanVector MST_external;
        KahanVector MST_internal;
        KahanVector ionic;
        KahanVector dbf;
        m.calculate_forces(kappa, &f_lhs[offset], &h_lhs[offset], qE,
		                   MST_external, MST_internal, dbf, ionic);
        
        std::cout << "Forces for mesh_id=" << m.get_id() << "\n";
        std::cout << " qE (patch-charge): " << (*qE).length() << " " << *qE 
		          << "\n";
        std::cout << " qE (MST internal): " << (*MST_internal).length()
		          << " " << *MST_internal <<"\n";
        std::cout << " dbf (delta_MST): " << (*dbf).length() << " " << *dbf 
				  << "\n";
        std::cout << " ionic: " << (*ionic).length() << " " << *ionic <<"\n";
        
        KahanVector total_patch = qE + dbf + ionic;
        KahanVector total_mst_int = MST_internal + dbf + ionic;
        KahanVector total_mst_ext = MST_external + ionic;
        std::cout << "Total forces (take your pick...)\n";
        std::cout << "total (patch-charge + dbf + ionic): "
		          << (*total_patch).length() << " " << *total_patch << "\n";
        std::cout << "total (MST internal + dbf + ionic): "
		          << (*total_mst_int).length() << " " << *total_mst_int << "\n";
        std::cout << "total (MST external + ionic): "
		          << (*total_mst_ext).length() << " " << *total_mst_ext <<"\n";
        //std::cout << "magnitude: " << *total_force.length() << "\n";
        
        //std::cout << "Force on mesh_id=" << m.instance_id << " (lib_id=" << m.lib_id << "): " << std::setprecision(10) << total_force << std::endl;
        offset += m.get_num_node_patches();
    }
    std::cout << std::flush;
}

void BEEP::write_matrix() const {
    unsigned int num_patches = get_total_patches();
    
    // precalc quad points
    std::vector< boost::shared_ptr<QuadList> > cache;
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        const BasicNodePatch& patch = get_patch(ii);
        patch.obtain_shared_quad_ptrs(cache);
    }
    
    std::ofstream fA,fB,fC,fD;
    fA.open("A.matrix", std::ios_base::out);
    fB.open("B.matrix", std::ios_base::out);
    fC.open("C.matrix", std::ios_base::out);
    fD.open("D.matrix", std::ios_base::out);
    
    std::cout << "Writing full BEM matrix... ("
	          << num_patches*num_patches*4 << " elements)" << std::endl;
    
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        const BasicNodePatch& src_patch = get_patch(ii);
        for (unsigned int jj=0; jj < num_patches; ++jj) {
            float A=0,B=0,C=0,D=0;
            LocalIntegrations integral;
            const BasicNodePatch& targ_patch = get_patch(jj);
            if (ii==jj) {
                singular_BEM_kernels(kappa, src_patch, A, B, C, D);
                B += fGeometricCorrection(src_patch.gc, src_patch.dielectric_ratio);
                C -= hGeometricCorrection(src_patch.gc, src_patch.dielectric_ratio);
            }
            else {
                integral.init(kappa, src_patch, targ_patch);
                A = integral.Apt;
                B = integral.Bpt;
                C = integral.Cpt;
                D = integral.Dpt;
            }
            fA << -A << " ";
            fB << B << " ";
            fC << -C << " ";
            fD << D << " ";
        }
        fA << "\n";
        fB << "\n";
        fC << "\n";
        fD << "\n";
    }

    fA.close();
    fB.close();
    fC.close();
    fD.close();
    
    std::ifstream fiA,fiB,fiC,fiD;
    fiA.open("A.matrix", std::ios_base::in);
    fiB.open("B.matrix", std::ios_base::in);
    fiC.open("C.matrix", std::ios_base::in);
    fiD.open("D.matrix", std::ios_base::in);
    
    std::ofstream final;
    final.open("matrix.dat", std::ios_base::out);
    
    const size_t BUFSIZE = 64*1024;
    boost::scoped_array<char> bufA(new char[BUFSIZE]);
    boost::scoped_array<char> bufB(new char[BUFSIZE]);
    
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        //assert(fiB.good() && fiA.good());
        fiA.getline(bufA.get(), BUFSIZE, '\n');
        fiB.getline(bufB.get(), BUFSIZE, '\n');
        final << std::string(bufB.get()) << " " << std::string(bufA.get())
		      << "\n";
    }
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        //assert(fiD.good() && fiC.good());
        fiC.getline(bufA.get(), BUFSIZE, '\n');
        fiD.getline(bufB.get(), BUFSIZE, '\n');
        final << std::string(bufB.get()) << " " << std::string(bufA.get())
		      << "\n";
    }
    fiA.close();
    fiB.close();
    fiC.close();
    fiD.close();
    final.close();
}

size_t BEEP::count_explicit_integrations() const {
    // dry run of precalc_neighbour_interactions to determine how much memory
    // is required (avoiding using all available memory then failing...)
    unsigned int total_local_ints=0;
    const FMM_BEM_Octree& fmm_octree = *fmm;
    return fmm_octree.calc_neighbourhood_interacts();

//     for (unsigned short level=fmm->get_top_level(); level <= fmm->get_bottom_level(); ++level)
//     {
//         FMM_BEM_Octree::NodeList& nodes = fmm->get_node_list(level);
//         for (FMM_BEM_Octree::NodeList::const_iterator it=nodes.begin(), end=nodes.end(); it != end; ++it)
//         {
//             // slightly cumbersome syntax- it's because NodeList is a
//             // sdt::map not a vector/list
//             const FMM_BEM_Octree::NodeT& node = *(it->second);
// 
//             // skip empty nodes and those which are under leaf nodes; only count
//             // nodes where we will actually evaluate BEM/FMM
//             if (node.isLeaf() == false || node.empty()) { continue; }
// 
//             const std::vector< BasicNodePatch* >& leaf_contents = node.get_contents();
//             boost::shared_ptr< std::vector<BasicNodePatch*> > neighbourlist = fmm->get_neighbourhood_contents(node);
//             total_local_ints += neighbourlist->size() * leaf_contents.size();
//         }
//     }
// 
//     return total_local_ints;
}

double BEEP::calc_local_neighbourhood_memory_requirement() const {
    size_t total_local_ints = count_explicit_integrations();
    const double megs = static_cast<double>(
					sizeof(LocalIntegrations)*total_local_ints) /(1024*1024);
    std::cout << total_local_ints << " local integrations --> " << megs
	          << " MB" << std::endl;
    return megs;
}

long BEEP::precalc_neighbour_interactions() {
    if (skipping_precalcs) return 0;

    long start_precalcs = myclock();
    
    //std::cout << "Precalculating neighbour interactions..." << std::endl;
    unsigned int num_patches = get_total_patches();
    
    // clear any previously stored local integration results
    local_integrations.clear();

    // create local integrations for the self-geometric interactions
    LintArray local_ints(new LocalIntegrations[num_patches]);
    local_integrations.push_back(LintArray_Size(local_ints, num_patches));
    
    long singular_clock = myclock();
    
#ifdef OPENCL

    // cauchy principle value part
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        const BasicNodePatch& np = get_patch(ii);
        double dielectric_ratio = np.get_dielectric_ratio();
        local_ints[ii].set(ii,
                           ii,
                           0,
                           fGeometricCorrection(np.gc, dielectric_ratio),
                           -hGeometricCorrection(np.gc, dielectric_ratio),
                           0);
    }
    
    // create local integrations for the self-geometric interactions
    LintArray gpu_local_ints(new LocalIntegrations[num_patches]);
    local_integrations.push_back(LintArray_Size(gpu_local_ints, num_patches));

    unsigned int chunksize = 128;
    unsigned int ii=0;
    while (ii < num_patches) {
  
        PatchPtrList patch_ptrs(new PPList);
        
        unsigned chunk_ctr=0;
        for ( ; ii+chunk_ctr < num_patches && chunk_ctr < chunksize;
		     ++chunk_ctr)
        {
            const BasicNodePatch& np = get_patch(ii+chunk_ctr);
            patch_ptrs->push_back(&np);
        }
        
        // singular part
        SingularBEM* ocl_singular_bem
				= new SingularBEM(patch_ptrs, kappa, &(gpu_local_ints[ii]));
        global_ocl_handler.add_work_to_queue(ocl_singular_bem);
        
        ii += chunk_ctr;
    }
    
    // block until all GPU work is complete
    //global_ocl_handler.wait_until_idle();
    
#else  // OPENCL
    for (unsigned int ii=0; ii < num_patches; ++ii) {
        const BasicNodePatch& np = get_patch(ii);
        double dielectric_ratio = np.get_dielectric_ratio();
        
        float A=0,B=0,C=0,D=0;
        //std::cout << ii << std::endl;
        singular_BEM_kernels(kappa, np, A, B, C, D);
        
        local_ints[ii].set(ii,
                           ii,
                           A,
                           B + fGeometricCorrection(np.gc, dielectric_ratio),
                           C - hGeometricCorrection(np.gc, dielectric_ratio),
                           D);
    }
#endif  // OPENCL
    std::cout << num_patches << " singular integrals took: "
	          << (myclock() - singular_clock) / 1000. << " ms" << std::endl;

    // figure out how much memory will be required
    calc_local_neighbourhood_memory_requirement();

#ifndef CACHE_GPU_RESULTS
    std::cout << "(Not caching the explicit-integration results anyway "
				 "though...)" << std::endl;
#endif
    
#ifdef CACHE_GPU_RESULTS

    long precalcs_clock = myclock();
    long running_ctr=0;
    
    // iterate over all leaf nodes of the fmm
    // for each leaf node, dump a list of the contents, and the neighbourhood
    // for each node patch of the leaf, integrate over all neighbours

    for (unsigned short level=fmm->get_top_level();
	     level <= fmm->get_bottom_level(); ++level)
    {
        FMM_BEM_Octree::NodeList& nodes = fmm->get_node_list(level);
        
        // quadrature_point cache
        std::vector<boost::shared_ptr<QuadList> > qp_cache;

        for (FMM_BEM_Octree::NodeList::const_iterator
				it=nodes.cbegin(), end=nodes.cend();
			 it != end; ++it)
        {
            // slightly cumbersome syntax- it's because NodeList is a
            // std::map not a vector/list
            const FMM_BEM_Octree::NodeT& node = *(it->second);

            // skip empty nodes and those which are under leaf nodes
            if (node.isLeaf() == false || node.empty()) continue;

            // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
            boost::shared_ptr<std::vector<BasicNodePatch*> > neighbourlist
				= fmm->get_neighbourhood_contents(node);
            size_t num_neighbours = neighbourlist->size();
            size_t ctr=0;
            size_t chunksize = static_cast<size_t>(ceil(static_cast<double>
				(BEM_EXPLICIT_CHUNKSIZE) / node.size()) );
            chunksize = (num_neighbours > chunksize)
				? chunksize : num_neighbours;

            while (ctr < num_neighbours) {
                boost::shared_ptr< std::vector<const BasicNodePatch*> >
					list_chunk(new std::vector<const BasicNodePatch*>);
                list_chunk->reserve(chunksize);

                for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
                {
                    ((*neighbourlist)[ctr])->obtain_shared_quad_ptrs(qp_cache);
                    list_chunk->push_back((*neighbourlist)[ctr++]);
                }

                size_t size = list_chunk->size() * node.size();
                LintArray results(new LocalIntegrations[size]);
                assert(results.get()); // check memory allocation didn't fail
                LintArray_Size arr_siz(results, size);
                local_integrations.push_back( arr_siz );

#ifdef OPENCL
                //std::cout << "Size of neighbourlist: " << neighbourlist->size() << std::endl;
                PatchPtrList cntnts(new PPList);
                cntnts->insert(cntnts->begin(), node.get_contents().begin(),
				                                node.get_contents().end());
                BEM_Resources* res_ptr
					= new BEM_Resources(cntnts, list_chunk,
					                    kappa, results.get());
                global_ocl_handler.add_work_to_queue(res_ptr);
#else  // OPENCL
                for (int src_ctr=0; src_ctr < node.get_contents().size();
					 ++src_ctr)
                {
                    const BasicNodePatch& src_patch
						= *(node.get_contents()[src_ctr]);
                    #pragma omp parallel for
                    for (int targ_ctr=0; targ_ctr < list_chunk->size();
					     ++targ_ctr)
                    {
                        const BasicNodePatch& targ_patch
							= *((*list_chunk)[targ_ctr]); // that's a bit ugly
                        
                        size_t idx = (src_ctr * list_chunk->size()) + targ_ctr;
                        results[idx].init(kappa, src_patch, targ_patch);
                    }
                }
                running_ctr += size;
                std::cout << running_ctr << "\r" << std::flush;
#endif  // OPENCL
            }
        }
        
#ifdef OPENCL
        // ensure that QuadPoints stay in scope for OpenCL
        OpenCL_WorkBlob* ocl_cache_ptr = new QuadPointCache(qp_cache);
        global_ocl_handler.add_work_to_queue(ocl_cache_ptr);
#endif
    }
    
    //std::cout << "Done " << local_integrations.size() << " groups of pair-wise local neighbour interactions" << std::endl;
    std::cout << "Precalc integrals took: "
	          << (myclock() - precalcs_clock) / 1000. << " ms" << std::endl;
#endif  // CACHE_GPU_RESULTS

    return (myclock() - start_precalcs);
}

void BEEP::evaluate_local_neighbours(double f_results[], double h_results[]) {

#ifndef CACHE_GPU_RESULTS
    for (unsigned short level=fmm->get_top_level();
		 level <= fmm->get_bottom_level(); ++level)
    {
        // quadrature_point cache
        std::vector<boost::shared_ptr<QuadList> > qp_cache;

        FMM_BEM_Octree::NodeList& nodes = fmm->get_node_list(level);
        for (FMM_BEM_Octree::NodeList::const_iterator
				it=nodes.cbegin(), end=nodes.cend();
			 it != end; ++it)
        {
            // slightly cumbersome syntax- it's because NodeList is a
            // std::map not a vector/list
            const FMM_BEM_Octree::NodeT& node = *(it->second);

            // skip empty nodes and those which are under leaf nodes
            if (node.isLeaf() == false || node.empty()) continue;

            // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
            boost::shared_ptr<std::vector<BasicNodePatch*> > neighbourlist
				= fmm->get_neighbourhood_contents(node);
            size_t num_neighbours = neighbourlist->size();
            size_t ctr=0;
            size_t chunksize=static_cast<size_t>(ceil(static_cast<double>
				(BEM_EXPLICIT_CHUNKSIZE) / node.size()) );
            chunksize = (num_neighbours > chunksize)
				? chunksize : num_neighbours;

            while (ctr < num_neighbours) {
                boost::shared_ptr< std::vector<const BasicNodePatch*> >
					list_chunk(new std::vector<const BasicNodePatch*>);
                list_chunk->reserve(chunksize);

                for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
                {
                    ((*neighbourlist)[ctr])->obtain_shared_quad_ptrs(qp_cache);
                    list_chunk->push_back((*neighbourlist)[ctr++]);
                }
#ifdef OPENCL
                PatchPtrList cntnts(new PPList);
                cntnts->insert(cntnts->begin(), node.get_contents().begin(),
				                                node.get_contents().end());
                
                //std::cout << "Size of neighbourlist: " << neighbourlist->size() << std::endl;
                BEM_OnDemand_Resources* res_ptr
					= new BEM_OnDemand_Resources(cntnts, list_chunk, kappa,
							 f_lhs.get(), h_lhs.get(), f_results, h_results);
                global_ocl_handler.add_work_to_queue(res_ptr);
#else  // OPENCL
                LocalIntegrations bem_kernels;

                #pragma omp parallel for private(bem_kernels)
                for (int src_ctr=0; src_ctr < node.get_contents().size();
					 ++src_ctr)
                {
                    const BasicNodePatch& src_patch
						= *(node.get_contents()[src_ctr]);
                    for (int targ_ctr=0; targ_ctr < list_chunk->size();
						 ++targ_ctr)
                    {
                        const BasicNodePatch& targ_patch
							= *((*list_chunk)[targ_ctr]); // that's a bit ugly
                        bem_kernels.init(kappa, src_patch, targ_patch);
                        bem_kernels.evaluate_local_contributions
							(f_lhs.get(), h_lhs.get(), f_results, h_results);
                    }
                }
#endif  // OPENCL
            }

        }
        
#ifdef OPENCL
        // ensure that QuadPoints stay in scope for OpenCL
        OpenCL_WorkBlob* ocl_cache_ptr = new QuadPointCache(qp_cache);
        global_ocl_handler.add_work_to_queue(ocl_cache_ptr);
#endif  // OPENCL
    }
#endif  // CACHE_GPU_RESULTS

#ifdef OPENCL
    // now wait for OpenCL to complete -- in case of precalcs still remaining!
    global_ocl_handler.wait_until_idle();
    sanity_check();
#endif // OPENCL

    // loop over the precalculated local integrations and multiply f/hvals
    // by the matrix elements
    for (std::vector<LintArray_Size>::const_iterator
			it=local_integrations.cbegin(), end=local_integrations.cend();
		 it != end; ++it)
    {
        const LintArray& array = it->first;
        const size_t& sz = it->second;
        //std::cout << "Evaluating: " << sz << " local ints, from mem: " << array.get() << std::endl;
        for (size_t ii=0; ii < sz; ++ii)
        {
            //std::cout << ii << " of " << sz << ": " << array[ii] << "\n";
            array[ii].evaluate_local_contributions
				(f_lhs.get(), h_lhs.get(), f_results, h_results);
        }
    }
    return;
}

void BEEP::sanity_check() {
    static bool done_sanity = false;
    if (done_sanity) { return ; }

    unsigned int num_patches = get_total_patches();
    double max_val = 1e6;
    bool force_recalc=false;

    do {
        // check local integrations for bad values
        // (can happen because of GPU memory glitches)
        for (std::vector<LintArray_Size>::const_iterator
				it=local_integrations.cbegin(), end=local_integrations.cend();
			 it != end; ++it)
        {
            const LintArray& array = it->first;
            const size_t& sz = it->second;
            //std::cout << "Evaluating: " << sz << " local ints, from mem: " << array.get() << std::endl;
            for (size_t ii=0; ii < sz; ++ii) {

                if (array[ii].src_global_idx >= num_patches
					|| array[ii].targ_global_idx >= num_patches) {
                    std::cout << "Bad index detected, repeating precalcs..."					          << std::endl;
                    force_recalc = true;
                    break;
                }
                if (array[ii].insane(max_val)) {
                    std::cout << "Detected a crazy value in near-field "
								 "integrations: " << array[ii] << std::endl;
                    const BasicNodePatch& src
						= get_patch(array[ii].src_global_idx);
                    const BasicNodePatch& targ
						= get_patch(array[ii].targ_global_idx);
                    if (&src == &targ) {
						// ill-conditioned patch- kill it
                        array[ii].set(0,0,0,0,0,0);
                    }
                    else {
                        array[ii].init(kappa, src, targ);
                    }
                }
            }

            if (force_recalc) break;
        }

        if (force_recalc) {
            precalc_neighbour_interactions();
    
#ifdef OPENCL
            global_ocl_handler.wait_until_idle();
#endif  // OPENCL
        }
    } while (force_recalc);
 
    done_sanity = true;
}

long BEEP::set_rhs() {
    using fmm::FMM_Octree_6FIG_ACCURACY;
    using fmm::EvalPt;
    
    long start_rhs_clock = myclock();

    // first get the bounding cube for the mesh
    Vector centre = bounding_cube_centre;
    double edge_length = bounding_cube_edge_length;
std::cout << "BEEP::set_rhs MAX FMM SIZE " << MAX_FMM_SIZE << std::endl;
    rhs_octree.reset(new FMM_Octree_6FIG_ACCURACY
								(MAX_FMM_SIZE, centre, edge_length) );
    
    std::vector<EvalPt*> eval_pts_rhs;
    size_t num_qps=0;
    // loop over each mesh instance and insert the charges into a big fat octree
    for (MeshInstanceList::const_iterator it=meshes.cbegin(), end=meshes.cend();
		 it != end; ++it)
    {
        const MeshInstance& mesh_instance = **it;
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
				np_it = mesh_instance.get_node_patches().cbegin(),
				np_end = mesh_instance.get_node_patches().cend();
			 np_it != np_end; ++np_it)
        {
            const BasicNodePatch& np = **np_it;
            boost::shared_ptr<QuadList> qual_pts = np.get_qualocation_points();
			// NB why not num_qps += qual_pts->cend() - qual_pts->cbegin() ??
            for (QuadList::const_iterator
					qp_it = qual_pts->cbegin(), qp_end=qual_pts->cend();
				 qp_it != qp_end; ++qp_it)
            {
                ++num_qps;
            }
        }
    }
    eval_pts_rhs.reserve(num_qps);

    unsigned int total_charges=0;
    // loop over each mesh instance and insert the charges into a big fat octree
    for (MeshInstanceList::const_iterator it=meshes.cbegin(), end=meshes.cend();
		 it != end; ++it)
    {
        const MeshInstance& mesh_instance = **it;
        for (std::vector< boost::shared_ptr<Charge> >::const_iterator
				ch_it = mesh_instance.get_charges().cbegin(),
				ch_end = mesh_instance.get_charges().cend();
			 ch_it != ch_end; ++ch_it)
        {
            const Charge& ch = **ch_it;
            rhs_octree->add(ch);
            ++total_charges;
        }

        for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
				np_it = mesh_instance.get_node_patches().cbegin(),
				np_end = mesh_instance.get_node_patches().cend();
			 np_it != np_end; ++np_it)
        {
            const BasicNodePatch& np = **np_it;
	    
#ifdef CENTROID_COLLOCATION
            eval_pts_rhs.push_back( new EvalPt( (Vector&) np.get_node() ) );
#else
            boost::shared_ptr<QuadList> qual_pts = np.get_qualocation_points();
            for (QuadList::const_iterator
					qp_it = qual_pts->cbegin(), qp_end=qual_pts->cend();
				 qp_it != qp_end; ++qp_it)
            {
                eval_pts_rhs.push_back( new EvalPt( Vector(qp_it->pt()) ) );
            }
#endif  // CENTROID_COLLOCATION
        }

    }

    std::cout << "total charges: " << total_charges << " eval points: "
	          << eval_pts_rhs.size() << std::endl;

    rhs_octree->solve(0.0);
#ifdef OPENCL
    rhs_octree->evaluate_many(eval_pts_rhs, global_ocl_handler);
#else
    rhs_octree->evaluate_many(eval_pts_rhs);
#endif // OPENCL

    // loop over node patches
    size_t ctr=0, ep_ctr=0;
    for (MeshInstanceList::const_iterator it=meshes.cbegin(), end=meshes.cend();
		 it != end; ++it)
    {
        const MeshInstance& mesh_instance = **it;
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
				np_it = mesh_instance.get_node_patches().cbegin(),
				np_end = mesh_instance.get_node_patches().cend();
			 np_it != np_end; ++np_it)
        {
            const BasicNodePatch& np = **np_it;
            f_rhs[ctr] = 0;
            h_rhs[ctr] = 0;

#ifdef CENTROID_COLLOCATION
            EvalPt* ep_ptr = eval_pts_rhs[ep_ctr++];
            f_rhs[ctr] += ep_ptr->get_potential() * ONE_OVER_4PI / Dsolvent;
            h_rhs[ctr] += ep_ptr->get_field().dot(np.get_normal())
						  * ONE_OVER_4PI / Dsolvent;
            delete ep_ptr;

#else  // CENTROID_COLLOCATION
            boost::shared_ptr<QuadList> qual_pts = np.get_qualocation_points();
            for (QuadList::const_iterator
					qp_it = qual_pts->cbegin(), qp_end=qual_pts->cend();
				 qp_it != qp_end; ++qp_it)
            {
                EvalPt* ep_ptr = eval_pts_rhs[ep_ctr++];
                f_rhs[ctr] += qp_it->weight() * ep_ptr->get_potential();
                h_rhs[ctr] += qp_it->weight()
							  * ep_ptr->get_field().dot(qp_it->normal());
                delete ep_ptr;
            }
            f_rhs[ctr] *= ONE_OVER_4PI / Dsolvent;
            h_rhs[ctr] *= ONE_OVER_4PI / Dsolvent;
#endif  // CENTROID_COLLOCATION

            //std::cout << ctr << ": " << f_rhs[ctr] << " " << h_rhs[ctr] << " " << "\n";
            ++ctr;
        }
    }
    eval_pts_rhs.clear(); // no longer contains valid pointers - deleted above
    std::cout << "RHS took: " << (myclock() - start_rhs_clock) / 1000.
	          << " ms" << std::endl;
    
    // the FMM has some internal profiling we can nab
    vanilla_fmm_timer += rhs_octree->get_timing_info();
    
    return (myclock() - start_rhs_clock); // return number of clock ticks taken
}

unsigned int BEEP::gmres(double residual_norm_stop_criteria, int max_iterations)
{
    // reset the fmm profile counters
    fmm->get_timing_info().zero();
    
    // DEBUG
    std::cout << "Starting GMRES." << std::endl;
    long start_gmres = myclock();
    long total_iter_timer = 0;

    unsigned int num_patches = get_total_patches();

    const int m=max_iterations; // retries
    const int n=2*num_patches; // dimensions of the problem
    boost::scoped_array<double> _b(new double[n]); // rhs
    boost::scoped_array<double> _x(new double[n]); // guess
    double* b = _b.get();
    double* x = _x.get();
    
    boost::scoped_array<double> preconditioner_rhs(new double[n]);
    boost::scoped_array<double> preconditioner_lhs(new double[n]);
    memset(preconditioner_rhs.get(),0,sizeof(double)*n);
    memset(preconditioner_lhs.get(),0,sizeof(double)*n);

	// stop when residual smaller than this multiple of rhs norm
    const double eps = residual_norm_stop_criteria;

    const bool detailed = true;
    //bool non_zero = false;

    // PREPARE RHS
    //std::cout << "RHS:\n";
    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        b[ctr] = f_rhs[ctr];
        b[ctr+num_patches] = h_rhs[ctr];

        //std::cout << f_rhs[ctr] << " " << h_rhs[ctr] << "\n";
    }

    // RESCALE?
    const double fscale = 1.0;
    const double hscale = 1.0;
    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        b[ctr] /= fscale;
        b[ctr+num_patches] /=hscale;
    }
    // same thing using BLAS
    //cblas_dscal(numCharges*2,1./fscale,&(b[0]),2);
    //cblas_dscal(numCharges*2 - 1,1./hscale,&(b[1]),2);

    // INIT INITIAL VALUES to what was set in input fh files
	bool preconditioned = meshes.init_library_fh_vals(x, num_patches);
//std::cout << "BEEP::gmres preconditioned " << preconditioned << std::endl;

#ifdef PRECONDITION 
std::cout << "BEEP::gmres PRECONDITION defined" << std::endl;
    // copy the x values into the preconditioner
    memcpy(preconditioner_lhs.get(), x, sizeof(double)*n);
    
    // If the preset fh vals are all zeroes, then it's not a very good
    // initial guess, would be better off with the rhs vector...
    if (preconditioned) {
        
        // pre-conditioning iteration
        do_bem_fmm_iteration(num_patches, preconditioner_lhs.get(),
		                                  preconditioner_rhs.get());
        
        // subtract the preconditioner_rhs from the actual rhs
        for (unsigned int np_ctr=0; np_ctr < num_patches; ++np_ctr) {
            b[np_ctr] -= preconditioner_rhs[np_ctr];
            b[np_ctr+num_patches] -= preconditioner_rhs[np_ctr+num_patches];
        }
    }
    
    // set lhs initial guess to rhs
    memcpy(x,b,sizeof(double)*n);
#endif  // PRECONDITION

    // smart pointers
    boost::scoped_array<double> _V(new double[n*(m+1)]);
    boost::scoped_array<double> _U(new double[m*(m+1)/2]);
    boost::scoped_array<double> _r(new double[n]);
    boost::scoped_array<double> _y(new double[m+1]);
    boost::scoped_array<double> _c(new double[m]);
    boost::scoped_array<double> _s(new double[m]);
    boost::scoped_array<double*> _v(new double*[m+1]);
    
    double* V = _V.get();
    double* U = _U.get();
    memset(V, 0, sizeof(double)*(n*(m+1)));
    memset(U, 0, sizeof(double)*(m*(m+1)/2));
    double* r = _r.get();
    double* y = _y.get();
    double* c = _c.get();
    double* s = _s.get();
    double** v = _v.get();
    
    for ( int i=0; i<=m; ++i ) v[i]=V+i*n;

	// Used to avoid incorrect nan results
	auto NEAR0 = [](double l) -> double { return (l==0?DBL_MIN:l); };

    int its=-1;
    {  // New scope -- this section appears to be imported - untouched
        double gmres_beta, h, rd, dd, nrm2b;
        int j, io, uij, u0j;
        nrm2b=cblas_dnrm2(n,b,1);

        io=0;
        do  { // "aussere Iteration
        ++io;
        std::cout << std::flush;
        total_iter_timer += do_bem_fmm_iteration(num_patches, x, r);
        std::cout << std::flush;

        cblas_daxpy(n,-1.,b,1,r,1);
        gmres_beta=cblas_dnrm2(n,r,1);
        cblas_dcopy(n,r,1,v[0],1);
        cblas_dscal(n,1./NEAR0(gmres_beta),v[0],1);

        y[0]=gmres_beta;
        j=0;
        uij=0;
        do { // innere Iteration j=0,...,m-1

            u0j=uij;

            double *xxx = v[j];
            double *rrr = v[j+1];

            total_iter_timer += do_bem_fmm_iteration(num_patches, xxx, rrr);

            //mult(v[j],v[j+1]);
            cblas_dgemv(CblasColMajor,CblasTrans,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);
            cblas_dgemv(CblasColMajor,CblasNoTrans,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);
            h=cblas_dnrm2(n,v[j+1],1);
            cblas_dscal(n,1./NEAR0(h),v[j+1],1);
            for ( int i=0; i<j; ++i ) { // rotiere neue Spalte
            double tmp = c[i]*U[uij]-s[i]*U[uij+1];
            U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
            U[uij]     = tmp;
            ++uij;
            }
            { // berechne neue Rotation
            rd     = U[uij];
            dd     = sqrt(rd*rd+h*h);
            c[j]   = rd/NEAR0(dd);
            s[j]   = -h/NEAR0(dd);
            U[uij] = dd;
            ++uij;
            }
            { // rotiere rechte Seite y (vorher: y[j+1]=0)
            y[j+1] = s[j]*y[j];
            y[j]   = c[j]*y[j];
            }
            ++j;
            if ( detailed )
            cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<y[j]<<"\t" << nrm2b*eps <<endl;

        } while ( j<m && fabs(y[j])>eps*nrm2b );  // >= precludes 0 = 0
        { // minimiere bzgl Y
			if (cblas_dasum(uij, U, 1) == 0 && cblas_dasum(j ,y, 1) == 0)
				;// Solve Uy'=y where U,y=0 ... continuity as y->0 => y'=y
			else
				cblas_dtpsv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,j,U,y,1);
            // korrigiere X
            cblas_dgemv(CblasColMajor,CblasNoTrans,n,j,-1.,V,n,y,1,1.,x,1);
        }
        //  } while ( fabs(y[j])>=eps*nrm2b );
        } while (false);

        // R"uckgabe: Zahl der inneren Iterationen
        its = m*(io-1)+j;
    }  // End scope

    // DEBUG
    std::cout << "GMRES reached convergence in " << its << " iterations."
	          << std::endl;

#ifdef PRECONDITION
    // Add preconditioned values back in
    if (preconditioned) {
        
        double total_E = 0.0;
        unsigned int offset=0;
        for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
				mit=meshes.cbegin(), mend=meshes.cend();
			 mit != mend; ++mit)
        {
            MeshInstance& m = **mit;
            if (m.isSilent()) continue;
            total_E += m.calculate_energy
								(kappa, &x[offset], &x[offset+num_patches]);
            offset += m.get_num_node_patches();
        }
        std::cout << "Total change in energy relative to preconditioner: "
		          << total_E << std::endl;
        
        // add the preconditioned values back on
        for (unsigned int np_ctr=0; np_ctr < num_patches; ++np_ctr) {
            x[np_ctr] += preconditioner_lhs[np_ctr];
            x[np_ctr+num_patches] += preconditioner_lhs[np_ctr+num_patches];
        }
    }
#endif  // PRECONDITION

    // get the results
    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        f_lhs[ctr] = x[ctr];
        h_lhs[ctr] = x[ctr+num_patches];
    }

    // normalise FMM timings
    fmm::TimeInfo& profiler = fmm->get_timing_info();
    profiler /= its+1;
    
    std::cout << "(actual GMRES algo took total of: "
	          << (myclock() - start_gmres - total_iter_timer)/1000. << " ms"
	          << std::endl;
    
    return its+1;
}

long BEEP::do_bem_fmm_iteration(
	unsigned int num_patches,
	double lhs[],
	double result[])
{
    long start = myclock();

    // get lhs from input
    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        f_lhs[ctr] = lhs[ctr];
        h_lhs[ctr] = lhs[ctr+num_patches];
    }

    // allocate memory for the explicit results
    boost::scoped_array<double> explicit_results_f(new double[num_patches]);
    boost::scoped_array<double> explicit_results_h(new double[num_patches]);
    memset(explicit_results_f.get(), 0, sizeof(double)*num_patches);
    memset(explicit_results_h.get(), 0, sizeof(double)*num_patches);

    // allocate memory for the FMM results
    boost::scoped_array<double> fmm_results_f(new double[num_patches]);
    boost::scoped_array<double> fmm_results_h(new double[num_patches]);
    memset(fmm_results_f.get(), 0, sizeof(double)*num_patches);
    memset(fmm_results_h.get(), 0, sizeof(double)*num_patches);

#ifdef CACHE_GPU_RESULTS

    // If the results are cached, do the fmm first to maximise the amount of
    // useful CPU work before we idle to wait on the GPU results
    // -- which will only need waiting for the first time round
    fmm->solve(kappa, beta0, f_lhs.get(), h_lhs.get(),
	           fmm_results_f.get(), fmm_results_h.get());
    wait_for_local_interactions();
    evaluate_local_neighbours
		(explicit_results_f.get(), explicit_results_h.get());

#else  // CACHE_GPU_RESULTS

    // if the results are not cached, need to kick them off on the GPU asap,
    // then do CPU fmm work in the foreground, then wait for GPU to complete.
    evaluate_local_neighbours
		(explicit_results_f.get(), explicit_results_h.get());
    fmm->solve(kappa, beta0, f_lhs.get(), h_lhs.get(),
	           fmm_results_f.get(), fmm_results_h.get());
    wait_for_local_interactions();

#endif  // CACHE_GPU_RESULTS

    // write rhs back
    //std::cout << "Matrix-Vector Result:\n";
    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        result[ctr] = explicit_results_f[ctr] + fmm_results_f[ctr];
        result[ctr + num_patches]
			= explicit_results_h[ctr] + fmm_results_h[ctr];
    }

#if 0 
    // statistics for explicit vs. fmm relative proportions --> total accuracy
    double min_prop_f=1e99;
    double max_prop_f=0;
    double total_prop_f=0;
    for (unsigned int ctr=0; ctr < num_patches; ++ctr)
    {
        double tmp = fabs(explicit_results_f[ctr]) + fabs(fmm_results_f[ctr]);
        double prop_f = fabs(fmm_results_f[ctr]) / tmp;
        total_prop_f += prop_f;
        if (prop_f > 0.5) {
            std::cout << "(" << ctr << ") Large prop_f=" << prop_f << ": (fmm/explicit/total) " << fmm_results_f[ctr] << " " << explicit_results_f[ctr] << " " << fmm_results_f[ctr] + explicit_results_f[ctr] << "\n";
        }
        min_prop_f = (prop_f < min_prop_f) ? prop_f : min_prop_f;
        max_prop_f = (prop_f > max_prop_f) ? prop_f : max_prop_f;
    }
    double mean_prop_f = total_prop_f / num_patches;
    std::cout << "f proportions min/max/ave: " << min_prop_f << " " << max_prop_f << " " << mean_prop_f << "\n";
    
    double min_prop_h=1e99;
    double max_prop_h=0;
    double total_prop_h=0;
    for (unsigned int ctr=0; ctr < num_patches; ++ctr)
    {
        double tmp = fabs(explicit_results_h[ctr]) + fabs(fmm_results_h[ctr]);
        double prop_h = fabs(fmm_results_h[ctr]) / tmp;
        total_prop_h += prop_h;
        if (prop_h > 0.5) {
            std::cout << "(" << ctr << ") Large prop_h=" << prop_h << ": (fmm/explicit/total) " << fmm_results_h[ctr] << " " << explicit_results_h[ctr] << " " << fmm_results_h[ctr] + explicit_results_h[ctr] << "\n";
        }
        min_prop_h = (prop_h < min_prop_h) ? prop_h : min_prop_h;
        max_prop_h = (prop_h > max_prop_h) ? prop_h : max_prop_h;
    }
    double mean_prop_h = total_prop_h / num_patches;
    std::cout << "h proportions min/max/ave: " << min_prop_h << " " << max_prop_h << " " << mean_prop_h << "\n";
    
#endif // if 0

    std::cout << "Time for BEM iteration: " << (myclock() - start) / 1000.
	          << " ms " << std::endl;

    return (myclock() - start);
}

long BEEP::init_fmm_octree() {
    // initialise an FMM octree for the mesh node patches -- subsequent calls
    // need to set the charge magnitudes at the patch locations (4 values)
    // then solve the FMM twice (for kappa=[smallnumber, approx. 0] and
    // kappa=[actual kappa value])
    // integrate over mesh

    long start_init_fmm = myclock();

    // first get the bounding cube for the mesh
    Vector centre = bounding_cube_centre;
    double edge_length = bounding_cube_edge_length;
    std::cout << centre << " " << edge_length << std::endl;

    // set the neighbourhood size
    unsigned int nbsize = (cmdline_bem_nbsize == -1)
		? DEFAULT_BEM_NEIGHBOURHOOD_SIZE : cmdline_bem_nbsize;
    
    fmm.reset(new FMM_BEM_Octree(nbsize, centre, edge_length));
    
    // for BEM-FMM we optimize the tree to have average number of neighbours,
    // so we relax the subdivision policy on octree nodes.
    // The usual octree divides an octree node as soon as the contents
    // reach the max_items_per_node limit.
    // NB: since we unset_strict_divisions we *must*
    // call optimize_neighbourhoods once we're done inserting
    // items.
    fmm->unset_strict_divisions();
    
    // Hack to limit the maximum octree depth.
    //fmm->set_max_depth(0);

    // insert patches (universe coords) into the octree
    for (MeshInstanceList::iterator mit=meshes.begin(), mend=meshes.end();
		 mit!=mend; ++mit)
    {
        MeshInstance& minst = **mit;
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::iterator
				np_it = minst.get_node_patches().begin(),
				np_end = minst.get_node_patches().end();
			 np_it != np_end; ++np_it)
        {
            // The FMM/BEM can accept shared pointers to BasicNodePatch objects
            // (since BasicNodePatch has a Vector-type interface)
            // Note: the fmm/bem tree will stash a copy of the shared pointer
            // so the stored objects will remain valid for the fmm tree lifetime
            fmm->add(*np_it);
        }
    }

    // fmm optimizations
    fmm->build_neighbourhoods();
    fmm->optimize_neighbourhoods();
    fmm->remove_empty_nodes();

    return (myclock() - start_init_fmm);
}

void BEEP::set_bounding_cube(
	const Vector& new_bc_centre,
	double new_edge_length)
{
    bounding_cube_centre = new_bc_centre;
    bounding_cube_edge_length = new_edge_length;
    manually_set_bounding_cube = true;
}

void BEEP::calculate_bounding_cube() {
    if (manually_set_bounding_cube == true) return;

    Vector max(-1e99,-1e99,-1e99); // this will be the 'top right' corner
    Vector min(1e99,1e99,1e99); // and this will be 'bottom left'
    Vector v;

    // loop over vertices in the mesh

    assert(meshes.size() > 0);
    for (MeshInstanceList::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit!=mend; ++mit)
    {
        const MeshInstance& m = **mit;
		double radius = m.get_radius();

        v = m.get_xyz_offset() + Vector(radius, radius, radius);
        max.x = v.x > max.x ? v.x : max.x;
        max.y = v.y > max.y ? v.y : max.y;
        max.z = v.z > max.z ? v.z : max.z;

        min.x = v.x < min.x ? v.x : min.x;
        min.y = v.y < min.y ? v.y : min.y;
        min.z = v.z < min.z ? v.z : min.z;

        v = m.get_xyz_offset() - Vector(radius, radius, radius);
        max.x = v.x > max.x ? v.x : max.x;
        max.y = v.y > max.y ? v.y : max.y;
        max.z = v.z > max.z ? v.z : max.z;

        min.x = v.x < min.x ? v.x : min.x;
        min.y = v.y < min.y ? v.y : min.y;
        min.z = v.z < min.z ? v.z : min.z;
    }

    // figure out the maximum edge length in x/y/z dimension
    Vector diff = max - min;
    bounding_cube_edge_length = diff.y > diff.x ? diff.y : diff.x;
    bounding_cube_edge_length = (diff.z > bounding_cube_edge_length)
		? diff.z : bounding_cube_edge_length;
    assert(bounding_cube_edge_length > 0.0);

    // centre of the cube is the mid point of the two extremities
    bounding_cube_centre = (max + min) / 2.0;
}

double BEEP::get_potential_at_point(const Vector& pt) const {
    bool inside_mesh = false;
    double potential = 0.0; // the return value

    // first need to figure out if the point is within a mesh instance
    for (MeshInstanceList::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        const MeshInstance& minst = **mit;
        if (minst.pt_is_internal(pt)) {
            inside_mesh = true;
            potential = minst.get_potential_at_internal_pt(pt);
            break;
        }
    }
    
    // if the above loop found that the point was within a mesh, then the
    // potential will now have been calculated and we can return
    if (inside_mesh) return potential;

    // otherwise we must calculate the total potential by integrating all
    // surface solutions
    for (MeshInstanceList::const_iterator
			mit=meshes.cbegin(), mend=meshes.cend();
		 mit != mend; ++mit)
    {
        const MeshInstance& minst = **mit;
        potential += minst.get_potential_at_external_pt(pt, kappa);
    }

    return potential;
}

void BEEP::write_opendx_xyz(
	const std::string& filename, 
	unsigned int x, 
	unsigned int y, 
	unsigned int z,
	const Vector& edges,
	const Vector& centre) const
{
    GridParms grid;
    grid.x_pts = x;
    grid.y_pts = y;
    grid.z_pts = z;
    grid.x_angstroms = edges.x;
    grid.y_angstroms = edges.y;
    grid.z_angstroms = edges.z;
    grid.origin_x = centre.x - (grid.x_angstroms/2.);
    grid.origin_y = centre.y - (grid.y_angstroms/2.);
    grid.origin_z = centre.z - (grid.z_angstroms/2.);

    write_opendx_grid(filename, grid);
}

void BEEP::write_opendx(
	const std::string& filename,
	unsigned int x,
	unsigned int y,
	unsigned int z) const
{
    GridParms grid;
    grid.x_pts = x;
    grid.y_pts = y;
    grid.z_pts = z;
    grid.x_angstroms = bounding_cube_edge_length;
    grid.y_angstroms = bounding_cube_edge_length;
    grid.z_angstroms = bounding_cube_edge_length;
    grid.origin_x = bounding_cube_centre.x - (grid.x_angstroms/2.);
    grid.origin_y = bounding_cube_centre.y - (grid.y_angstroms/2.);
    grid.origin_z = bounding_cube_centre.z - (grid.z_angstroms/2.);
    
    write_opendx_grid(filename, grid);
}

void BEEP::write_opendx_grid(
	const std::string& filename,
	const GridParms& grid) const
{
    
    std::ofstream buf(filename.c_str());
    buf << grid.OpenDX_Preamble();

    unsigned int x = grid.x_pts;
    unsigned int y = grid.y_pts;
    unsigned int z = grid.z_pts;
    
    boost::scoped_array<double> vals(new double[z]);
    
    int ii=0;
    for (unsigned int xctr=0; xctr < x; ++xctr) {
        for (unsigned int yctr=0; yctr < y; ++yctr) {
            #pragma omp parallel for
            for (unsigned int zctr=0; zctr < z; ++zctr) {
                Vector xyz = grid.GridIdxToPos(xctr, yctr, zctr);
                double pot = get_potential_at_point(xyz)
							 * beep_constants::convert_potential_to_kT_per_e;
                vals[zctr] = pot;
            }

            for (unsigned int zctr=0; zctr < z; ++zctr) {
                buf << boost::format("% 13.6e ") % vals[zctr]; 
                //buf << vals[zctr] << " ";
                if (++ii == 3) { ii=0; buf << "\n";  }
            }
        }
    }
    buf << grid.OpenDX_Suffix();
    buf.close();
}

