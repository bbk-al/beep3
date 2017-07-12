/*      Author: david fallaize    Created on: 21 Jul 2010 */

/*! \file mesh.cpp
 * \brief This module implements the reference mesh class.
 */

#include "../common/math_vector.h"
#include "mesh.h"
#include "bem_kernels.h"
#include <gsl/gsl_linalg.h>
#include "../fmm/fmm_octree.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <cassert>
#include <exception>
#include "constants.h"
#include "gts_utils.h"
#include "off_utils.h"
#include <boost/scoped_ptr.hpp>
#include <boost/functional/hash.hpp>

namespace fs = boost::filesystem;	// Easier to swap to std::filesystem in '17

// NB PRETRIPTR this would have copied triangle_ptrs vector pointing to same
// triangles;  so here the shared_ptrs are copied with same result but safer
// on destruction.  It is not clear why a Mesh copy would ever be needed.
Mesh::Mesh(const Mesh& other) {
    vertices.reserve(other.vertices.size());
    triangles.reserve(other.triangles.size());
    node_patches.reserve(other.node_patches.size());
    charges.reserve(other.charges.size());

	// NB2 Then all this would have happened by default? So delete constructor?
	// Which would also allow a default move constructor with !PRETRIPTR?
    vertices.insert
		(vertices.end(), other.vertices.begin(), other.vertices.end());
    triangles.insert
		(triangles.end(), other.triangles.begin(), other.triangles.end());
    node_patches.insert(node_patches.end(),
					    other.node_patches.begin(), other.node_patches.end());
    charges.insert(charges.end(), other.charges.begin(), other.charges.end());
#ifndef PREHYDROPHOBIC
    allCharges.insert(allCharges.end(),
					  other.allCharges.begin(), other.allCharges.end());
    allChargesMap.insert(other.allChargesMap.cbegin(),
						 other.allChargesMap.cend());
	nppc.insert(nppc.end(), other.nppc.begin(), other.nppc.end());
#endif // PREHYDROPHOBIC

    net_charge = other.net_charge;
    centre = other.centre;
    radius = other.radius;
    done_energy_precalcs = other.done_energy_precalcs;
}

#ifdef PRETRIPTR
Mesh::~Mesh() {
    // delete the dynamically allocated triangles
    for (size_t ii=0; ii < triangle_ptrs.size(); ++ii) {
        delete triangle_ptrs[ii];
    }
}
#endif // PRETRIPTR

Mesh::Mesh(const std::string& mesh_filename,
           const std::string& xyzqr_filename,
           bool force_planar
		  ) : done_energy_precalcs(false)		
{
	fs::path mesh_path{mesh_filename};
    init_mesh(mesh_path);
	fs::path xyzqr_path{xyzqr_filename};
    init_charges(xyzqr_path);
    
    // default to charge centroid...
#ifdef DELETED  // This is already done in init_charges and this does nothing!
    calculate_charges_centroid();
#endif
    
    if (force_planar) {
        create_bezier_triangles<Triangle>();
    }
    else {
        create_bezier_triangles<PNG1_Triangle>();
    }
 
    create_node_patches();

    // set maximum radius of any surface point from the centre
    radius = calculate_radius();

    init_energy_precalcs();
}

void Mesh::init(const std::string& mesh_filename,
                bool force_planar)
{
    // decompress the mesh tarball to temp folder
    try {
		MeshTarball mesh_tar{mesh_filename};

        // get the paths of the extracted files
        fs::path mesh_filename = mesh_tar.get_mesh_filename();
        fs::path xyzqr_filename = mesh_tar.get_xyzqr_filename();
        fs::path centre_filename = mesh_tar.get_centre_filename();
        fs::path energies_filename = mesh_tar.get_energies_filename();

        // init the mesh using the pre-calculated files (nodes and quads etc.)
        init_mesh(mesh_filename);
        init_charges(xyzqr_filename);
        init_centre(centre_filename);

        //std::vector<HybridBezierTriangle> bezier_patches;
        if (force_planar) {
            create_bezier_triangles<Triangle>();
        }
        else {
            create_bezier_triangles<PNG1_Triangle>();
        }
        create_node_patches();
	
    	// read energy precalcs from data file
        if (force_planar) {
             init_energy_precalcs();
        }
        else {
            read_energy_precalcs(energies_filename);    
        }
        
        // read the fh values (if they are defined in the mesh-tar-gzip file)
        try {
            // this will throw an exception if the fh_file doesn't exist
            fs::path fh_filename = mesh_tar.get_fh_filename();
            
            // if haven't thrown an exception by this point then
            // can load fh values from file
            init_fh_vals(fh_filename);
        }
        catch (MeshTarball::MeshTarball_Exception) {
            // no FH values in the .mtz file
            // so nothing to do!
        }

        // calculate the maximum radius of any surface point from centre
        radius = calculate_radius();

    }
    catch (MeshTarball::MeshTarball_Exception) {
        throw MeshError();
    }

//	std::stringstream buf;
//	buf << "@kinemage\n";
//	for (std::vector<HybridBezierTriangle>::const_iterator it=bezier_patches.begin(), end=bezier_patches.end(); it != end; ++it)
//	{
//		const HybridBezierTriangle& bez = *it;
//		bez.kinemage(buf);
//	}
//	std::ofstream fout("mesh.kin");
//	fout << buf.str();
//	fout.close();
}

std::vector<Charge> Mesh::get_ecm_charges(
	const std::string& ecm_filename,
	const Quaternion& rot_to_universe,
	const Vector& xyz_offset) const
{
    std::vector<Charge> ecm_charges;
    
    Charge::read_charges_from_file(ecm_filename, ecm_charges);
    for (std::vector<Charge>::iterator
			it=ecm_charges.begin(), end=ecm_charges.end();
		 it != end; ++it)
    {
        it->change_coordinate_frame(centre, rot_to_universe, xyz_offset);
    }
    
	return ecm_charges;
}

void Mesh::init_centre(const fs::path& centre_filename) {
    // open the file
    fs::ifstream centre_file;
    centre_file.open(centre_filename);
    assert(centre_file.good());

    centre_file >> centre.x >> centre.y >> centre.z;
    //std::cout << "Read centre: " << centre << std::endl;

    centre_file.close();
}

void Mesh::init_mesh(const fs::path& mesh_filename) {
    
	if (mesh_filename.extension() == ".off") {
        // Call a utility function to read the actual file
        OFF_Surface::read_mesh_from_file(mesh_filename.string(),
                                         vertices,
                                         triangles);
        
    } 
	else if (mesh_filename.extension() == ".gts") {
        // Call a utility function to read the actual file
        GTS_Surface::read_mesh_from_file(mesh_filename.string(),
                                         vertices,
                                         triangles);
    }
    else {
        std::cerr << "Mesh filename " << mesh_filename
		          << " does not have either a GTS or OFF extension."
		          << std::endl;
        throw MeshError();
    }
    
    // set the total planar area (flat triangles)
    total_planar_area = 0;
    for (std::vector<Triangle>::const_iterator
			it=triangles.cbegin(), end=triangles.cend();
         it != end; ++it)
    {
        total_planar_area += it->get_planar_area();
    }
    
    // generally useful for vertices to have coherent vertex indices
    reorder_vertex_triangles();
    
    return;
}

void Mesh::init_charges(const fs::path& xyzqr_filename) {
    // read the charges using utility function in Charge class
#ifdef PREHYDROPHOBIC
    net_charge = Charge::read_charges_from_file(xyzqr_filename, charges);
#else
    net_charge = Charge::read_charges_from_file(xyzqr_filename, allCharges,
												false);
	// Build a key map to find the charge indices later and
	// copy the non-zero ones into charges
	unsigned int ctr = 0;
	for (std::vector<Charge>::const_iterator
			ch = allCharges.cbegin(), ch_end = allCharges.cend();
		 ch != ch_end; ++ch)
	{
		// Despite its name, exclude charges with no charge, hydrophobicity or
		// LJ parameters - these will be HETATMs such as water
		if (ch->charge != 0 || ch->hydrophobicity != 0 ||
			(ch->sigma != 0 && ch->epsilon != 0)) {
			allChargesMap[*ch] = ctr++;

			// And for backward compatibility include only charges in this list
			if (ch->charge != 0) {
				Charge c(*ch);
				charges.push_back(c);
			}
		}
	}
	
#endif // PREHYDROPHOBIC

    //std::cout << "Added " << charges.size() << " charges." << std::endl;
    centre = calculate_charges_centroid();
    //std::cout << centre << std::endl;
}

double Mesh::calculate_radius() {
    radius = 0.0;

    for (std::vector<BasicNodePatch>::const_iterator
			np_it=node_patches.cbegin(), np_end=node_patches.cend();
         np_it != np_end; ++np_it)
    {
        double dist = (*np_it - centre).length();
        radius = (dist > radius) ? dist : radius;
    }

    return radius;
}

double Mesh::calculate_volume() const {
    double total_volume = 0.0;
    for (std::vector<Triangle>::const_iterator
			tri_it = triangles.cbegin(), tri_end=triangles.cend();
		 tri_it != tri_end; ++tri_it)
    {
        Vector P = tri_it->v1 - centre;
        Vector Q = tri_it->v2 - centre;
        Vector R = tri_it->v3 - centre;
        
        double pv =  P.x*Q.y*R.z + P.y*Q.z*R.x + P.z*Q.x*R.y
		           - P.x*Q.z*R.y - P.y*Q.x*R.z - P.z*Q.y*R.x;
        total_volume += pv;
    }
    return total_volume / 6.;
}

#if 0
void Mesh::init_ellipsoid_shape_def(const std::string& ellipsoid_filename)
{
    // open the file
    std::ifstream ellipsoid_file;
    ellipsoid_file.open(ellipsoid_filename.c_str());
    assert(ellipsoid_file.good());

    double x,y,z,r;
    ellipsoid_file >> x >> y >> z;
    ellipsoid_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    Vector centre(x,y,z);

    Vector radii(0,0,0);

    ellipsoid_file >> x >> y >> z >> r;
    ellipsoid_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    Vector v1(x,y,z);
    radii.x = r;

    ellipsoid_file >> x >> y >> z >> r;
    ellipsoid_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    Vector v2(x,y,z);
    radii.y = r;

    ellipsoid_file >> x >> y >> z >> r;
    ellipsoid_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    Vector v3(x,y,z);
    radii.z = r;

    shape_defs.push_back(EllipsoidShapeDefinition(centre, v1, v2, v3, radii));

    ellipsoid_file.close();

}

Vector Mesh::calculate_planar_E_from_spharms(const Vector& v, const Vector& origin, const SpharmHolder& spf) const
{
    Vector planar_E(0,0,0);

    planar_E = planar_E - spf.dphi(v, origin);
    planar_E = planar_E - spf.dtheta(v, origin);

    return planar_E;
}

void Mesh::write_fh(const std::string& filename)
{
    std::ofstream fh_vals;
    fh_vals.open(filename.c_str(), std::ios::out);
    fh_vals << std::setprecision(10);

    unsigned int ctr=0;
    for (std::vector<BasicNodePatch>::const_iterator it=node_patches.begin(), end=node_patches.end();
        it != end;
        ++it)
    {
        const BasicNodePatch& np = *it;
        fh_vals << np.f*4171.8 << " " << np.h *4171.8 << "\n";
    }
    fh_vals.close();

    return;
}

Vector Mesh::calc_boundary_force(const Offsets& offs, double Dprotein, double Dsolvent, double kappa, const SpharmHolder& spf, const SpharmHolder& sph) const
{
    Vector f(0,0,0);
    double epsilon = Dsolvent / Dprotein;

    // offs gives the offsets of node/triangle/charge indices which identify a given molecule within
    // the entire mesh

    unsigned int first_idx = offs.vertex_numbering_offset;
    unsigned int end_idx = first_idx + offs.num_vertices;

    for (unsigned int idx=first_idx; idx < end_idx; ++idx)
    {
        const BasicNodePatch& np = node_patches[idx];

        Vector E_planar = calculate_planar_E_from_spharms(np.node, offs.origin, spf);
        //double fval = spf.evaluate(np.node, offs.origin);
        //double hval = sph.evaluate(np.node, offs.origin);
        double hval = np.h;

        // TODO:: += operator for Vectors !!
        Vector Eout = E_planar + (np.normal * -hval);
        Vector Ein = E_planar + (np.normal * -hval*epsilon);

        //f = f + force_from_E(E, tri.normal * tri.area, kappa, meanf);
        f = f - np.normal * np.area * (Dsolvent - Dprotein) * Eout.dot(Ein);
//         f = f + np.normal * np.area * (np.normal.dot(force_from_E(Ein, np.normal, 0.0, fval)) -
//                                        np.normal.dot(force_from_E(Eout, np.normal, kappa, fval)));

    }
//     unsigned int first_idx = offs.triangle_numbering_offset;
//     unsigned int end_idx = first_idx + offs.num_triangles;
//
//     for (unsigned int idx=first_idx; idx < end_idx; ++idx)
//     {
//         const Triangle& tri = triangles[idx];
//         Vector E_planar = f_at_verts_to_E(
//         Vector E_planar = calculate_planar_E_from_spharms(np.node, offs.origin, spf);
//
//         double hval = np.h;
//
//         // TODO:: += operator for Vectors !!
//         Vector Eout = E_planar + (np.normal * -hval);
//         Vector Ein = E_planar + (np.normal * -hval*epsilon);
//
//         //f = f + force_from_E(E, tri.normal * tri.area, kappa, meanf);
//         f = f - np.normal * np.area * (80.0 - 2.0) * Eout.dot(Ein);
// //         f = f + np.normal * np.area * (np.normal.dot(force_from_E(Ein, np.normal, 0.0, fval)) -
// //                                        np.normal.dot(force_from_E(Eout, np.normal, kappa, fval)));
//
//     }
    return f / 2.0;

}

Vector Mesh::force_from_E(const Vector& E, const Vector& normal, double kappa, double f) const
{
    const double En = E.dot(normal);
    const double half_E2_kf2 = 0.5 * (E.dot(E) + (kappa*f)*(kappa*f));

    Vector force(0,0,0);
    force.x = E.x*En - half_E2_kf2*normal.x;
    force.y = E.y*En - half_E2_kf2*normal.y;
    force.z = E.z*En - half_E2_kf2*normal.z;

    //std::cout << "elemental force: " << force << std::endl;
    return force;
}

Vector Mesh::calculate_MST_force(const Offsets& offs, double Dsolvent, double kappa) const
{

    unsigned int first_idx = offs.triangle_numbering_offset;
    unsigned int end_idx = first_idx + offs.num_triangles;

//     for (int ii=0; ii < node_patches.size(); ++ii)
//     {
//         std::cout << "fht (ii=" << ii << ") " << node_patches[ii].alt_normal.length()*node_patches[ii].h << std::endl;
//     }

    int ctr=0;
    Vector F(0,0,0);
    for (unsigned int idx=first_idx; idx < end_idx; ++idx)
    {
        const Triangle& tri = triangles[idx];

        // Create a 3x3 matrix containing the geometry of the triangle
        double a_data[9];
        double rhs[3];

        const BasicNodePatch& np1 = node_patches[tri.v1_idx];
        const BasicNodePatch& np2 = node_patches[tri.v2_idx];
        const BasicNodePatch& np3 = node_patches[tri.v3_idx];
        const Vector& n = tri.normal;

        Vector alt_n = (np1.alt_normal.normalised() +
                        np2.alt_normal.normalised() +
                        np3.alt_normal.normalised()) / 3.0;
        double interp_H = (np1.alt_normal.length()*np1.h +
                        np2.alt_normal.length()*np2.h +
                        np3.alt_normal.length()*np3.h) / 3.0;
        Vector axis1 = np2.node - np1.node;
        Vector axis2 = np3.node - np1.node;
        a_data[0] = axis1.x;
        a_data[1] = axis1.y;
        a_data[2] = axis1.z;
        a_data[3] = axis2.x;
        a_data[4] = axis2.y;
        a_data[5] = axis2.z;
        a_data[6] = alt_n.x;
        a_data[7] = alt_n.y;
        a_data[8] = alt_n.z;
/*        a_data[0] = axis1.x;
        a_data[1] = axis2.x;
        a_data[2] = alt_n.x;
        a_data[3] = axis1.y;
        a_data[4] = axis2.y;
        a_data[5] = alt_n.y;
        a_data[6] = axis1.z;
        a_data[7] = axis2.z;
        a_data[8] = alt_n.z;*/

        //std::cout << "A(0,ii) " << ctr++ << " " << a_data[0] << " " << a_data[3] << " " << a_data[6] << std::endl;

        rhs[0] = np2.f - np1.f;
        rhs[1] = np3.f - np1.f;
        rhs[2] = interp_H;

        //std::cout << "rhs " << rhs[0] << " " << rhs[1] << " " << rhs[2] << std::endl;
        gsl_matrix_view m = gsl_matrix_view_array(a_data, 3, 3);
        gsl_vector_view b = gsl_vector_view_array(rhs, 3);
        gsl_vector *x = gsl_vector_alloc(3);
        int s;
        gsl_permutation * p = gsl_permutation_alloc(3);
        gsl_linalg_LU_decomp (&m.matrix, p, &s);
        gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
        Vector E(-gsl_vector_get(x,0),
                -gsl_vector_get(x,1),
                -gsl_vector_get(x,2));

        //std::cout << "E on triangle @ " << tri.centre << " #"  << ctr++ << " is: " << E << std::endl;

        double mean_f = (np1.f + np2.f + np3.f) / 3.0;
        Vector inc_F = force_from_E(E, tri.normal, kappa, mean_f)*tri.area;
        //std::cout << "Force on element: " << inc_F << " " << tri.area << "\n";
        F = F + inc_F;

        // free GSL allocated memory
        gsl_permutation_free(p);
        gsl_vector_free(x);

    }

    return F*Dsolvent;

}
#endif // if 0

void Mesh::calculate_vertex_normals() {

    // loop over all vertex objects in the mesh
    for (unsigned int vctr=0; vctr < vertices.size(); ++vctr) {
        Vertex& vertex = vertices[vctr];
        //BasicNodePatch& np = node_patches[vctr];

        double a_data[9];
        for (unsigned short ii=0; ii < 9; ++ii) {
            a_data[ii] = 0;
        }
        Vector total_n(0,0,0);

        for (std::vector<unsigned int>::const_iterator
				connect_it=vertex.triangle_indices.cbegin(),
				connect_end=vertex.triangle_indices.cend();
             connect_it != connect_end; ++connect_it)
        {
            const Triangle& tri = triangles[*connect_it];
            const Vector& n = tri.get_planar_normal();
            total_n = total_n + n;

            for (unsigned short ii=0; ii < 3; ++ii) {
                for (unsigned short jj=0; jj < 3; ++jj) {
                    a_data[ii*3 + jj] += n[ii]*n[jj];
                }
            }

        }
        double b_data[3] = { total_n.x, total_n.y, total_n.z };
        gsl_matrix_view m = gsl_matrix_view_array(a_data, 3, 3);
        gsl_vector_view b = gsl_vector_view_array(b_data, 3);
        gsl_vector *x = gsl_vector_alloc(3);
        int s;
        gsl_permutation * p = gsl_permutation_alloc(3);
        gsl_linalg_LU_decomp (&m.matrix, p, &s);
        double detA = gsl_linalg_LU_det(&m.matrix, s);

        if (detA > 1e-6) {
            gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

            // TODO:: check determinant of matrix before solving
            // if det(A) small then vector will be rubbish

            vertex.normal = Vector(gsl_vector_get(x,0),
                                gsl_vector_get(x,1),
                                gsl_vector_get(x,2));
        }
        else {
            vertex.normal = total_n / vertex.triangle_indices.size();
        }

        //std::cout << "Alt normal for patch " << ctr++ << ": " << np.alt_normal.normalised() << " " << np.normal << std::endl;

        gsl_permutation_free(p);
        gsl_vector_free(x);

        vertex.normal.normalise();

    }
}

void Mesh::init_fh_vals(const fs::path& fh_filename) {
  
    // open the file
    fs::ifstream fh_file;
    fh_file.open(fh_filename);
    assert(fh_file.good());

    for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
    {
        BasicNodePatch& np = *it;
        fh_file >> np.f >> np.h;
    }

    fh_file.close();
}

void Mesh::init_energy_precalcs() {

    const unsigned int num_patches = node_patches.size();
#ifdef OPENCL
    boost::scoped_ptr<OpenCL_Handler> opencl_handler(new OpenCL_Handler);
    const unsigned int NB_SIZE = 2500;
#else  // OPENCL
    const unsigned int NB_SIZE = 80;
#endif  // OPENCL

    double net_charge = 0.0;

    // first get the bounding cube for the mesh
    Vector centre;
    double edge_length;
    get_bounding_cube(centre, edge_length);
#ifndef PREHYDROPHOBIC
//TODO this should be in its own sub-method to avoid distraction

	// Octree for hydrophobicities
	// This achieves O(logN) for charges, but alas not for node patches.
	boost::scoped_ptr<Octree<Node<Charge>, Charge> >
		hytree(new Octree< Node<Charge>, Charge >(10, centre, edge_length));
    // add all charges into the octree
    for (std::vector<Charge>::const_iterator
			ch_it=allCharges.cbegin(), ch_end=allCharges.cend();
		 ch_it != ch_end; ++ch_it)
    {
        const Charge& ch = *ch_it;
		hytree->add(ch);
    }
	// Find each node's nearest charge and assign hydrophobicity to the patch
	auto dist = [](const Vector& v, const Charge& c) noexcept -> double
		{ return (v - c).length() - c.radius; };
    for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
    {
		// Do this relative to vertex;  arguably should be to centroid
		const Charge& c = hytree->get_nearest(it->get_node(), dist);
		it->ch_idx = allChargesMap.at(c);	// Exception if c not present
		if (it->ch_idx >= nppc.size()) nppc.resize(it->ch_idx+1);
		nppc[it->ch_idx]++;				// Number of NP's for this charge
		it->hydrophobicity = c.hydrophobicity;  // Stored on np for smoothing
    }
	// Smooth values across all node patches
	// TODO is this actually useful;  if not, lose np hydrophobicity
	std::vector<double> merged;
	unsigned int ctr = 0;
	merged.reserve(get_num_node_patches());
    for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
    {
		// Get list of indices of mesh triangles for which this node is a vertex
		const std::vector<unsigned int>& tri
			= vertices[it->vertex_idx].get_triangle_indices();
		// Obtain weighted average of neighbouring node patch hydrophobicities
		// Easiest to double count each node's contribution (triangles link)
		// And by construction, node patch indices match vertex indices...
		merged[ctr] = it->hydrophobicity;
		Charge& ch = allCharges[it->ch_idx];
		double datom = (ch - it->get_node()).length();
		double ratom = ch.radius;
		double sow = 1.0;
		if (datom > ratom) {
			sow = 2 * ratom / datom;	// count double as others are
			merged[ctr] *= sow;
			for (std::vector<unsigned int>::const_iterator
					t=tri.begin(), tend=tri.end();
				t != tend; ++t)
			{
				// Add up the hydrophobicities from the neighbouring nodes
				const Triangle& tr = triangles[*t];
				unsigned int tvi[] = {tr.v1_idx != ctr ? tr.v1_idx : tr.v2_idx,
									tr.v3_idx != ctr ? tr.v3_idx : tr.v2_idx};
				for (int i = 0; i < 2; i++) {
					BasicNodePatch& nbr = node_patches[tvi[i]];
					Charge& ch = allCharges[nbr.ch_idx];	// new reference
					datom = (ch - nbr.get_node()).length();
					ratom = ch.radius;
					double w = ratom /
						(datom + (nbr.get_node() - it->get_node()).length());
					merged[ctr] += nbr.hydrophobicity * w;
					sow += w;
				}
			}
		}
		merged[ctr] /= sow;  // sum of weights
		ctr++;
	}
	// Play the merged values back into the nodes
	ctr = 0;
    for (auto& np: node_patches) {
		np.hydrophobicity = merged[ctr];
		ctr++;
	}

	// Octree for FMM
#endif // PREHYDROPHOBIC
    boost::scoped_ptr<fmm::FMM_Octree_6FIG_ACCURACY>
		fmm(new fmm::FMM_Octree_6FIG_ACCURACY(NB_SIZE, centre, edge_length));

    // loop over each mesh instance and insert the charges into a big fat octree
    for (std::vector<Charge>::const_iterator
			ch_it=charges.cbegin(), ch_end=charges.cend();
		 ch_it != ch_end; ++ch_it)
    {
        const Charge& ch = *ch_it;
        if (ch.get_charge() != 0.0) {
			fmm->add(ch);
            net_charge += ch.get_charge();
        }
    }

    // solve the FMM
    fmm->solve(0.0);

    // create list of evaluation points
    const size_t CHUNKSIZE = 10000; //64*1024;
    boost::shared_array<fmm::EvalPt_2ndDerivs>
		eval_points(new fmm::EvalPt_2ndDerivs[CHUNKSIZE]);
    assert(eval_points.get() != nullptr);
    std::vector<fmm::EvalPt_2ndDerivs*> ep_vec;
    ep_vec.reserve(CHUNKSIZE);

    unsigned int patch_ctr=0;
    std::vector<BasicNodePatch>::iterator
		it=node_patches.begin(), end=node_patches.end();
    double total_weight = 0.0;

    while (it != end) {
        ep_vec.clear();
        
        // this locks the energy quad points until we've finished with them
        // as long as these shared_ptrs exist, the underlying quadpoints won't
        // be deleted --> won't need recalculating for the other loop below
        std::vector<boost::shared_ptr<QuadList> > eqp_cache;
        
        std::vector<BasicNodePatch>::iterator result_it=it;
        size_t ep_ctr=0;
        boost::shared_ptr<QuadList> energy_quad_points
			= it->get_energy_quad_points();
        for (; it != end
				&& (ep_ctr + energy_quad_points->size()) < CHUNKSIZE;
			 ++it)
        {
			// prevent auto-deletion of quad points!
            eqp_cache.push_back(energy_quad_points);
            
            for (std::vector<QuadPoint>::const_iterator
					qp_it=energy_quad_points->cbegin(),
					qp_end=energy_quad_points->cend();
				 qp_it != qp_end; ++qp_it)
            {
                eval_points[ep_ctr].reset();
                eval_points[ep_ctr].pt() = Vector(qp_it->pt());
                ++ep_ctr;

				// avoid overrunning arrays
                assert(ep_ctr < CHUNKSIZE);
            }
            
            // have to do this here to provide peek-ahead of the
            // energy_quad_points list (the list of quad points is generated
            // on-demand, and only cached as long as the pointer remains in
            // existence, so need to prevent it dropping out of scope whilst
            // we're using it... Slightly annoying I know, but otherwise
            // we have to prestore all those pesky quad points and they
            // actually add up to quite a lot of memory y'know.)
            if (it+1 != end) {
                energy_quad_points = (it+1)->get_energy_quad_points();
            }
        }

        //ep_vec.resize(ep_ctr);
        for (size_t setter=0; setter < ep_ctr; ++setter) {
            ep_vec.push_back(&(eval_points[setter]));
        }

        // evaluate the FMM at the evaluation points
        // - results end up within the EvalPt struct
        // NB: this call may invoke some OpenCL stuff to speed it up
#ifdef OPENCL
        fmm->evaluate_many(ep_vec, *opencl_handler);
#else
        fmm->evaluate_many(ep_vec);
#endif  // OPENCL

        ep_ctr=0;
        for (; result_it != it; ++result_it) {
            boost::shared_ptr<QuadList> quad_points
				= result_it->get_energy_quad_points();
            result_it->energy_coefficient_f = 0;
            result_it->energy_coefficient_h = 0;
            KahanVector fres;
            KahanVector hres;
            for (unsigned int ii=0; ii < quad_points->size(); ++ii) {
                const QuadPoint& qp = (*quad_points)[ii];
                fmm::EvalPt_2ndDerivs& ep = eval_points[ep_ctr++];

                result_it->energy_coefficient_f
					+= ep.get_field().dot( qp.normal() ) * qp.weight()
						* ONE_OVER_4PI;
                result_it->energy_coefficient_h
					+= ep.get_potential() * qp.weight() * ONE_OVER_4PI;

                const GradField3x3& f2 = ep.get_field2();

				// f coeff
                double xx = Vector(f2(0,0), f2(0,1), f2(0,2))
							.dot( qp.normal() ) * qp.weight() * ONE_OVER_4PI;
                double yy = Vector(f2(1,0), f2(1,1), f2(1,2))
							.dot( qp.normal() ) * qp.weight() * ONE_OVER_4PI;
                double zz = Vector(f2(2,0), f2(2,1), f2(2,2))
							.dot( qp.normal() ) * qp.weight() * ONE_OVER_4PI;
                fres += Vector(xx,yy,zz);

				// h coeff
                hres += ep.get_field() * qp.weight() * ONE_OVER_4PI;

                total_weight += qp.weight();
            }

            result_it->force_coefficient_f = *fres;
            result_it->force_coefficient_h = *hres;
 
            // progress counter
            ++patch_ctr;
            //if (patch_ctr % 100 == 0 || patch_ctr == num_patches) { std::cout << patch_ctr << " / " << num_patches << "\r" << std::flush; }
        }
    }

    // renormalise the energy/force coefficients - enforce Gauss Law and
    // zero force summation
    renormalise_energy_precalcs();

    // Done precalcs
    done_energy_precalcs = true;
}
    
void Mesh::renormalise_energy_precalcs()
{
    // Gauss' Law says surface integral of E.n is equal to enclosed charge.
    // Since E.n is the f_energy_coefficient, then the sum total_f should equal
    // the net_charge calculated earlier. Adjust the values to achieve equality.
    double total_f=0;
    for (std::vector<BasicNodePatch>::const_iterator
			np_it=node_patches.cbegin(), np_end=node_patches.cend();
		 np_it != np_end; ++np_it)
    {
        const BasicNodePatch& np = *np_it;
        total_f -= np.energy_coefficient_f;
    }
    double gauss_error = net_charge - total_f;
    for (std::vector<BasicNodePatch>::iterator
			np_it=node_patches.begin(), np_end=node_patches.end();
		 np_it != np_end; ++np_it)
    {
        BasicNodePatch& np = *np_it;
        np.energy_coefficient_f
			-= (np.get_bezier_area() / total_bezier_area) * gauss_error;
    }
    
//     Vector zero_force(0,0,0);
//     for (std::vector<BasicNodePatch>::const_iterator np_it=node_patches.begin(), np_end=node_patches.end(); np_it != np_end; ++np_it)
//     {
//         const BasicNodePatch& np = *np_it;
//         zero_force += (np.force_coefficient_h*np.energy_coefficient_f - np.force_coefficient_f*np.energy_coefficient_h) / np.get_bezier_area();
//     }
//     for (std::vector<BasicNodePatch>::iterator np_it=node_patches.begin(), np_end=node_patches.end(); np_it != np_end; ++np_it)
//     {
//         BasicNodePatch& np = *np_it;
//         Vector half_err = zero_force * np.get_bezier_area() * np.get_bezier_area() / (2.0*total_bezier_area);
//         np.force_coefficient_h -= half_err / np.energy_coefficient_f;
//         np.force_coefficient_f += half_err / np.energy_coefficient_h;
//     }
}

void Mesh::read_energy_precalcs(
	const fs::path& energies_filename)
{
    // open the file
    fs::ifstream energies_file;
    energies_file.open(energies_filename);
    assert(energies_file.good());

	std::string str;
	std::istringstream is;
    for (std::vector<BasicNodePatch>::iterator
			it=node_patches.begin(), end=node_patches.end();
		 it != end; ++it)
    {
        BasicNodePatch& np = *it;
#ifndef PREHYDROPHOBIC
		// Default the values that might not be present
		np.hydrophobicity = 0.0;
		np.ch_idx = allCharges.size();
#endif // PREHYDROPHOBIC
		std::getline(energies_file, str);
		is.clear();	// Needed to reset after reading past end of string
		is.str(str);
        is >> np.energy_coefficient_f >> np.energy_coefficient_h 
           >> np.force_coefficient_f.x >> np.force_coefficient_f.y
		   >> np.force_coefficient_f.z 
           >> np.force_coefficient_h.x >> np.force_coefficient_h.y
#ifdef PREHYDROPHOBIC
		   >> np.force_coefficient_h.z;
#else
		   >> np.force_coefficient_h.z
		   >> np.hydrophobicity >> np.ch_idx;
		if (np.ch_idx >= nppc.size()) nppc.resize(np.ch_idx+1);
		nppc[np.ch_idx]++;				// Number of NP's for this charge
#endif // PREHYDROPHOBIC
    }

    energies_file.close();
    done_energy_precalcs = true;
}

KahanVector Mesh::calculate_qe_force(
	double Dprotein,
	double Dsolvent,
	const double fvals[],
	const double hvals[]) const
{
    assert(done_energy_precalcs == true);
    
    const double epsilon_ratio = Dsolvent / Dprotein;
    
    KahanVector qe_force; // this does compensated addition on the += operator
    unsigned int ctr=0;
    
    for (std::vector<BasicNodePatch>::const_iterator
			nit=node_patches.cbegin(), nend=node_patches.cend();
		 nit != nend; ++nit)
    {
        const BasicNodePatch& np = *nit;
        
        // compensated addition of the qe force component
        qe_force += np.force_coefficient_h*epsilon_ratio*hvals[ctr];
        qe_force -= np.force_coefficient_f*fvals[ctr];
        //std::cout << "forces: " << np.force_coefficient_f << " " << np.force_coefficient_h << " fh: " << fvals[ctr] << " " << hvals[ctr] << " qe: " << *qe_force << std::endl;
        ++ctr;
    }
    
    qe_force *= beep_constants::convert_force_to_kJ_per_mol_Angstrom;
    return qe_force;
}

void Mesh::calculate_surface_integral_forces(
	double kappa,
	double Dprotein,
	double Dsolvent,
	const double fvals[],
	const double hvals[],
	KahanVector& MST_external,
	KahanVector& MST_internal,
	KahanVector& dbf,
	KahanVector& ionic) const
{
    for (unsigned int tctr=0; tctr < triangles.size(); ++tctr) {
        // this may in fact be a PNG1, Bezier, or planar triangle
        // virtual function calls in calculate_force_components
        // will take advantage of the parametric representation of
        // each type, so don't need to worry here about what type
        // it actually is.
#ifdef PRETRIPTR
        const Triangle& tri = dynamic_cast<Triangle&>(*(triangle_ptrs[tctr])); 
#else
        const Triangle& tri = *(triangle_ptrs[tctr]);
#endif // PRETRIPTR

        const Vertex& v1 = vertices[tri.v1_idx];
        const Vector& n1 = v1.normal;
        double f1 = fvals[tri.v1_idx];
        double h1 = hvals[tri.v1_idx];

        const Vertex& v2 = vertices[tri.v2_idx];
        const Vertex& n2 = v2.normal;
        double f2 = fvals[tri.v2_idx];
        double h2 = hvals[tri.v2_idx];

        const Vertex& v3 = vertices[tri.v3_idx];
        const Vertex& n3 = v3.normal;
        double f3 = fvals[tri.v3_idx];
        double h3 = hvals[tri.v3_idx];
        
        //std::cout << "f/h vals over triangle: " << Vector(f1,f2,f3) << " " << Vector(h1,h2,h3) << std::endl;

        calculate_force_components(tri, Vector(f1,f2,f3), Vector(h1,h2,h3),
			Dprotein, Dsolvent, kappa, BasicTriangle::gauss_legendre_pts_1(),
			0, MST_external, MST_internal, dbf, ionic);
    }
    
    MST_external *= beep_constants::convert_force_to_kJ_per_mol_Angstrom;
    MST_internal *= beep_constants::convert_force_to_kJ_per_mol_Angstrom;
    dbf *= beep_constants::convert_force_to_kJ_per_mol_Angstrom;
    ionic *= beep_constants::convert_force_to_kJ_per_mol_Angstrom;
}

void Mesh::calculate_forces(
	double kappa,
	double Dprotein,
	double Dsolvent,
	Vector& qE,
	Vector& _MST_external,
	Vector& _MST_internal,
	Vector& _dbf,
	Vector& _ionic) const
{
    assert(node_patches.size() == vertices.size());

    KahanVector kMST_ext;
    KahanVector kMST_int;
    KahanVector kdbf;
    KahanVector kionic;
    
    double* fvals = new double[vertices.size()];
    double* hvals = new double[vertices.size()];

    size_t ctr=0;
    for (std::vector<BasicNodePatch>::const_iterator
			it=node_patches.cbegin(), end=node_patches.cend();
		 it != end; ++it)
    {
        const BasicNodePatch& np = *it;
        fvals[ctr] = np.f;
        hvals[ctr] = np.h;
        ++ctr;
    }

    qE = *(calculate_qe_force(Dprotein, Dsolvent, fvals, hvals));
    calculate_surface_integral_forces(kappa, Dprotein, Dsolvent, fvals, hvals,
	                                  kMST_ext, kMST_int, kdbf, kionic);

    _MST_external = *kMST_ext;
    _MST_internal = *kMST_int;
    _dbf = *kdbf;
    _ionic = *kionic;
    
    delete[] fvals;
    delete[] hvals;
}

double Mesh::calculate_energy(
	double kappa,
	double Dprotein,
	double Dsolvent) const
{
    assert(node_patches.size() == vertices.size());

    double* fvals = new double[vertices.size()];
    double* hvals = new double[vertices.size()];

    size_t ctr=0;
    for (std::vector<BasicNodePatch>::const_iterator
			it=node_patches.cbegin(), end=node_patches.cend();
		 it != end; ++it)
    {
        const BasicNodePatch& np = *it;
        fvals[ctr] = np.f;
        hvals[ctr] = np.h;
        ++ctr;
    }

    double retval = calculate_energy(kappa, Dprotein, Dsolvent, fvals, hvals);

    delete[] fvals;
    delete[] hvals;
    
    return retval;
}

double Mesh::calculate_energy(
	double kappa,
	double Dprotein,
	double Dsolvent,
	const double fvals[],
	const double hvals[]) const
{

    assert(done_energy_precalcs == true);

    const unsigned int num_patches = node_patches.size();
    const double epsilon_ratio = Dsolvent / Dprotein;
    double energy=0.0;
    double kahan=0.0;

    for (unsigned int ctr=0; ctr < num_patches; ++ctr) {
        const BasicNodePatch& np = node_patches[ctr];
        const double& f = fvals[ctr];
        const double& h = hvals[ctr];
        {  // New scope
            double energy_frag = (h*epsilon_ratio*np.energy_coefficient_h
								  - f*np.energy_coefficient_f);
            double y = energy_frag - kahan;
            double t = energy + y;
            kahan = (t - energy) - y;
            energy = t;
        }  // End scope
    }

    return energy * 0.5 * beep_constants::convert_energy_to_kj_per_mol;
}

std::string Mesh::kinemage_node_patches() const {
    std::ostringstream buf;
    buf << "@vectorlist {patch_edges} off\n";
    for (std::vector<Vertex>::const_iterator
			it=vertices.cbegin(), end=vertices.cend();
		 it != end; ++it)
    {
        buf << it->kinemage_edges(triangles);
        break;
    }
    
    buf << "@trianglelist {node_patches}\n";
    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr) {
        const Vertex& v = vertices[np_ctr];

        for (std::vector<unsigned int>::const_iterator
				it=v.triangle_indices.cbegin(), end=v.triangle_indices.cend();
             it != end; ++it)
        {
            const Triangle& tri = triangles[*it];

            Vector a;
            Vector b;
            try {
                tri.get_anti_clockwise_vertices_from(v, a, b);

                Vector mid_a = (v + a) / 2.0;
                Vector mid_b = (v + b) / 2.0;
                Vector centre = tri.get_planar_centroid();

                // kinemage a quadrilateral with the fh values
                buf << "X " << v.x << " " << v.y << " " << v.z << " "
                    << mid_a.x << " " << mid_a.y << " " << mid_a.z << " "
                    << centre.x << " " << centre.y << " " << centre.z << "\n";
                buf << "X " << v.x << " " << v.y << " " << v.z << " "
                    << centre.x << " " << centre.y << " " << centre.z << " "
                    << mid_b.x << " " << mid_b.y << " " << mid_b.z << "\n";
            }
            catch (Vertex::BadVertex &bv) {
                std::cout << "Vertex " << v << " not a member of Triangle "
			              << tri << ".  Puzzling ... " << std::endl;
            }
        }
        break;
    }
    
    buf << "@vectorlist {qual_pts} color=grey off\n";
    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr) {
        boost::shared_ptr<QuadList> qual_pts
			= node_patches[np_ctr].get_qualocation_points();
        for (QuadList::const_iterator
				it=qual_pts->cbegin(), end=qual_pts->cend();
             it != end; ++it)
        {
            const QuadPoint& qp = *it;
            Vector a = qp.pt();
            Vector b = a + qp.normal() * qp.weight();
            buf << "{}P " << a.x << " " << a.y << " " << a.z << " "
			              << b.x << " " << b.y << " " << b.z << "\n";
        }
        break;
    }
    
    buf << "@spherelist {qual_pts_spheres} color=blue off\n";
    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr) {
        boost::shared_ptr<QuadList> qual_pts
			= node_patches[np_ctr].get_qualocation_points();
        for (QuadList::const_iterator
				it=qual_pts->cbegin(), end=qual_pts->cend();
             it != end; ++it)
        {
            const QuadPoint& qp = *it;
            Vector a = qp.pt();
            buf << "r=0.01 {} " << a.x << " " << a.y << " " << a.z << "\n";
        }
        break;
    }
    
    buf << "@spherelist {quad_pts_spheres} color=blue off\n";
    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr)
    {
        boost::shared_ptr<QuadList> quad_pts
			= node_patches[np_ctr].get_quadrature_points();
        for (QuadList::const_iterator
				it=quad_pts->cbegin(), end=quad_pts->cend();
             it != end; ++it)
        {
            const QuadPoint& qp = *it;
            Vector a = qp.pt();
            buf << "r=0.01 {} " << a.x << " " << a.y << " " << a.z << "\n";
        }
        break;
    }
    
    buf << "@vectorlist {quad_pts} color=grey off\n";
    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr) {
        boost::shared_ptr<QuadList> quad_pts
			= node_patches[np_ctr].get_quadrature_points();
        for (QuadList::const_iterator
				it=quad_pts->cbegin(), end=quad_pts->cend();
             it != end; ++it)
        {
            const QuadPoint& qp = *it;
            Vector a = qp.pt();
            Vector b = a + qp.normal() * fabs(qp.weight());
            buf << "{}P " << a.x << " " << a.y << " " << a.z << " "
			              << b.x << " " << b.y << " " << b.z << "\n";
        }
        break;
    }
    
    return buf.str();
}

std::string Mesh::kinemage_fh_vals(
	double fscale,
	double hscale,
	int num_colours) const
{
    std::ostringstream buf_f;
    std::ostringstream buf_h;

    buf_f << "@trianglelist {fvals}\n";
    buf_h << "@trianglelist {hvals} off\n";

    for (size_t np_ctr=0; np_ctr < node_patches.size(); ++np_ctr) {
        const BasicNodePatch& np = node_patches[np_ctr];
        const Vertex& v = vertices[np_ctr];
        std::string f_name;
        std::string h_name;

        {  // New scope
            std::stringstream s;
            int f_idx = static_cast<int>
					(round(static_cast<double>(num_colours)*np.f / fscale));
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

        {  // New scope
            std::stringstream s;
            int h_idx = static_cast<int>
					(round(static_cast<double>(num_colours)*np.h / hscale));
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

        for (std::vector<unsigned int>::const_iterator
				it=v.triangle_indices.cbegin(), end=v.triangle_indices.cend();
             it != end; ++it)
        {
            const Triangle& tri = triangles[*it];

            Vector a;
            Vector b;
            try {
                tri.get_anti_clockwise_vertices_from(v, a, b);

                Vector mid_a = (v + a) / 2.0;
                Vector mid_b = (v + b) / 2.0;
                Vector centre = tri.get_planar_centroid();

                // kinemage a quadrilateral with the fh values
                buf_f << "X " << v.x << " " << v.y << " " << v.z << " "
                      << mid_a.x << " " << mid_a.y << " " << mid_a.z << " "
                      << f_name << " "
                      << centre.x << " " << centre.y << " " << centre.z << "\n";
                buf_f << "X " << v.x << " " << v.y << " " << v.z << " "
                      << centre.x << " " << centre.y << " " << centre.z << " "
                      << f_name << " "
                      << mid_b.x << " " << mid_b.y << " " << mid_b.z << "\n";

                // kinemage a quadrilateral with the fh values
                buf_h << "X " << v.x << " " << v.y << " " << v.z << " "
                      << mid_a.x << " " << mid_a.y << " " << mid_a.z << " "
                      << h_name << " "
                      << centre.x << " " << centre.y << " " << centre.z << "\n";
                buf_h << "X " << v.x << " " << v.y << " " << v.z << " "
                      << centre.x << " " << centre.y << " " << centre.z << " "
                      << h_name << " "
                      << mid_b.x << " " << mid_b.y << " " << mid_b.z << "\n";
            }
            catch (Vertex::BadVertex &bv) {
                std::cout << "Vertex " << v << " not a member of Triangle "
				          << tri << ".  Puzzling ... " << std::endl;
            }
        }
    }

    buf_f << buf_h.str();

    return buf_f.str();
}

void Mesh::create_node_patches() {
	node_patches.reserve(vertices.size());

	double planar_area = 0;
	total_bezier_area = 0.0;
#ifdef PRETRIPTR
	for (std::vector<BasicTriangle*>::const_iterator
#else
	for (std::vector<std::shared_ptr<Triangle>>::const_iterator
#endif // PRETRIPTR
			it=triangle_ptrs.cbegin(), end=triangle_ptrs.cend();
		 it != end; ++it)
	{
		planar_area += (**it).get_planar_area();
		total_bezier_area += (**it).get_area();
	}

	// instantiates umbrella-shaped Node Patches
	unsigned int ctr=0;
	for (auto& v: vertices)
		node_patches.push_back(BasicNodePatch(ctr++,*this));

	std::cout << "Planar area: " << total_planar_area
	          << " vs. bezier area: " << total_bezier_area << std::endl;
}

Vector Mesh::calculate_charges_centroid() {
	Vector sum(0,0,0);
	if (charges.size() > 0) {
		// iterate over all charges and find centre
		for (std::vector<Charge>::const_iterator
				chit=charges.cbegin(), chend=charges.cend();
			 chit != chend; ++chit)
		{
			sum += chit->position();
		}

		sum /= charges.size();
	}

	return sum;
}

Vector Mesh::calculate_patches_centroid() {
	Vector centroid(0,0,0);
	for (std::vector<BasicNodePatch>::const_iterator
			np_it=node_patches.cbegin(), np_end=node_patches.cend();
		 np_it != np_end; ++np_it)
	{
		const BasicNodePatch& np = *np_it;
		centroid += np.get_centroid() * np.get_bezier_area();            
	}
	
	if (total_bezier_area > 0) {
		centroid /= total_bezier_area;
	}
	else {
		std::cerr << "Error: Surface area of the mesh appears to be zero!"
		          << std::endl;
		throw std::exception();
	}
	
	return centroid;
}


#if 0

void Mesh::precalc_peak_splitting(boost::python::list& fvals,
                                boost::python::list& hvals,
                                double kappa,
                                double Dext,
                                double Dint)
{
    set_vertex_normals();

    assert(len(fvals) == len(hvals));
    assert(len(fvals) == vertices.size());
    boost::shared_array<double> ff(new double[len(fvals)]);
    boost::shared_array<double> hh(new double[len(hvals)]);
    boost::shared_array<double> kahan_ff(new double[len(fvals)]);
    boost::shared_array<double> kahan_hh(new double[len(hvals)]);

    const double epsilon = Dext / Dint;
    const double fconst = 2.0 * epsilon / (Dext*(1.0 + epsilon));
    const double hconst = 2.0 * epsilon / (Dext*(1.0 + epsilon));

    // init the f/h accumulators
    for (unsigned int vertex_idx=0; vertex_idx < vertices.size(); ++vertex_idx)
    {
        double& fpart = ff[vertex_idx];
        double& hpart = hh[vertex_idx];
        double& kahan_f = kahan_ff[vertex_idx];
        double& kahan_h = kahan_hh[vertex_idx];
        kahan_f = 0.0;
        kahan_h = 0.0;

        const Vertex& v = vertices[vertex_idx];
        const Vector& n0 = v.normal;

        double pot=0.0;
        Vector field(0,0,0);
        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end();
                ch_it != ch_end;
                ++ch_it)
        {
            double pot_term = ch_it->charge * upt(v, ch_it->position, kappa);
            pot += pot_term;

            Vector r = ch_it->position - v;
            double r2 = r.length2();
            r.normalise();
            double field_term = ch_it->charge * ONE_OVER_4PI / r2;
            field.x += field_term * r.x;
            field.y += field_term * r.y;
            field.z += field_term * r.z;
        }
        fpart = pot * fconst * ::fGeometricCorrection(0.5, epsilon);
        hpart = field.dot(n0) * hconst * ::hGeometricCorrection(0.5, epsilon);

    }

    double magic_area = get_magic_area();
    int tri_progress_ctr=0;
    for (std::vector<Triangle>::const_iterator tri_it=triangles.begin(),
            tri_end = triangles.end();
            tri_it != tri_end;
            ++tri_it)
    {
        const Triangle& tri = *tri_it;

        //int max_cached_n = START_N;
        int max_cached_n = static_cast<int>(ceil(sqrt(tri.area / magic_area)));

        boost::shared_array<QuadPoint> quad_points = precalc_gl_points_with_vals(tri, max_cached_n, charges, kappa, Dext, Dint);
        for (unsigned int vertex_idx=0; vertex_idx < vertices.size(); ++vertex_idx)
        {
            const Vertex& v = vertices[vertex_idx];
            double& fpart = ff[vertex_idx];
            double& hpart = hh[vertex_idx];
            double& kahan_f = kahan_ff[vertex_idx];
            double& kahan_h = kahan_hh[vertex_idx];

            {
                double fnum = calculate_Bpt_Apt(v, tri, quad_points, max_cached_n, charges, kappa, Dext, Dint);
                double y = fnum - kahan_f;
                double t = fpart + y;
                kahan_f = (t - fpart) - y;
                fpart = t;
            }
            {
                double hnum = calculate_Dpt_Cpt(v, v.normal, tri, quad_points, max_cached_n, charges, kappa, Dext, Dint);
                double y = hnum - kahan_h;
                double t = hpart + y;
                kahan_h = (t - hpart) - y;
                hpart = t;
            }
        }
        std::cout << "Integrated over  " << ++tri_progress_ctr << " surface triangles out of " << triangles.size() << " (max GL recursions: " << max_cached_n << ")" << std::endl;
    }

    for (unsigned int ctr=0; ctr < vertices.size(); ++ctr)
    {
        fvals[ctr] = ff[ctr];
        hvals[ctr] = hh[ctr];
    }

    return;
}

void Mesh::precalc_rhs(boost::python::list& fvals,
                    boost::python::list& hvals,
                    double Dext)
{
    set_vertex_normals();

    assert(len(fvals) == len(hvals));
    assert(len(fvals) == vertices.size());

    // init the f/h accumulators
    for (unsigned int vertex_idx=0; vertex_idx < vertices.size(); ++vertex_idx)
    {
        double frhs = 0.0;
        double hrhs = 0.0;

        const Vertex& v = vertices[vertex_idx];
        const Vector& n0 = v.normal;

        double kahan_gpk=0.0, kahan_dgpk_dn=0.0;
        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end();
            ch_it != ch_end;
            ++ch_it)
        {
            {
                double num = ch_it->charge * Gpt(v, ch_it->position);
                double y = num - kahan_gpk;
                double t = frhs + y;
                kahan_gpk = (t - frhs) - y;
                frhs = t;
            }
            {
                double num = ch_it->charge * dGpt_dn(ch_it->position,v,n0);
                double y = num - kahan_dgpk_dn;
                double t = hrhs + y;
                kahan_dgpk_dn = (t - hrhs) - y;
                hrhs = t;
            }
        }

        fvals[vertex_idx] = frhs / Dext;
        hvals[vertex_idx] = hrhs / Dext;
    }

    return;
}

Vector Mesh::calc_reaction_field_force(const Offsets& offs,
                                    double kappa,
                                    double Dprotein,
                                    double Dsolvent,
                                    const SpharmHolder& sph) const
{
    Vector force(0,0,0);

    unsigned int first_np_idx = offs.vertex_numbering_offset;
    unsigned int num_patches = offs.num_vertices;
    unsigned int first_ch_idx = offs.charge_numbering_offset;
    unsigned int num_charges = offs.num_charges;

    double epsilon_ratio = Dsolvent / Dprotein;

    // loop over all charge positions
    for (unsigned int ch_idx=first_ch_idx; ch_idx < first_ch_idx+num_charges; ++ch_idx)
    {
        const Vector& pos = charges[ch_idx].position;
        const double& charge = charges[ch_idx].charge;

        Vector fpt(0,0,0);
        for (unsigned int np_idx=first_np_idx; np_idx < first_np_idx+num_patches; ++np_idx)
        {
            const BasicNodePatch& np = node_patches[np_idx];
            Vector r = np.node - pos;
            //double hval = sph.evaluate(r, offs.origin);
            double hval = np.h;
            double r3 = r.length2();
            r3 *= sqrt(r3);
            fpt = fpt + (r * hval * np.area / r3);
//             for (std::vector<Triangle>::const_iterator tri_it=np.triangles.begin(), tri_end=np.triangles.end();
//                  tri_it != tri_end;
//                  ++tri_it)
//             {
//                 Vector r = pos - tri_it->centre;
//                 double hval = sph.evaluate(r, offs.origin);
//                 double r3 = r.length2();
//                 r3 *= sqrt(r3);
//                 fpt = fpt + (r * hval * tri_it->area / r3);
//             }
        }
        force = force + fpt*charge*(1.0 - epsilon_ratio);
    }

    return force*ONE_OVER_4PI;
}
#endif // if 0

// Add a mesh to the list
boost::shared_ptr<ListedMesh>
MeshList::add(const std::string& filename, bool force_planar)
{
	// add mesh to library
	boost::shared_ptr<ListedMesh>
		mesh_ptr(new ListedMesh(size(), filename, force_planar));
	push_back( mesh_ptr );
	return mesh_ptr;
}

void MeshList::reset_energy_precalcs() {
    for (MeshList::iterator mit=begin(), mend=end(); mit != mend; ++mit) {
		(**mit).init_energy_precalcs();
    }
}


