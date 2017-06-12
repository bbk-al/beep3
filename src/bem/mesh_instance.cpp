/*
* mesh_instance.cpp
*
*  Created on: 27 Jul 2010
*      Author: david
*/

#include "mesh_instance.h"
#include <iostream>
#include <iomanip>



MeshInstance::MeshInstance(unsigned int mesh_lib_id,
                            unsigned int mesh_instance_id,
                            const MeshList& mesh_library,
                            const Vector& offset,
                            const Quaternion& rot,
                            double protein_dielectric,
                            double solvent_dielectric,
                            unsigned int num_quad_points,
                            unsigned int num_qual_points,
                            bool _silent) :
            lib_id(mesh_lib_id),
            instance_id(mesh_instance_id),
            silent(_silent),
            xyz_offset(offset),
            rotation(rot),
            quad_points_per_triangle(num_quad_points),
            qual_points_per_triangle(num_qual_points)
            
{
    // rather essential- copy the shared_ptr to the underlying mesh type!
    assert(lib_id <= mesh_library.size());
    mesh_ptr = mesh_library[lib_id];
    
    init();
    set_dielectrics(protein_dielectric, solvent_dielectric);
}

void MeshInstance::init()
{
    const Mesh& mesh = *mesh_ptr;

    // set BasicNodePatches
    patches.clear();
    const std::vector<BasicNodePatch>& library_node_patches = mesh.get_node_patches();
    for (std::vector<BasicNodePatch>::const_iterator it=library_node_patches.begin(), end=library_node_patches.end(); it != end; ++it)
    {
        // copy the node patches from the reference mesh
        // hahaha the wonders of polymorphism- create the NodePatch but immediately store it as a pointer to the base class
        // (i really really hope the virtual destructor works correctly...)
        boost::shared_ptr<BasicNodePatch> ptr(new NodePatch(*it, *this)); 
        patches.push_back(ptr);
    }

    // set Charges
    charges.clear();
    const std::vector<Charge>& library_charges = mesh.get_charges();
    for (std::vector<Charge>::const_iterator it=library_charges.begin(), end=library_charges.end(); it != end; ++it)
    {
        boost::shared_ptr<Charge> ptr(new Charge(*it, mesh.centre, rotation, xyz_offset));
        charges.push_back(ptr);
    }
    
    // set radius
    radius = mesh.get_radius();

    // create a simple octree of the mesh instance (for collision / grid checks)
    mesh_tree.reset(new Octree<Node<BasicNodePatch>,BasicNodePatch>(10, xyz_offset, radius*4.));
    for (PatchList::iterator it=patches.begin(), end=patches.end(); it != end; ++it)
    {
        mesh_tree->add(*it);
    }
}

bool MeshInstance::pt_is_internal(const Vector& pt) const
{
    // if the point is outside of the maximum radius then cannot be
    // inside the mesh instance
    if ((pt - xyz_offset).length() >= radius) { return false; }

    // ok so it's within the maximum radius, still not necessarily
    // internal -- find the nearest node patch and compare this 
    // point to the normal vector of the patch
    const BasicNodePatch& nearest_np = mesh_tree->get_nearest(pt);
    if ((pt - nearest_np).dot(nearest_np.get_normal()) > 0)
    {
        return false;
    }

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

    if (charge_fmm_tree.get() == NULL)
    {
        charge_fmm_tree.reset(new fmm::FMM_Octree_6FIG_ACCURACY(10000, xyz_offset, radius * 4.));
        for (std::vector<boost::shared_ptr<Charge> >::iterator charge_it=charges.begin(), charge_end=charges.end();
             charge_it != charge_end; 
             ++charge_it)
        {
            charge_fmm_tree->add(*charge_it);
        }
        charge_fmm_tree->solve(0.0); // solve for kappa=0, i.e. non-screened Coulomb potential
    }
    pot = charge_fmm_tree->calculate_potential(pt);
    pot *= ONE_OVER_4PI / Dprotein;

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator np_it=patches.begin(), np_end=patches.end(); 
         np_it != np_end; 
         ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double Gpt_val = Gpt(pt, np);
        double dGpt_val = dGpt_dn(pt, np, n);
        accum += np.get_bezier_area() * (-Gpt_val*np.h*np.dielectric_ratio + dGpt_val*np.f);
    }

    return pot+accum;
}

double MeshInstance::get_potential_at_external_pt(const Vector& pt, double kappa) const
{
    assert(this->pt_is_internal(pt) == false);

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator np_it=patches.begin(), np_end=patches.end(); 
         np_it != np_end; 
         ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double upt_val = upt(pt, np, kappa);
        double dupt_val = dupt_dn(pt, np, n, kappa);
        accum += np.get_bezier_area() * (dupt_val*np.f - upt_val*np.h);
       
    }
    
    return accum;
}
double MeshInstance::calculate_energy(double kappa,
                                      double fvals[],
                                      double hvals[]) const

{
    // Dsolvent and Dprotein should have been set by the set_dielectrics function call
    double E = mesh_ptr->calculate_energy(kappa, Dprotein, Dsolvent, fvals, hvals);
    std::cout << "Energy for mesh " << instance_id << " (lib_id=" << lib_id << ") = " << std::setprecision(10) << E << std::endl;
    return E;
}

Vector MeshInstance::calculate_force(double kappa,
                                     double fvals[],
                                     double hvals[]) const

{
    // Dsolvent and Dprotein should have been set by the set_dielectrics function call
    KahanVector MST_ext,MST_int,dbf,ionic;
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent, fvals, hvals, MST_ext, MST_int, dbf, ionic);
    Vector force(*MST_ext);
    force.apply_rotation(rotation);
    return force;
}

void MeshInstance::calculate_forces(double kappa,
                                    double fvals[],
                                    double hvals[],
                                    KahanVector& qE,
                                    KahanVector& MST_ext,
                                    KahanVector& MST_int,
                                    KahanVector& dbf,
                                    KahanVector& ionic) const

{
    // Dsolvent and Dprotein should have been set by the set_dielectrics function call
    qE = mesh_ptr->calculate_qe_force(Dprotein, Dsolvent, fvals, hvals);
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent, fvals, hvals, MST_ext, MST_int, dbf, ionic);
    
    (qE).apply_rotation(rotation);
    (MST_int).apply_rotation(rotation);
    (MST_ext).apply_rotation(rotation);
    (ionic).apply_rotation(rotation);
    (dbf).apply_rotation(rotation);
    
    return;
}

void MeshInstance::kinemage_fh_vals(double fscale, double hscale, int num_colours, std::ostringstream& buf_f, std::ostringstream& buf_h) const
{
    
    for (size_t np_ctr=0; np_ctr < patches.size(); ++np_ctr)
    {
        const NodePatch& np = dynamic_cast<NodePatch&>(*(patches[np_ctr]));
        
        std::string f_name;
        std::string h_name;

        // figure out the colour name corresponding to the fval
        {
            std::stringstream s;
            int f_idx = static_cast<int>(round(static_cast<double>(num_colours)*np.f / fscale));
            if (f_idx < 1){
                f_idx = abs(f_idx);
                int fcol = f_idx <= num_colours ? f_idx : num_colours+1;
                s << "red" << "_" << fcol;
            } else {
                int fcol = f_idx <= num_colours ? f_idx : num_colours+1;
                s << "blue" << "_" << fcol;
            }
            f_name = s.str();
        }

        // figure out the colour name corresponding to the fval
        {
            std::stringstream s;
            int h_idx = static_cast<int>(round(static_cast<double>(num_colours)*np.h / hscale));
            if (h_idx < 1){
                h_idx = abs(h_idx);
                int hcol = h_idx <= num_colours ? h_idx : num_colours+1;
                s << "red" << "_" << hcol;
            } else {
                int hcol = h_idx <= num_colours ? h_idx : num_colours+1;
                s << "blue" << "_" << hcol;
            }
            h_name = s.str();
        }

        Vector local_centre = mesh_ptr->get_centre();
        std::vector<PointNormal> edge_points;
        np.get_edge_points(edge_points);
        for (std::vector<PointNormal>::const_iterator it=edge_points.begin(), next_it, end=edge_points.end(); it!=end; ++it)
        {
            next_it = it+1;
            if (next_it == end) { next_it = edge_points.begin(); }
            
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

    return;
}
