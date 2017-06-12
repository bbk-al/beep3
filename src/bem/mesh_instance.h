/*
* mesh_instance.h
*
*  Created on: 27 Jul 2010
*      Author: david
*/

#ifndef MESH_INSTANCE_H_
#define MESH_INSTANCE_H_

#include "mesh.h"
#include "node_patch.h"
#include "../common/math_vector.h"
#include <boost/shared_ptr.hpp>
#include "../fmm/octree.h"
#include "../fmm/fmm_octree.h"

class MeshInstance {

public:

    friend class BEEP;
    
    MeshInstance(unsigned int mesh_lib_id,
                 unsigned int mesh_instance_id,
                 const MeshList& mesh_library,
                 const Vector& offset,
                 const Quaternion& rot,
                 double protein_dielectric,
                 double solvent_dielectric,
                 unsigned int num_quad_points,
                 unsigned int num_qual_points,
                 bool _silent=false);
                
    MeshInstance(const MeshInstance& other) : 
        lib_id(other.lib_id), 
        instance_id(other.instance_id), 
        silent(other.silent), 
        Dprotein(other.Dprotein), 
        Dsolvent(other.Dsolvent), 
        xyz_offset(other.xyz_offset), 
        rotation(other.rotation), 
        radius(other.radius),
        mesh_ptr(other.mesh_ptr),
        quad_points_per_triangle(other.quad_points_per_triangle),
        qual_points_per_triangle(other.qual_points_per_triangle)
    {
        patches.insert(patches.end(), other.patches.begin(), other.patches.end());
        charges.insert(charges.end(), other.charges.begin(), other.charges.end());
    }
    
    ~MeshInstance() {}

    void init();
    inline unsigned int get_num_node_patches() const { return patches.size(); }
    inline unsigned int get_num_charges() const { return charges.size(); }
     
    // accessors for patches/charges
    BasicNodePatch& get_node_patch(size_t idx) const { return *(patches[idx]); }
    const Charge& get_charge(size_t idx) const { return *(charges[idx]); }

    double calculate_energy(double kappa,
                            double fvals[],
                            double hvals[]) const;
    Vector calculate_force(double kappa,
                           double fvals[],
                           double hvals[]) const;
    void calculate_forces(double kappa,
                          double fvals[],
                          double hvals[],
                          KahanVector& qE,
                          KahanVector& MST_ext,
                          KahanVector& MST_int,
                          KahanVector& dbf,
                          KahanVector& ionic) const;
                          
    Vector get_h_squared() const
    {
        KahanVector h_squared;
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator it=patches.begin(), end=patches.end(); it != end; ++it)
        {
            const BasicNodePatch& np = **it;
            h_squared += np.get_normal() * np.h * np.h;
        }
        return *h_squared;
    }


    inline std::vector<Charge> get_ecm_charges(const std::string& ecm_filename) const {
        return mesh_ptr->get_ecm_charges(ecm_filename, rotation, xyz_offset);
    }
                           
    inline void set_dielectrics(double protein_dielectric, double solvent_dielectric)
    {
        Dprotein = protein_dielectric; // need this later for calculating energies
        Dsolvent = solvent_dielectric;
        double epsilon = solvent_dielectric / protein_dielectric;
        for (PatchList::iterator nit=patches.begin(), nend=patches.end(); nit != nend; ++nit)
        {
            (**nit).set_dielectric_ratio(epsilon);
        }

    }
    
    inline const Mesh& get_ref_mesh() const {
        return *mesh_ptr; 
    }
    
    inline bool isSilent() const { return silent; }
    inline void set_silent() { silent = true; }
    inline void unset_silent() { silent = false; }
    
    inline const Vector& get_xyz_offset() const { return xyz_offset; }
    inline const Quaternion& get_rotation() const { return rotation; }
    inline double get_radius() const { return radius; }
    inline const Vector& get_local_centre() const { return mesh_ptr->get_centre(); }
    
    void kinemage_fh_vals(double fscale, double hscale, int num_colours, std::ostringstream& fbuf, std::ostringstream& hbuf) const;
    bool pt_is_internal(const Vector& pt) const;
    double get_potential_at_internal_pt(const Vector& pt) const;
    double get_potential_at_external_pt(const Vector& pt, double kappa) const;

    inline unsigned int get_quad_points_per_triangle() const { return quad_points_per_triangle; }
    inline unsigned int get_qual_points_per_triangle() const { return qual_points_per_triangle; }

    inline void set_quad_points_per_triangle(unsigned int quad_points) { 
        quad_points_per_triangle = quad_points; 
        set_quad_points();
    }
    
    inline void set_qual_points_per_triangle(unsigned int qual_points) { 
        qual_points_per_triangle = qual_points; 
        set_qual_points();
    }
    
private:

    inline void set_quad_points() 
    {
        // set the number of qualocation / quadrature points on each node patch
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::iterator it=patches.begin(), end=patches.end(); it != end; ++it)
        {
            (*it)->set_quad_points_per_triangle(quad_points_per_triangle);
        }
    }
    
    inline void set_qual_points()
    {
        for (std::vector< boost::shared_ptr<BasicNodePatch> >::iterator it=patches.begin(), end=patches.end(); it != end; ++it)
        {
            (*it)->set_qual_points_per_triangle(qual_points_per_triangle);
        }
    }
    
    const unsigned int lib_id;
    const unsigned int instance_id;
    bool silent;

    double Dprotein; // The dielectric constant of this protein
    double Dsolvent; // The dielectric constant of this protein
    Vector xyz_offset;
    Quaternion rotation;
    double radius;

    // shared_ptr to the underlying reference Mesh
    boost::shared_ptr<Mesh> mesh_ptr;
    
    std::vector< boost::shared_ptr<BasicNodePatch> > patches;
    mutable std::vector< boost::shared_ptr<Charge> > charges;

    boost::scoped_ptr< Octree< Node<BasicNodePatch>, BasicNodePatch > > mesh_tree;
    mutable boost::scoped_ptr< fmm::FMM_Octree_6FIG_ACCURACY > charge_fmm_tree;
    
    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;

};

typedef std::vector< boost::shared_ptr<MeshInstance> > MeshInstanceList;

#endif /* MESH_INSTANCE_H_ */

