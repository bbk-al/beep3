/*      Author: david    Created on: 27 Jul 2010*/
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file mesh_instance.h
 * \brief This module declares the mesh instance and list classes.
 *
 * The MeshInstance class adds positions and orientations to Mesh objects,
 * which are stored as references in MeshInstances.  Many MeshInstance
 * objects may refer to the same Mesh object.
 *
 *	This module contains the following public classes:
 *	- MeshInstance -- the mesh instance class
 *	- MeshInstanceList -- a list of MeshInstance
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

	//! Primary constructor
    MeshInstance(unsigned int mesh_lib_id,
                 unsigned int mesh_instance_id,
                 MeshList& mesh_library,
                 const Vector& offset,
                 const Quaternion& rot,
                 double protein_dielectric,
                 double solvent_dielectric,
                 unsigned int num_quad_points,
                 unsigned int num_qual_points,
                 bool _silent=false);
                
	inline MeshInstance(const MeshInstance& other); //! Copy constructor
	inline MeshInstance(MeshInstance&& other);		//! Move constructor
#if 0
//NB
	MeshInstance& operator=(const MeshInstance&) = default;	//! Copy assignment
	MeshInstance& operator=(MeshInstance&&) = default;		//! Move assignment
#endif

	// Destructor
    virtual ~MeshInstance() = default;

	// get_ methods
	unsigned int get_id() const { return instance_id; }
    const Vector& get_xyz_offset() const { return xyz_offset; }
    const Quaternion& get_rotation() const { return rotation; }
    double get_radius() const { return radius; }
    const Vector& get_local_centre() const { return mesh_ptr->get_centre(); }
    const Mesh& get_ref_mesh() const { return *mesh_ptr; }

    std::vector< boost::shared_ptr<BasicNodePatch> >&
    get_node_patches()
	{
		return patches;
	}
    const std::vector< boost::shared_ptr<BasicNodePatch> >&
	get_node_patches() const
	{
		//NB For more complex cases:
		// return const_cast<RET &>(
		//			static_cast<const CLASS &>(*this).get_xxx() );
		return patches;
	}

    std::vector< boost::shared_ptr<Charge> >& get_charges() { return charges; }
    const std::vector< boost::shared_ptr<Charge> >&
    get_charges() const
	{
		return charges;
	}

    unsigned int get_quad_points_per_triangle() const {
		return quad_points_per_triangle;
	}
    unsigned int get_qual_points_per_triangle() const {
		return qual_points_per_triangle;
	}
    
	// get_ index methods
	//NB choice of const/non-const is set in pybeep: could add variations here
    BasicNodePatch& get_node_patch(size_t idx) const { return *(patches[idx]); }
    const Charge& get_charge(size_t idx) const { return *(charges[idx]); }

	// get counts
    inline unsigned int get_num_node_patches() const { return patches.size(); }
    inline unsigned int get_num_charges() const { return charges.size(); }

	// booleans
    bool isSilent() const { return silent; }
    bool pt_is_internal(const Vector& pt) const;

	// get - other
    inline Vector get_h_squared() const;

    inline std::vector<Charge>
	get_ecm_charges(const std::string& ecm_filename) const
	{
        return mesh_ptr->get_ecm_charges(ecm_filename, rotation, xyz_offset);
    }
    double get_potential_at_internal_pt(const Vector& pt) const;
    double get_potential_at_external_pt(const Vector& pt, double kappa) const;
                           
	// set methods
    void set_silent() { silent = true; }
    void unset_silent() { silent = false; }

    inline void
	set_dielectrics(double protein_dielectric, double solvent_dielectric);

    inline void set_quad_points_per_triangle(unsigned int quad_points);
    inline void set_qual_points_per_triangle(unsigned int qual_points);
	void set_unique_patch_id(unsigned int patch_ctr);

	// Move the instance
	int move(const Vector& translate, const Quaternion& rotate);
     
	// Calculations
    double calculate_energy(
		double kappa,
		double fvals[],
		double hvals[]) const;

    Vector calculate_force(
		double kappa,
		double fvals[],
		double hvals[]) const;

    void calculate_forces(
		double kappa,
		double fvals[],
		double hvals[],
		KahanVector& qE,
		KahanVector& MST_ext,
		KahanVector& MST_int,
		KahanVector& dbf,
		KahanVector& ionic) const;
                          
    void kinemage_fh_vals(
		double fscale,
		double hscale,
		int num_colours,
		std::ostringstream& fbuf,
		std::ostringstream& hbuf) const;

//NB is const really valid?
	bool init_fh_vals(double *x, unsigned int& xctr,
                      unsigned int offset, bool preconditioned) const;
	unsigned int reset_fh_vals(const boost::shared_array<double>& f_lhs, 
	                           const boost::shared_array<double>& h_lhs,
							   unsigned int start_ctr);
    
private:

    void init();

    inline void set_quad_points();
    inline void set_qual_points();
	void reset_mesh_tree();

	// Attributes
    bool silent;

    double Dprotein; // The dielectric constant of this protein  //NB ??
    double Dsolvent; // The dielectric constant of the solvent  //NB ??
    Vector xyz_offset;		// Centre of mesh instance
    Quaternion rotation;	// Net rotation relative to reference mesh
    double radius;

    const unsigned int instance_id;
#ifdef __DELETED__
public:
    boost::shared_ptr<Mesh> mesh_ptr; // shared_ptr to the reference Mesh
private:
#else
    boost::shared_ptr<ListedMesh> mesh_ptr; // shared_ptr to the reference Mesh
#endif
    
#ifdef __DELETED__
public:
    std::vector< boost::shared_ptr<BasicNodePatch> > patches;
private:
#endif
    mutable std::vector< boost::shared_ptr<Charge> > charges;

    boost::scoped_ptr< Octree< Node<BasicNodePatch>, BasicNodePatch > >
		mesh_tree;
    mutable boost::scoped_ptr< fmm::FMM_Octree_6FIG_ACCURACY > charge_fmm_tree;
    
    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;
};

#ifdef __DELETED__
typedef std::vector< boost::shared_ptr<MeshInstance> > MeshInstanceList;
#else
class MeshInstanceList : public std::vector< boost::shared_ptr<MeshInstance> >
{
public:

	MeshInstanceList() {}

	// Library methods
	boost::shared_ptr<ListedMesh>
		addMesh(const std::string& filename, bool force_planar)
	{
		return library.add(filename, force_planar);
	}
	size_type librarySize() { return library.size(); }
#ifndef __DELETED__
	Mesh& getMesh(unsigned int id) { return *(library[id]); }
	void reset_energy_precalcs() { library.reset_energy_precalcs(); }
#endif // ! __DELETED__

	// Instance methods
	boost::shared_ptr<MeshInstance>
		add(unsigned int mesh_lib_id, unsigned int mesh_instance_id,
            const Vector& offset, const Quaternion& rotation,
            double Dprotein, double Dsolvent,
            unsigned int num_quad_points, unsigned int num_qual_points,
            bool _silent=false);
	boost::shared_ptr<MeshInstance>
	move(unsigned int mesh_instance_id,
		 const Vector& offset, const Quaternion& rotation,
		 double Dprotein, double Dsolvent, bool _silent=false);

//NB const?
	bool init_library_fh_vals(double *x, unsigned int offset);
	void reset_library_fh_vals(const boost::shared_array<double>& f_lhs,
					           const boost::shared_array<double>& h_lhs);

	size_t get_total_patches() const;

private:
	MeshList library;
};
#endif

// inline methods
inline MeshInstance::MeshInstance(const MeshInstance& other) : 
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
//NB emplace?
	patches.insert(patches.end(), other.patches.begin(), other.patches.end());
	charges.insert(charges.end(), other.charges.begin(), other.charges.end());
}
    
//NB should use std::move here? And emplace?
inline MeshInstance::MeshInstance(MeshInstance&& other) : 
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
//NB emplace?
	patches.insert(patches.end(), other.patches.begin(), other.patches.end());
	charges.insert(charges.end(), other.charges.begin(), other.charges.end());
}
    
inline Vector MeshInstance::get_h_squared() const {
	KahanVector h_squared;
	for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
			it=patches.cbegin(), end=patches.cend();
		 it != end; ++it)
	{
		const BasicNodePatch& np = **it;
		h_squared += np.get_normal() * np.h * np.h;
	}
	return *h_squared;
}

inline void MeshInstance::set_dielectrics(
	double protein_dielectric,
	double solvent_dielectric)
{
	Dprotein = protein_dielectric; // need this later for calculating energies
	Dsolvent = solvent_dielectric;
	double epsilon = solvent_dielectric / protein_dielectric;
	for (PatchList::iterator nit=patches.begin(), nend=patches.end();
		 nit != nend; ++nit)
	{
		(**nit).set_dielectric_ratio(epsilon);
	}
}

inline void
MeshInstance::set_quad_points_per_triangle(unsigned int quad_points)
{ 
	quad_points_per_triangle = quad_points; 
	set_quad_points();
}

inline void
MeshInstance::set_qual_points_per_triangle(unsigned int qual_points)
{
	qual_points_per_triangle = qual_points; 
	set_qual_points();
}

inline void MeshInstance::set_quad_points() {
	// set the number of qualocation / quadrature points on each node patch
	for (std::vector< boost::shared_ptr<BasicNodePatch> >::iterator
			it=patches.begin(), end=patches.end();
		 it != end; ++it)
	{
		(*it)->set_quad_points_per_triangle(quad_points_per_triangle);
	}
}

inline void MeshInstance::set_qual_points() {
	for (std::vector< boost::shared_ptr<BasicNodePatch> >::iterator
			it=patches.begin(), end=patches.end();
		 it != end; ++it)
	{
		(*it)->set_qual_points_per_triangle(qual_points_per_triangle);
	}
}

#endif /* MESH_INSTANCE_H_ */

