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
#include <unordered_set>

class MeshInstanceList;

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

#ifdef PREHYDROPHOBIC
    std::vector< boost::shared_ptr<Charge> >& get_charges() { return charges; }
    const std::vector< boost::shared_ptr<Charge> >& get_charges() const {
		return charges;
	}
#else
    const std::vector<boost::shared_ptr<Charge>>&
		get_charges(bool all = false) const
	{
		return all ? allCharges : charges;
	}
    std::vector<boost::shared_ptr<Charge>>& get_charges(bool all = false) {
		return all ? allCharges : charges;
	}
#endif // PREHYDROPHOBIC

    unsigned int get_quad_points_per_triangle() const {
		return quad_points_per_triangle;
	}
    unsigned int get_qual_points_per_triangle() const {
		return qual_points_per_triangle;
	}
    
	// get_ index methods
	//NB choice of const/non-const is set in pybeep: could add variations here
    BasicNodePatch& get_node_patch(size_t idx) const { return *(patches[idx]); }
#ifdef PREHYDROPHOBIC
	const Charge& get_charge(size_t idx) const { return *(charges[idx]); }
#else
    const Charge& get_charge(size_t idx, bool all = false) const {
		return all ? *(allCharges[idx]) : *(charges[idx]);
	}
#endif // PREHYDROPHOBIC

	// get counts
    unsigned int get_num_node_patches() const { return patches.size(); }
#ifdef PREHYDROPHOBIC
	unsigned int get_num_charges() const { return charges.size(); }
#else
    unsigned int get_num_charges(bool all = false) const {
		return all ? allCharges.size() : charges.size();
	}
#endif // PREHYDROPHOBIC

	// Local types
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
	using pt_set = std::unordered_set<Vector, hashFunc, equalsFunc>;
	// Default map of hit vertices - must be reset by caller of pt_triangle
	static pt_set pt_hitVertices;

	// booleans
    bool isSilent() const { return silent; }
    bool pt_is_internal(const Vector& pt,
					pt_set& hitVertices = MeshInstance::pt_hitVertices) const;

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
#define MI_MOVE
#ifndef MI_MOVE
	[[deprecated("Method not working - use MeshInstanceList::move()")]]
#endif
	MeshInstance& move(const Vector& translate, const Quaternion& rotate);

#ifndef PREHYDROPHOBIC
	// Non-electrostatic energy
	void update_energy(MeshInstanceList& mis);
	void revert_energy(MeshInstanceList& mis);
#endif // PREHYDROPHOBIC
     
	// Calculations
    double calculate_energy(
#ifndef PREHYDROPHOBIC
		bool electrostatics,	// if false, only do non-electrostatic energy
#endif // PREHYDROPHOBIC
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
                          
#ifndef PREHYDROPHOBIC
	// For backward compatibility only
	[[deprecated("Use kinemage_vals() instead with defaulted scales")]]
#endif
	void kinemage_fh_vals(
		double fscale,
		double hscale,
		int num_colours,
		std::ostringstream& buf_f,
		std::ostringstream& buf_h) const;

	void kinemage_vals(
		int num_colours,
#ifndef PREHYDROPHOBIC
		std::ostringstream& buf_hy,
		std::ostringstream& buf_s,
		std::ostringstream& buf_e,
		std::ostringstream& buf_he,
		std::ostringstream& buf_lj,
#endif // PREHYDROPHOBIC
		std::ostringstream& buf_f,
		std::ostringstream& buf_h) const;

	// Allow caller to specify scales, to support backward compatibility...
	template<unsigned int S>
	void kinemage_vals(
		int num_colours,
#ifndef PREHYDROPHOBIC
		std::ostringstream& buf_hy,
		std::ostringstream& buf_s,
		std::ostringstream& buf_e,
		std::ostringstream& buf_he,
		std::ostringstream& buf_lj,
#endif // PREHYDROPHOBIC
		std::ostringstream& buf_f,
		std::ostringstream& buf_h,
		double (&scale)[S]) const;	// S must match buf arg count

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
#if 1 // TESTING BOTH pt_is_internal METHODS else delete and re-enable DELETED
//#ifdef DELETED
	void reset_mesh_tree();
#endif // DELETED
#ifndef PREHYDROPHOBIC
	unsigned int resize_pot(unsigned int n);
	void calculate_boundary_energy(	//! Non-electrostatic energies
		BasicNodePatch& np,			//!\param The patch to calculate for
		const MeshInstance& omi,	//!\param The impinging MeshInstance
		std::vector<unsigned int>& track,	//!\param omi.patches.size(), all 0
		pt_set& hitVertices = pt_hitVertices); //!\param Detect vertex repeats
#endif // PREHYDROPHOBIC

	// Attributes
    bool silent;

    double Dprotein; // The dielectric constant of this protein  //NB ??
    double Dsolvent; // The dielectric constant of the solvent  //NB ??
    Vector xyz_offset;		// Centre of mesh instance
    Quaternion rotation;	// Net rotation relative to reference mesh
    double radius;

    const unsigned int instance_id;
    boost::shared_ptr<ListedMesh> mesh_ptr; // shared_ptr to the reference Mesh
    
    std::vector< boost::shared_ptr<BasicNodePatch> > patches;
    mutable std::vector< boost::shared_ptr<Charge> > charges;
#ifndef PREHYDROPHOBIC
    mutable std::vector< boost::shared_ptr<Charge> > allCharges;
#endif // PREHYDROPHOBIC

#ifdef DELETED // redundant - pt_is_internal uses a different method
    boost::scoped_ptr< Octree< Node<BasicNodePatch>, BasicNodePatch > >
		mesh_tree;
#endif // DELETED
    mutable boost::scoped_ptr< fmm::FMM_Octree_6FIG_ACCURACY > charge_fmm_tree;
    
    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;

#ifndef PREHYDROPHOBIC
	// Non-initialised  i.e. reset when used
	std::vector<double> hepot[2];	// hydrophobic potential (current or old)
	std::vector<double> ljpot[2];	// LJ potential (current or old)
	unsigned int psel;			// To select which pots to use (current or old)
	double hea;					// Hydrophobic area in range
#endif // PREHYDROPHOBIC

#ifdef USING_MAXD  // maxd no longer required -- see mesh_instance.cpp-old
	// maximum dimension of a triangle for x and y directions
	std::unique_ptr<double[]> maxd;
#endif // maxd
};

#ifdef DELETED
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
	Mesh& getMesh(unsigned int id) { return *(library[id]); }
	void reset_energy_precalcs() { library.reset_energy_precalcs(); }

	// Instance methods
	boost::shared_ptr<MeshInstance>
		add(unsigned int mesh_lib_id, unsigned int mesh_instance_id,
            const Vector& offset, const Quaternion& rotation,
            double Dprotein, double Dsolvent,
            unsigned int num_quad_points, unsigned int num_qual_points,
            bool _silent=false);
#ifdef MI_MOVE
	[[deprecated("No support for hydrophobicity - use MeshInstance::move()")]]
#endif
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
#ifndef PREHYDROPHOBIC
		psel(other.psel),
#endif // PREHYDROPHOBIC
#ifdef USING_MAXD  // maxd no longer required
		maxd(std::make_unique<double[]>(2)), // Each copy must have own version
#endif // maxd
        quad_points_per_triangle(other.quad_points_per_triangle),
        qual_points_per_triangle(other.qual_points_per_triangle)
{
	patches.insert(patches.end(), other.patches.begin(), other.patches.end());
	charges.insert(charges.end(), other.charges.begin(), other.charges.end());
#ifndef PREHYDROPHOBIC
	allCharges.insert(allCharges.end(),
					  other.allCharges.begin(), other.allCharges.end());
	for (int s = 0; s < 2; s++) {
		hepot[s].insert(hepot[s].end(),
						other.hepot[s].begin(), other.hepot[s].end());
		ljpot[s].insert(ljpot[s].end(),
						other.ljpot[s].begin(), other.ljpot[s].end());
	}
#endif // PREHYDROPHOBIC
}
    
//TODO should use std::move here
inline MeshInstance::MeshInstance(MeshInstance&& other) : 
        instance_id(other.instance_id), 
        silent(other.silent), 
        Dprotein(other.Dprotein), 
        Dsolvent(other.Dsolvent), 
        xyz_offset(other.xyz_offset), 
        rotation(other.rotation), 
        radius(other.radius),
        mesh_ptr(other.mesh_ptr),
#ifndef PREHYDROPHOBIC
		psel(other.psel),
#endif // PREHYDROPHOBIC
#ifdef USING_MAXD  // maxd no longer required
		maxd(std::move(other.maxd)),	// Ok to move values
#endif // maxd
        quad_points_per_triangle(other.quad_points_per_triangle),
        qual_points_per_triangle(other.qual_points_per_triangle)
{
//NB emplace?
//TODO move constructor should do better than this - default?
	patches.insert(patches.end(), other.patches.begin(), other.patches.end());
	charges.insert(charges.end(), other.charges.begin(), other.charges.end());
#ifndef PREHYDROPHOBIC
	allCharges.insert(allCharges.end(),
					  other.allCharges.begin(), other.allCharges.end());
	for (int s = 0; s < 2; s++) {
		hepot[s].insert(hepot[s].end(),
						other.hepot[s].begin(), other.hepot[s].end());
		ljpot[s].insert(ljpot[s].end(),
						other.ljpot[s].begin(), other.ljpot[s].end());
	}
#endif // PREHYDROPHOBIC
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

