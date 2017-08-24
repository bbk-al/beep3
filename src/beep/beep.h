/* Author: david fallaize		Created on: 23 Jul 2010 */

/*! \file beep.h
 * \brief This module declares the BEEP class for boundary element
 *	electrostatic potential calculations.
 *
 *	Two types of usage are anticipated:
 *	-- file-based:  a full configuration is written to a file and the
 *	   energy for that configuration is calculated.
 *	-- programmatic:  a program is written that directly interacts with BEEP
 *	   through the interfaces available.
 *	Calculation times will generally dominate overall processing times.
 *
 *	Example use:
 *	-- file-based: configuration files produced for each required calculation.
 *		Run either by BEEP executable or from a simple program:
 *		BEEP b(config, quad_points, qual_points, nbsize, kappa, use_planar);
 *		b.solve(1e-6, 100);
 *		b.calculate_energies();
 *		b.calculate_forces();
 *		b.write_fh(file);
 *	-- programmatic: custom use of the "API" to maintain the scenario, such as:
 *		BEEP b(solvent_dielectric, screening_kappa, quads, quals,
 *		       nbsize, bool planar);
 *		b.load_library_mesh(mtz);
 *		b.insert_mesh_instance(...);
 *		b.move_mesh_instance(...);
 *		... as with static case
 *
 *	This module contains the following public classes:
 *	- BEEP -- provides for scenario construction from configuration
 *		files and calculation of electrostatic	potentials and
 *		overall energy.
 *	- RunInfo -- provides calculation scales and timing data.
 */

#ifndef BEEP_H_
#define BEEP_H_

#include "../common/math_vector.h"
#include <vector>
#include <fstream>
#include "../bem/mesh.h"
#include "../bem/mesh_instance.h"
#include "../bem/local_integrations.h"
#include "../bem/config_file.h"
#include <boost/shared_ptr.hpp>
#include "../fmm/fmm_bem_hybrid.h"
#include "../fmm/fmm_time_info.h"
#include "../bem/opendx_utils.h"

// these set the number of qualocation and quadrature points for each
// triangle of the node-patch : this controls the discretization and level
// of detail in numerical integrations.  The settings for these two can
// be overriden in the xml config file using quad= or qual= in the <BEEP>
// tag (for global settings) or for individual mesh instances using quad=
// and qual= in the <instance> tags.
#define DEFAULT_QUADS 0
#define DEFAULT_QUALS 4

#ifdef OPENCL
#include "../opencl/opencl_handler.h"
#include "../bem/opencl_bem.h"
//#define CACHE_GPU_RESULTS
#define BEM_EXPLICIT_CHUNKSIZE 1024*1024
#define MAX_FMM_SIZE 250000
#define DEFAULT_BEM_NEIGHBOURHOOD_SIZE 3000

#else

#define CACHE_GPU_RESULTS
#define BEM_EXPLICIT_CHUNKSIZE 1024*1024*1024
#define MAX_FMM_SIZE 80
#define DEFAULT_BEM_NEIGHBOURHOOD_SIZE 1500

#endif

/*! \class RunInfo
 * \brief Provides calculation size and timing data.
 */
class RunInfo
{
    
public:

    inline std::string str() const
    {
        std::ostringstream buf;
        buf << num_patches << "," 
            << neighbourhood_size << "," 
            << time_total << "," 
            << time_rhs << "," 
            << time_gmres << ","
            << time_precalcs << "," 
            << num_explicits << ","
            << num_iterations <<  ","
            << bem_fmm_per_iter << "," 
            << rhs_fmm << "\n";
        return buf.str();
    }
    
    fmm::TimeInfo bem_fmm_per_iter;
    fmm::TimeInfo rhs_fmm;
    unsigned int num_patches;
    unsigned int neighbourhood_size;
    unsigned int num_iterations;
    double time_total;
    double time_rhs;
    double time_gmres;
    double time_precalcs;
    size_t num_explicits;
    
};

/*! \class BEEP
 * \brief Provides the top level molecular scenario definition and solver.
 */
class BEEP
{

	//typedef fmm::FMM_BEM_Octree_2TERMS FMM_BEM_Octree;
	//typedef fmm::FMM_BEM_Octree_4TERMS FMM_BEM_Octree;
	// 3FIGURE ACCURACY:
    typedef fmm::FMM_BEM_Octree<9,9,67,BasicNodePatch> FMM_BEM_Octree;
	//6FIGURE ACCURACY:
	//typedef fmm::FMM_BEM_Octree<18,18,300,BasicNodePatch> FMM_BEM_Octree;

public:
	//! Constructor from configuration for static programs
    BEEP(const ConfigFile& config,
         int cmdline_quad_points,
         int cmdline_qual_points,
         int cmdline_nbsize,
         double cmdline_kappa,
         bool force_planar);

	//! Constructor from named configuration file for static programs
    BEEP(const std::string& config_filename, bool planar, bool read_fh);

	//! Constructor for dynamic programs
    BEEP(double solvent_dielectric, double screening_kappa,
	     int cmdline_quads, int cmdline_quals,
	     int cmdline_nbsize, bool planar) :
	  manually_set_bounding_cube(false), 
	  beta0(1e-10),
	  Dsolvent(solvent_dielectric), 
	  kappa(screening_kappa), 
	  cmdline_bem_nbsize(cmdline_nbsize), 
	  cmdline_quad_points_per_triangle(cmdline_quads), 
	  cmdline_qual_points_per_triangle(cmdline_quals), 
	  cmdline_kappa(screening_kappa), 
	  force_planar(planar),
	  resolve(true),
	  previd(1),	// Will be invalid on BEEP creation (as required)
	  prevloc(Vector(0,0,0)),
	  unrot(Quaternion(1,0,0,0)),
	  skipping_precalcs(false)
	{ }

//TODO: NB neither default nor copy constructor (others?) cover these:
#if 0
    Vector bounding_cube_centre;
    double bounding_cube_edge_length;
    std::vector<LintArray_Size> local_integrations;
    std::unique_ptr<FMM_BEM_Octree> fmm;   // init_fmm_octree
    std::unique_ptr<fmm::FMM_Octree_6FIG_ACCURACY> rhs_octree;  // set_rhs
    fmm::TimeInfo vanilla_fmm_timer;
    fmm::TimeInfo bem_fmm_timer;
#ifdef OPENCL
    OpenCL_Handler global_ocl_handler;
#endif
#endif  // if 0
    
	//! Copy constructor
	BEEP(const BEEP& other) : 
	  manually_set_bounding_cube(false), 
	  beta0(other.beta0), 
	  Dsolvent(other.Dsolvent), 
	  kappa(other.kappa), 
	  cmdline_bem_nbsize(other.cmdline_bem_nbsize), 
	  cmdline_quad_points_per_triangle(other.cmdline_quad_points_per_triangle), 
	  cmdline_qual_points_per_triangle(other.cmdline_qual_points_per_triangle),
	  cmdline_kappa(other.cmdline_kappa), 
	  force_planar(other.force_planar),
	  resolve(other.resolve),
	  previd(other.previd),
	  prevloc(other.prevloc),
	  unrot(other.unrot),
	  skipping_precalcs(false),
	  meshes(other.meshes)
	{}
    
	//NB Destructor is default anyway
    // ~BEEP() {};
    
//Python interface methods:
	//! Initialise from configuration object
    void init(const ConfigFile& config);
 
	//! Load prepared meshes into library from mesh-tar-zip
    Mesh& load_library_mesh(const std::string& mtz_filename);

	//! Access library mesh
	Mesh& get_library_mesh(unsigned int id) { return meshes.getMesh(id); }

	//! Clear mesh instances in range [start, end)
    void clear_mesh_instances
		(unsigned int start = 0, int end = -1);

	//! Create a new mesh instance from library mesh
    MeshInstance& insert_mesh_instance(
		unsigned int id,			//! Library (reference) mesh id
	    const Vector& location,		//! Offset vector to apply to PDB structure
		const Quaternion& rotation,	//! Rotation to apply to PDB structure
	    double dielectric);			//! Protein dielectric to apply
	//! Move a mesh instance
	MeshInstance& move_mesh_instance(
		unsigned int instance_id,
	    const Vector& location,
		const Quaternion& rotation,
		double dielectric);
	//! Check if a point is inside a MeshInstance (other than skip_id)
	int get_instance_id(const Vector& pt, int skip_id = -1) const;
	//! Create a kinemage file
	//! \param filename string name of file to write to
    void create_kinemage(const std::string& filename,
	                     double fscale, double hscale, int num_colours,
	                     const std::string& preamble) const;
#ifdef FUTURE
	//! Preferred version - scales are worked out by this call
    void create_kinemage(const std::string& filename,
	                     int num_colours, const std::string& preamble) const;
#endif // FUTURE
    
    RunInfo solve(double gmres_stop_criteria, int max_iterations);
    std::string benchmark();

    double calculate_energies();
    void calculate_forces();
	size_t reset_fh_vals();  // returns as get_total_patches()
    void reset_library_fh_vals();
    void write_fh(const std::string& output_filename);

    inline size_t get_total_patches() const;
    unsigned int get_bem_neighbourhood_size() const;
    inline void set_bem_neighbourhood_size(unsigned int nbsize) {
		cmdline_bem_nbsize = static_cast<int>(nbsize);
	}
    const BasicNodePatch& get_patch(size_t idx) const {
        return *get_patch_ptr(idx);
    }

    size_t get_num_mesh_instances() const { return meshes.size(); }
    inline MeshInstance& get_mesh_instance(size_t idx) const {
		return *(meshes[idx]);
	}
    void set_bounding_cube(const Vector& new_bc_centre, double new_edge_length);
    double get_potential_at_point(const Vector& pt) const;
    
    void write_opendx_xyz(const std::string& filename, 
                          unsigned int x, unsigned int y, unsigned int z,
                          const Vector& edges, const Vector& centre) const;
    void write_opendx(const std::string& filename,
	                  unsigned int x, unsigned int y, unsigned int z) const;
    void write_matrix() const;
    inline void skip_precalcs() { skipping_precalcs = true; }
    inline void set_quad_points_per_triangle(unsigned int quad_points);
    inline void set_qual_points_per_triangle(unsigned int qual_points);

//Not Python but still public:

    long init_fmm_octree();
    long precalc_neighbour_interactions();
    void sanity_check();
    unsigned int gmres(double residual_norm_stop_criteria, int max_iterations);
    long do_bem_fmm_iteration(unsigned int, double[], double[]);

    long set_rhs();
    
    void write_opendx_grid(const std::string& filename, const GridParms& grid)
		const;

private:

    size_t count_explicit_integrations() const;
    double calc_local_neighbourhood_memory_requirement() const;
    void calculate_bounding_cube();

    inline boost::shared_ptr<BasicNodePatch> get_patch_ptr(size_t idx) const;

    void evaluate_local_neighbours(double f_results[], double h_results[]);
    
    inline void wait_for_local_interactions();

    bool manually_set_bounding_cube;
    Vector bounding_cube_centre;
    double bounding_cube_edge_length;

    const double beta0;

    //MeshList mesh_library;
    MeshInstanceList meshes;

	// reset as needed in solve and benchmark (indirectly)
    std::vector<LintArray_Size> local_integrations;

    std::unique_ptr<FMM_BEM_Octree> fmm;   // init_fmm_octree
    std::unique_ptr<fmm::FMM_Octree_6FIG_ACCURACY> rhs_octree;  // set_rhs

	// reset in solve (and benchmark):
    boost::shared_array<double> f_rhs;  // set in set_rhs
    boost::shared_array<double> h_rhs;  // set in set_rhs
    boost::shared_array<double> f_lhs;	// set in gmres (and benchmark)
    boost::shared_array<double> h_lhs;	// set in gmres (and benchmark)

	//NB These are constants!
    double Dsolvent;  // Solvent dielectric
    double kappa;	  // Inverse Debye screening length
    
    // command line parameters -- these override xml config settings
    int cmdline_bem_nbsize;
    int cmdline_quad_points_per_triangle;
    int cmdline_qual_points_per_triangle;
    double cmdline_kappa;
    bool force_planar;
    bool skipping_precalcs;

	// Flag to track whether solve() is required
	bool resolve;
	// Reversion data
	unsigned int previd;
	Vector prevloc;
	Quaternion unrot;

	// Timers
    fmm::TimeInfo vanilla_fmm_timer;
    fmm::TimeInfo bem_fmm_timer;

#ifdef OPENCL
    OpenCL_Handler global_ocl_handler;
#endif
};

// Globally inlined functions
//NB some of these look improbable!

inline unsigned int BEEP::get_bem_neighbourhood_size() const { 
	unsigned int nbsize = (cmdline_bem_nbsize != -1)
		? (static_cast<unsigned int>(cmdline_bem_nbsize))
		: DEFAULT_BEM_NEIGHBOURHOOD_SIZE;
	return nbsize;
}

inline size_t BEEP::get_total_patches() const {
	return meshes.get_total_patches();
}

inline boost::shared_ptr<BasicNodePatch> BEEP::get_patch_ptr(size_t idx) const
{
	unsigned int ctr=0;
	for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator
			it=meshes.cbegin(), end=meshes.cend();
		 it != end; ++it)
	{
		const MeshInstance& m = **it;
		ctr += m.get_num_node_patches();

		if (ctr > idx) {
			unsigned int idx_offset = m.get_num_node_patches() - (ctr - idx);
			return m.get_node_patches()[idx_offset];
		}
	}

	// should not get here!
	throw std::exception();
	//return NodePatch();
}

inline void BEEP::set_quad_points_per_triangle(unsigned int quad_points) { 
	cmdline_quad_points_per_triangle = quad_points;
	
	// enforce new value in mesh instances and mesh library instances
	for (unsigned int ii=0; ii < meshes.size(); ++ii) {
		meshes[ii]->set_quad_points_per_triangle(quad_points);
	}
}

inline void BEEP::set_qual_points_per_triangle(unsigned int qual_points) { 
	cmdline_qual_points_per_triangle = qual_points; 
	
	// enforce new value in mesh instances and mesh library instances
	for (unsigned int ii=0; ii < meshes.size(); ++ii) {
		meshes[ii]->set_qual_points_per_triangle(qual_points);
	}
}

inline void BEEP::wait_for_local_interactions() {

#ifdef OPENCL
	bool precalcs_completed = (global_ocl_handler.pending() == 0);
	if (!precalcs_completed) {
		//std::cout << "Waiting for local interactions (GPU precalcs) to complete" << std::endl;
		while (true) {
			if (global_ocl_handler.pending() == 0) break;
			boost::this_thread::sleep(boost::posix_time::microseconds(1));
		}
		//std::cout << "Local interaction pre-calcs are complete" << std::endl;
	}
	precalcs_completed = true;
#endif  // OPENCL
}

#endif /* BEEP_H_ */

