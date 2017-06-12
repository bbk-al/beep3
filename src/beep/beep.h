/*
* beep.h
*
*  Created on: 23 Jul 2010
*      Author: david
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

// these set the number of qualocation and quadrature points for each triangle of
// the node-patch : this controls the discretization and level of detail in numerical
// integrations.  The settings for these two can be overriden in the xml config file
// using quad= or qual= in the <BEEP> tag (for global settings) or for individual 
// mesh instances using quad= and qual= in the <instance> tags.
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

class BEEP {

    //typedef fmm::FMM_BEM_Octree_2TERMS FMM_BEM_Octree;
    //typedef fmm::FMM_BEM_Octree_4TERMS FMM_BEM_Octree;
    typedef fmm::FMM_BEM_Octree<9,9,67,BasicNodePatch> FMM_BEM_Octree; // 3FIGURE ACCURACY
    //typedef fmm::FMM_BEM_Octree<18,18,300,BasicNodePatch> FMM_BEM_Octree; //6FIGURE ACCURACY

public:

    BEEP(double solvent_dielectric, double screening_kappa, int cmdline_quads, int cmdline_quals, int cmdline_nbsize, bool planar) :
       manually_set_bounding_cube(false), 
       beta0(1e-10), Dsolvent(solvent_dielectric), 
       kappa(screening_kappa), 
       cmdline_bem_nbsize(cmdline_nbsize), 
       cmdline_quad_points_per_triangle(cmdline_quads), 
       cmdline_qual_points_per_triangle(cmdline_quals), 
       cmdline_kappa(screening_kappa), 
       force_planar(planar),
       skip_precalcs(false)
    {}
    
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
      skip_precalcs(false)
    {
        mesh_library.insert(mesh_library.end(), other.mesh_library.begin(), other.mesh_library.end());
        meshes.insert(meshes.end(), other.meshes.begin(), other.meshes.end());
    }
    
    BEEP(const ConfigFile& config, 
         int cmdline_quad_points, 
         int cmdline_qual_points,
         int cmdline_nbsize,
         double cmdline_kappa,
         bool force_planar);
         
    BEEP(const std::string& config_filename, bool planar, bool read_fh);
    ~BEEP() {};
    
    void write_matrix() const;
    void init(const ConfigFile& config);
    
    Mesh& py_load_library_mesh(const std::string& mtz_filename);
    void py_clear_mesh_instances();
    MeshInstance& py_insert_mesh_instance(unsigned int, const Vector&, const Quaternion&, double dielectric);
    void py_kinemage(const std::string& filename, double fscale, double hscale, int num_colours, const std::string& preamble) const;
    
    RunInfo solve(double gmres_stop_criteria, int max_iterations);
    std::string benchmark();
    long set_rhs();

    long init_fmm_octree();
    long precalc_neighbour_interactions();
    void sanity_check();
    unsigned int gmres(double residual_norm_stop_criteria, int max_iterations);
    long do_bem_fmm_iteration(unsigned int, double[], double[]);
    double calculate_energies();
    void calculate_forces();
    void write_fh(const std::string& output_filename);
    void reset_library_fh_vals();
    inline void hack_skip_precalcs() { skip_precalcs = true; }

    inline const BasicNodePatch& get_patch(size_t idx) const
    {
        return *get_patch_ptr(idx);
    }

    inline size_t get_total_patches() const
    {
        unsigned int total_np=0;
        for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator it=meshes.begin(), end=meshes.end(); it != end; ++it)
        {
            const MeshInstance& m = **it;
            total_np += m.get_num_node_patches();
        }
        return total_np;
    }
    
    inline size_t get_num_mesh_instances() const { return meshes.size(); }
    inline MeshInstance& get_mesh_instance(size_t idx) const { return *(meshes[idx]); }
    void set_bounding_cube(const Vector& new_bc_centre, double new_edge_length);
    double get_potential_at_point(const Vector& pt) const;
    
    void write_opendx_xyz(const std::string& filename, 
                          unsigned int x, 
                          unsigned int y, 
                          unsigned int z,
                          const Vector& edges,
                          const Vector& centre) const;
    void write_opendx(const std::string& filename, unsigned int x, unsigned int y, unsigned int z) const;
    void write_opendx_grid(const std::string& filename, const GridParms& grid) const;

    inline void set_bem_neighbourhood_size(unsigned int nbsize) { cmdline_bem_nbsize = static_cast<int>(nbsize); }
    inline unsigned int get_bem_neighbourhood_size() const { 
        unsigned int nbsize = (cmdline_bem_nbsize != -1) ? (static_cast<unsigned int>(cmdline_bem_nbsize)) : DEFAULT_BEM_NEIGHBOURHOOD_SIZE;
        return nbsize;
    }
    
    inline void set_quad_points_per_triangle(unsigned int quad_points) { 
        cmdline_quad_points_per_triangle = quad_points;
        
        // enforce new value in mesh instances and mesh library instances
        for (unsigned int ii=0; ii < meshes.size(); ++ii) {
            meshes[ii]->set_quad_points_per_triangle(quad_points);
        }
    }
    
    inline void set_qual_points_per_triangle(unsigned int qual_points) { 
        cmdline_qual_points_per_triangle = qual_points; 
        
        // enforce new value in mesh instances and mesh library instances
        for (unsigned int ii=0; ii < meshes.size(); ++ii) {
            meshes[ii]->set_qual_points_per_triangle(qual_points);
        }
    }
    
private:

    size_t count_explicit_integrations() const;
    double calc_local_neighbourhood_memory_requirement() const;
    void calculate_bounding_cube();

    inline boost::shared_ptr<BasicNodePatch> get_patch_ptr(size_t idx) const
    {
        unsigned int ctr=0;
        for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator it=meshes.begin(), end=meshes.end(); it != end; ++it)
        {
            const MeshInstance& m = **it;
            ctr += m.get_num_node_patches();

            if (ctr > idx)
            {
                unsigned int idx_offset = m.get_num_node_patches() - (ctr - idx);
                return m.patches[idx_offset];
            }

        }

        // should not get here!
        throw std::exception();
        //return NodePatch();
    }

    void evaluate_local_neighbours(double f_results[], double h_results[]);
    
    inline void wait_for_local_interactions()
    {

#ifdef OPENCL
        bool precalcs_completed = (global_ocl_handler.pending() == 0);
        if (!precalcs_completed)
        {
            //std::cout << "Waiting for local interactions (GPU precalcs) to complete" << std::endl;
            while (true)
            {
                if (global_ocl_handler.pending() == 0) {
                    break;
                }
                boost::this_thread::sleep(boost::posix_time::microseconds(1));
            }
            //std::cout << "Local interaction pre-calcs are complete" << std::endl;
        }
        precalcs_completed = true;
#endif

        return;
    }

    bool manually_set_bounding_cube;
    Vector bounding_cube_centre;
    double bounding_cube_edge_length;

    const double beta0;

    MeshList mesh_library;
    MeshInstanceList meshes;

    std::vector<LintArray_Size> local_integrations;

    std::auto_ptr<FMM_BEM_Octree> fmm;
    std::auto_ptr<fmm::FMM_Octree_6FIG_ACCURACY> rhs_octree;

    boost::shared_array<double> f_rhs;
    boost::shared_array<double> h_rhs;
    boost::shared_array<double> f_lhs;
    boost::shared_array<double> h_lhs;

    double Dsolvent;
    double kappa;
    
    // command line parameters -- these override xml config settings
    int cmdline_bem_nbsize;
    int cmdline_quad_points_per_triangle;
    int cmdline_qual_points_per_triangle;
    double cmdline_kappa;
    bool force_planar;

    fmm::TimeInfo vanilla_fmm_timer;
    fmm::TimeInfo bem_fmm_timer;
    bool skip_precalcs;

#ifdef OPENCL
    OpenCL_Handler global_ocl_handler;
#endif

};

#endif /* BEEP_H_ */

