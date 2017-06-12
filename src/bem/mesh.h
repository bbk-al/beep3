/*
* mesh.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef MESH_H_
#define MESH_H_

#include <boost/shared_array.hpp>

#include "../common/math_vector.h"
#include "bem_kernels.h"
#include "node_patch.h"
#include "vertex.h"
#include "triangle.h"
#include "edge.h"
#include "../common/charge.h"
#include "spharm.h"
#include <exception>
#include <boost/shared_ptr.hpp>
#include "bezier.h"
#include "png1.h"
#include <fstream>
#include "mesh_tarball.h"

#ifdef __CHARMC__
#include <pup_stl.h>
#endif

class MeshInstance; // fwd decl.
class BEEP;

class Mesh
{

public:

    friend class MeshInstance;
    friend class BEEP;
    friend class BasicNodePatch;

    Mesh() :
        centre(Vector(0,0,0)),
        total_planar_area(0),
        total_bezier_area(0),
        net_charge(0),
        radius(0),
        done_energy_precalcs(false),
        quad_points_per_triangle(0),
        qual_points_per_triangle(0)
       {}

    Mesh(const std::string& fname, bool force_planar=false) : done_energy_precalcs(false) {
        init(fname, force_planar);
    }
    Mesh(const std::string& mesh_filename,
        const std::string& xyzq_filename,
        bool force_planar=false);

    // copy constructor
    Mesh(const Mesh& other);

    virtual ~Mesh();

#ifdef __CHARMC__
    virtual void pup(PUP::er &p) {

        p | node_patches;
        p | charges;
        p | centre;
        p | radius;

    }
#endif

    // init a Mesh object from tarball filename
    void init(const std::string&, bool force_planar=false);

    // exception that might be thrown if mesh cannot be init'd
    class MeshError : public std::exception {
    public:
        MeshError() : std::exception() {}
    };

    inline const std::vector<BasicNodePatch>& get_node_patches() const { return node_patches; }
    inline const std::vector<Charge>& get_charges() const { return charges; }
    inline std::vector<BasicNodePatch>& get_node_patches() { return node_patches; }
    inline std::vector<Charge>& get_charges() { return charges; }

    std::vector<Charge> get_ecm_charges(const std::string& ecm_filename, const Quaternion&, const Vector&) const;

    const BasicNodePatch& get_node_patch(unsigned int index) const
    {
        if (index >= node_patches.size()) { throw std::out_of_range("index out of range for mesh"); }

        return node_patches[index];
    }

    const Triangle& get_triangle(unsigned int index) const
    {
        if (index >= triangles.size()) { throw std::out_of_range("index out of range for mesh"); }

        return triangles[index];
    }

    const Charge& get_charge(unsigned int index) const
    {
        if (index >= charges.size()) { throw std::out_of_range("index out of range for mesh"); }

        return charges[index];
    }

    const Vertex& get_vertex(unsigned int index) const
    {
        if (index >= vertices.size()) { throw std::out_of_range("index out of range for mesh"); }

        return vertices[index];
    }


    // this function returns the centre and edge length of a cube such that all vertices defined in the Mesh lie within the cube.
    // (note that there is no rotational optimization involved here -- so this is probably not the smallest volume cube in which
    // the mesh can possibly fit.  This is meant to be used only to figure out the scaling required to fit this Mesh into a unit
    // Octree cube.)
    void get_bounding_cube(Vector &centre, double &edge_length) const
    {
        Vector max; // this will be the 'top right' corner
        Vector min; // and this will be 'bottom left'

        get_bounding_cube_limits(max, min);

        // figure out the maximum edge length in x/y/z dimension
        Vector diff = max - min;
        edge_length = diff.y > diff.x ? diff.y : diff.x;
        edge_length = diff.z > edge_length ? diff.z : edge_length;
        //edge_length = 1.0;

        // centre of the cube is the mid point of the two extremities
        centre = (max + min) / 2.0;

        return;
    }

    void get_bounding_cube_limits(Vector &max, Vector& min) const
    {
        // loop over node patch points in the mesh
        for (std::vector<BasicNodePatch>::const_iterator it=node_patches.begin(), end=node_patches.end();
            it != end;
            ++it)
        {
            const Vector& v = *it;
            if (it == node_patches.begin())
            {
                max = v;
                min = v;
            }
            else
            {
                max.x = v.x > max.x ? v.x : max.x;
                max.y = v.y > max.y ? v.y : max.y;
                max.z = v.z > max.z ? v.z : max.z;

                min.x = v.x < min.x ? v.x : min.x;
                min.y = v.y < min.y ? v.y : min.y;
                min.z = v.z < min.z ? v.z : min.z;
            }
        }
        
        // loop over charges in the mesh
        for (std::vector<Charge>::const_iterator it=charges.begin(), end=charges.end();
            it != end;
            ++it)
        {
            const Vector& v = *it;
            if (it == charges.begin() && node_patches.size() == 0)
            {
                max = v;
                min = v;
            }
            else
            {
                max.x = v.x > max.x ? v.x : max.x;
                max.y = v.y > max.y ? v.y : max.y;
                max.z = v.z > max.z ? v.z : max.z;

                min.x = v.x < min.x ? v.x : min.x;
                min.y = v.y < min.y ? v.y : min.y;
                min.z = v.z < min.z ? v.z : min.z;
            }
        }
        
        return;
    }

    void add_vertex(const Vector& v, const Vector& vn)
    {
        vertices.push_back(Vertex(v,vn));

    }

    void define_triangle(const unsigned int v1,
                        const unsigned int v2,
                        const unsigned int v3)
    {
        unsigned int t_idx = triangles.size();
        triangles.push_back(Triangle(vertices, v1,v2,v3));

        vertices[v1].add_triangle(t_idx);
        vertices[v2].add_triangle(t_idx);
        vertices[v3].add_triangle(t_idx);

        //const Triangle& t = triangles[t_idx];

        return;
    }

    inline unsigned int get_num_vertices() const { return vertices.size(); }
    inline unsigned int get_num_triangles() const { return triangles.size(); }
    inline unsigned int get_num_node_patches() const { return node_patches.size(); }
    inline unsigned int get_num_charges() const { return charges.size(); }
    inline unsigned int len() const { return get_num_node_patches(); }

    std::vector<Triangle>& get_triangles() { return triangles; }
    const std::vector<Triangle>& get_triangles() const { return triangles; }

    void calculate_vertex_normals();

    std::string kinemage_node_patches() const;
    std::string kinemage_fh_vals(double fscale, double hscale, int num_colours) const;

    void reorder_vertex_triangles()
    {
        // reorder the vertices for coherence
        for (std::vector<Vertex>::iterator it=vertices.begin(), end=vertices.end();
            it != end;
            ++it)
        {
            it->reorder(triangles);
        }
    }

    // create Bezier curved patches for each triangle element
    // These will be used to create quasi-curved node patches
    // CurvedTriangleType should be Triangle (for planar); HybridBezierPatch or PNG1_Triangle for curved
    template<typename CurvedTriangleType>
    void create_bezier_triangles()
    {
        triangle_ptrs.reserve(triangles.size());
        for (std::vector<Triangle>::const_iterator tri_it=triangles.begin(), tri_end=triangles.end(); tri_it != tri_end; ++tri_it)
        {
            triangle_ptrs.push_back(new CurvedTriangleType(*tri_it));
        }

    }

    // Create BasicNodePatch objects from the triangulated surface mesh
    void create_node_patches()
    {
        node_patches.reserve(vertices.size());

        double planar_area = 0;
        total_bezier_area = 0.0;
        for (std::vector<BasicTriangle*>::const_iterator it=triangle_ptrs.begin(), end=triangle_ptrs.end(); it != end; ++it)
        {
            planar_area += (**it).get_planar_area();
            total_bezier_area += (**it).get_area();
        }

        // instantiates umbrella-shaped Node Patches
        unsigned int ctr=0;
        for (std::vector<Vertex>::iterator it=vertices.begin(), end=vertices.end();
            it != end;
            ++it)
        {
            node_patches.push_back(BasicNodePatch(ctr++,*this));
        }

        std::cout << "Planar area: " << total_planar_area  << " vs. bezier area: " << total_bezier_area << std::endl;

        return;
    }

    Vector calculate_charges_centroid()
    {
        // This function sets the centre of the mesh by finding centroid of charges.
        // So obviously you should have defined the charges in the mesh first...!
        // (Better: use the init_centre() function to get the centre of the mesh
        // from an external file -- which lets you use e.g. hydropro to get a better
        // estimate of where the real centre of the molecule is, from a diffusey point
        // of view).
        Vector sum(0,0,0);
        if (charges.size() > 0)
        {
            // iterate over all charges and find centre
            for (std::vector<Charge>::const_iterator chit=charges.begin(), chend=charges.end(); chit != chend; ++chit)
            {
                sum += chit->position();
            }

            sum /= charges.size();
        }

        return sum;
    }

    Vector calculate_patches_centroid()
    {
        Vector centroid(0,0,0);
        for (std::vector<BasicNodePatch>::const_iterator np_it=node_patches.begin(), np_end=node_patches.end(); np_it != np_end; ++np_it)
        {
            const BasicNodePatch& np = *np_it;
            centroid += np.get_centroid() * np.get_bezier_area();            
        }
        
        if (total_bezier_area > 0)
        {
            centroid /= total_bezier_area;
        }
        else
        {
            std::cerr << "Error: Surface area of the mesh appears to be zero!" << std::endl;
            throw std::exception();
        }
        
        return centroid;
    }

    inline double get_radius() const { return radius; }

    double calculate_energy(double kappa,
                            double Dprotein,
                            double Dsolvent,
                            const double fvals[],
                            const double hvals[]) const;
    
    KahanVector calculate_qe_force(double Dprotein,
                              double Dsolvent,
                              const double fvals[],
                              const double hvals[]) const;
    void calculate_surface_integral_forces(double kappa,
                                           double Dprotein,
                                           double Dsolvent,
                                           const double fvals[],
                                           const double hvals[],
                                           KahanVector& MST_external,
                                           KahanVector& MST_internal,
                                           KahanVector& dbf,
                                           KahanVector& ionic) const;
                                             
    double py_calculate_energy(double kappa, double Dprotein, double Dsolvent) const;
    void py_calculate_forces(double kappa, double Dprotein, double Dsolvent, Vector& qE, Vector& MST_ext, Vector& MST_int, Vector& dbf, Vector& ionic) const;

    inline const Vector& get_centre() const { return centre; }
    inline void set_centre(const Vector& new_centre) {
        centre = new_centre;
    }

    inline void set_dielectric_ratio(double Dsolvent, double Dprotein)
    {
        double epsilon_ratio = Dsolvent/Dprotein;
        for (std::vector<BasicNodePatch>::iterator it=node_patches.begin(), end=node_patches.end(); it != end; ++it)
        {
            it->set_dielectric_ratio(epsilon_ratio);
        }
    }

    double calculate_volume() const;
    void init_energy_precalcs(); // must be public so that python module can use it

    inline void set_quad_points_per_triangle(unsigned int quad_points) { 
        quad_points_per_triangle = quad_points; 
        for (std::vector<BasicNodePatch>::iterator it=node_patches.begin(), end=node_patches.end(); it != end; ++it)
        {
            it->set_quad_points_per_triangle(quad_points);
        }
    }
    
    inline void set_qual_points_per_triangle(unsigned int qual_points) { 
        qual_points_per_triangle = qual_points; 
        for (std::vector<BasicNodePatch>::iterator it=node_patches.begin(), end=node_patches.end(); it != end; ++it)
        {
            it->set_qual_points_per_triangle(qual_points);
        }
    }
    
    inline unsigned int get_quad_points_per_triangle() const { return quad_points_per_triangle; }
    inline unsigned int get_qual_points_per_triangle() const { return qual_points_per_triangle; }
    
private:

    void init_mesh(const std::string& mesh_filename);
    void init_charges(const std::string& xyzq_filename);
    void init_fh_vals(const std::string& fh_filename);
    void init_centre(const std::string& centre_filename);
    void read_energy_precalcs(const std::string& energies_filename);
    void renormalise_energy_precalcs();
    

    double calculate_radius();

    // this has to be public so that the python module can mess with it
    Vector centre;

    std::vector<Triangle> triangles;
    std::vector<BasicTriangle*> triangle_ptrs;
    std::vector<Vertex> vertices;
    std::vector<BasicNodePatch> node_patches;
    std::vector<Charge> charges;

    double total_planar_area;
    double total_bezier_area;
    double net_charge;
    double radius;
    bool done_energy_precalcs;
    
    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;

};

typedef std::vector< boost::shared_ptr<Mesh> > MeshList;

#endif /* MESH_H_ */
