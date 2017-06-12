/*
* node_patch.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef NODE_PATCH_H_
#define NODE_PATCH_H_

#include "../common/math_vector.h"
#include "../common/octree_indexer.h"
#include "spharm.h"
#include "quad_point.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "bezier.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/thread/mutex.hpp>

extern double fGeometricCorrection(double geometry_correction, double epsilon);
extern double hGeometricCorrection(double geometry_correction, double epsilon);

// fwd decl.
class Mesh;
class MeshInstance;
class Vertex;
class Triangle;

//#define CENTROID_COLLOCATION

class BasicNodePatch : public Vector
{

public:

    friend class Mesh;

    BasicNodePatch() : Vector(0,0,0), mesh_ptr(NULL), idx(0), vertex_idx(0), quad_points_per_triangle(0), qual_points_per_triangle(0), galerkin_points_per_triangle(0) {
        mutex_ptr = new boost::mutex;
        
    }
    virtual ~BasicNodePatch() {
        delete mutex_ptr;
    }
    
    BasicNodePatch(const BasicNodePatch& other);
    BasicNodePatch(unsigned int vertex_idx, const Mesh& mesh);

    // const accessors for private member data
    inline void set_dielectric_ratio(double val) { dielectric_ratio = val; }
    inline double get_dielectric_ratio() const { return dielectric_ratio; }
    inline const Vector& get_node() const { return node; }
    inline const Vector& get_centroid() const { return centroid; }
    inline const Vector& get_normal() const { return normal; }
    inline const Vector& get_alt_normal() const { return alt_normal; }
    inline double get_planar_area() const { return planar_area; }
    inline double get_bezier_area() const { return bezier_area; }
    
    virtual const Mesh& get_ref_mesh() const  { return *mesh_ptr; }
    
    inline void obtain_shared_quad_ptrs(std::vector<boost::shared_ptr<QuadList> >& dump_here) const { 
        dump_here.push_back(get_quadrature_points());
        dump_here.push_back(get_qualocation_points());
        return;
    }

    inline void set_quad_points_per_triangle(unsigned int quad_points) { quad_points_per_triangle = quad_points;  }
    inline void set_qual_points_per_triangle(unsigned int qual_points) { qual_points_per_triangle = qual_points; }
    inline void set_galerkin_points_per_triangle(unsigned int galerkin_points) { galerkin_points_per_triangle = galerkin_points; }

    inline boost::shared_ptr<QuadList> get_single_qual_pt() const { return single_qual_pt; }
    inline boost::shared_ptr<QuadList> get_single_quad_pt() const { return single_quad_pt; }

    boost::shared_ptr<QuadList> generate_points(const QuadratureRule& rule, int subdivides=0) const;

    inline static const QuadratureRule& get_rule(int order)
    {
        const QuadratureRule* rule;
        switch (order)
        {
        case(1):
            rule = &(BasicTriangle::gauss_legendre_pts_1());
            break;
        case(4):
            rule = &(BasicTriangle::gauss_legendre_pts_4());
            break;
        case(6):
            rule = &(BasicTriangle::triangle_order_6());
            break;
        case(7):
            rule = &(BasicTriangle::gauss_legendre_pts_7());
            break;
        case(16):
            rule = &(BasicTriangle::gauss_legendre_pts_16());
            break;
        default:
            std::cerr << "Bad rule specified: " << order << std::endl;
            throw std::exception();
        }
        return *rule;
    }
    
    inline const QuadratureRule& get_qualocation_rule() const { return get_rule(qual_points_per_triangle); }
    inline const QuadratureRule& get_quadrature_rule() const  { return get_rule(quad_points_per_triangle); }
    inline const QuadratureRule& get_galerkin_rule() const  { return get_rule(galerkin_points_per_triangle); }
    
    inline boost::shared_ptr<QuadList> get_qualocation_points() const 
    { 
   
        // create a shared pointer to the existing list (*if* it exists)
        boost::shared_ptr<QuadList> qp_ptr = weak_qual_ptr.lock();
        if (!qp_ptr)
        {
            // get mutex lock then recheck weak pointer
            boost::mutex::scoped_lock lock(*mutex_ptr);
            qp_ptr = weak_qual_ptr.lock();
            if (qp_ptr) { return qp_ptr; }
            
            // if qual rule is zero, return a qual-point for the whole node patch
            if (qual_points_per_triangle == 0) { 
                weak_qual_ptr = single_qual_pt; 
                return single_qual_pt;
            }
            
            // if the weak pointer is empty/expired then need to regenerate
            // a list of quallocation points
            qp_ptr = generate_points(get_qualocation_rule());
            weak_qual_ptr = qp_ptr;
            qual_pts = qp_ptr; // permanently store it
            
            // renormalize quallocation weights: they should sum to 1.0
            // not the total area
            for (int ii=0; ii < qual_pts->size(); ++ii)
            {
                QuadPoint& qp = (*qual_pts)[ii];
                qp.weight() = qp.weight() / get_bezier_area();
            }
        }
        
        return qp_ptr;
        
    }  
        
    inline boost::shared_ptr<QuadList> get_quadrature_points() const
    {
        // create a shared pointer to the existing list (*if* it exists)
        boost::shared_ptr<QuadList> qp_ptr = weak_quad_ptr.lock();
        if (!qp_ptr)
        {
            
            // get mutex lock then recheck weak pointer
            boost::mutex::scoped_lock lock(*mutex_ptr);
            qp_ptr = weak_quad_ptr.lock();
            if (qp_ptr) { return qp_ptr; }

            // if quad rule is zero, return a quad-point for the whole node patch
            if (quad_points_per_triangle == 0) { 
                weak_quad_ptr = single_quad_pt; 
                return single_quad_pt;
            }

            // if the weak pointer is empty/expired then need to regenerate
            // a list of quallocation points
            qp_ptr = generate_points(get_quadrature_rule());
            weak_quad_ptr = qp_ptr;
        }
        
        return qp_ptr;
    }

    inline  boost::shared_ptr<QuadList> get_galerkin_points() const
    {
        // create a shared pointer to the existing list (*if* it exists)
        boost::shared_ptr<QuadList> qp_ptr = weak_galerkin_ptr.lock();
        if (!qp_ptr)
        {
            // get mutex lock then recheck weak pointer
            boost::mutex::scoped_lock lock(*mutex_ptr);
            qp_ptr = weak_galerkin_ptr.lock();
            if (qp_ptr) { return qp_ptr; }
        
            // if the weak pointer is empty/expired then need to regenerate
            // a list of quallocation points
            qp_ptr = generate_points(BasicTriangle::gauss_legendre_pts_4(),0);
            weak_galerkin_ptr = qp_ptr;
        }
        
        return qp_ptr;
    }
    
    // pure virtual as it's very important that the derived classes handle any changes in 
    // coordinate frame correctly
    virtual void change_coordinate_frame(QuadList& qps) const { return; }
    
    inline boost::shared_ptr<QuadList> get_energy_quad_points() const { return get_galerkin_points(); }
    
    void get_edge_points(std::vector<PointNormal>& pts) const;

    // Python- not so great at automagically handling C++ style inheritance :-(
    inline Vector& py_vector() { return static_cast<Vector&>(*this); }

    inline double get_charge(unsigned short charge_idx) const
    {
        return get_fmm_charge_equiv(charge_idx);
    }
    
    Vector calculate_force(double f,
                            double h,
                            double epsilon_int,
                            double epsilon_ext,
                            double kappa) const;
    
    double dielectric_ratio;
    double f;
    double h;
    double energy_coefficient_f;
    double energy_coefficient_h;
    Vector force_coefficient_f;
    Vector force_coefficient_h;
    double gc;

    inline unsigned int get_idx() const { return idx; }
    inline void set_idx(unsigned int new_idx) { idx = new_idx; }
    
    void calc_solid_angle();
    
protected:

    inline double get_fmm_charge_equiv(unsigned short equiv_charge_idx) const
    {

        switch(equiv_charge_idx)
        {
        case 0:
            return bezier_area*f*normal.x;
            break;
        case 1:
            return bezier_area*f*normal.y;
            break;
        case 2:
            return bezier_area*f*normal.z;
            break;
        case 3:
            return bezier_area*h;
            break;
        
        // epsilon weighted values
        case 4:
            return bezier_area*f*normal.x / dielectric_ratio;
            break;
        case 5:
            return bezier_area*f*normal.y / dielectric_ratio;
            break;
        case 6:
            return bezier_area*f*normal.z / dielectric_ratio;
            break;
        case 7:
            return bezier_area*h / dielectric_ratio;
            break;

        // epsilon weighted values for kappa=0 (Gpt kernels) (and non-epsilon weighted h vals)
        case 8:
            return bezier_area*f*normal.x / dielectric_ratio;
            break;
        case 9:
            return bezier_area*f*normal.y / dielectric_ratio;
            break;
        case 10:
            return bezier_area*f*normal.z / dielectric_ratio;
            break;
        case 11:
            return bezier_area*h;
            break;

            
        default:
            // should not get here
            throw std::exception();
        }

        return 0.0;
    }

    const Mesh* mesh_ptr;
    unsigned int idx; // unique number of this node patch
    unsigned int vertex_idx; // vertex index within a Mesh object -- need to know this to build quadrature points on the fly
    Vector node;
    Vector centroid;
    Vector normal;
    Vector alt_normal;
    double bezier_area;
    double planar_area;
    
    boost::shared_ptr<QuadList> single_qual_pt;
    boost::shared_ptr<QuadList> single_quad_pt;

    unsigned int quad_points_per_triangle;
    unsigned int qual_points_per_triangle;
    unsigned int galerkin_points_per_triangle;
    
    
private:
    
    mutable boost::weak_ptr<QuadList> weak_quad_ptr;
    mutable boost::weak_ptr<QuadList> weak_qual_ptr;
    mutable boost::weak_ptr<QuadList> weak_galerkin_ptr;
    
    // permanently store qual_pts (more efficient when using
    // OpenCL to compute local integrations, as these are
    // requested repeatedly rather than just once; less efficient
    // if not using OpenCL, but not that huge an overhead).
    mutable boost::shared_ptr<QuadList> qual_pts;

    mutable boost::mutex* mutex_ptr;
    
};

// This NodePatch class is an actual NodePatch located in universe coordinates- it is 
// linked conceptually to a MeshInstance which holds the rotation/xyz_offset for the
// node patch.
class NodePatch : public BasicNodePatch
{
public:
    
    NodePatch(const BasicNodePatch& other, const MeshInstance& mesh_instance, unsigned int _idx=0);
    NodePatch(const NodePatch& other) : BasicNodePatch(other), minst(other.minst) {}; // copy c'tor
    virtual ~NodePatch() {}

    // gets the reference mesh via the underlying MeshInstance pointer
    virtual const Mesh& get_ref_mesh() const;
    virtual void change_coordinate_frame(QuadList& qps) const;
    
private:
    
    const MeshInstance& minst;
    
};

typedef std::vector< boost::shared_ptr<BasicNodePatch> > PatchList;


#endif /* NODE_PATCH_H_ */
