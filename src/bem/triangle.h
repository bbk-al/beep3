/*
* triangle.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../common/math_vector.h"
#include "quad_point.h"
#include <vector>
#include <iostream>
#include <ostream>
#include "vertex.h"
#include <gsl/gsl_linalg.h>

class BasicTriangle;
static const double ONE_THIRD = 1./3.;

Vector get_field(const BasicTriangle& tri, const Vector& fvals, const Vector& hvals, double u, double v);
Vector MST_from_E(const Vector& Eo,
                    const Vector& Ei,
                    const Vector& norm,
                    double kappa,
                    double f,
                    double dielectric);

double solve_delta(const Vector& Eo,
                   const Vector& Ei,
                   const Vector& norm);

Vector dbf_from_E(const Vector& Eo,
                    const Vector& Ei,
                    const Vector& norm,
                    double kappa,
                    double f,
                    double epsilon_int,
                    double epsilon_ext);
Vector dbf_from_MST(const Vector& Eo,
                    const Vector& Ei,
                    const Vector& norm,
                    double epsilon_int,
                    double epsilon_ext);
Vector dbf_from_delta(double delta,
                        const Vector& Ei,
                        const Vector& norm,
                        double epsilon_int,
                        double epsilon_ext);
                        
Vector ionic_from_E(const Vector& Eo,
                    const Vector& Ei,
                    const Vector& norm,
                    double kappa,
                    double f,
                    double epsilon_int,
                    double epsilon_ext);                    
void calculate_force_components(const BasicTriangle& tri,
                                Vector fvals,
                                Vector hvals,
                                double epsilon_int,
                                double epsilon_ext,
                                double kappa,
                                const QuadratureRule& rule,
                                unsigned int subdivides,
                                KahanVector& MST_external,
                                KahanVector& MST_internal,
                                KahanVector& dbf,
                                KahanVector& ionic);
Vector f_at_verts_to_E(const Vector& vertex1, const Vector& vertex2, const Vector& vertex3,
                       double vf1, double vf2, double vf3);

                            
static double tri6x[] = {0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877};
static double tri6y[] = {0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031};
static double tri6w[] = {0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667};

static double gl1x[] = {1./3.};
static double gl1y[] = {1./3.};
static double gl1w[] = {1.};

static double gl4x[] = {1./3., 0.2, 0.6, 0.2};
static double gl4y[] = {1./3., 0.6, 0.2, 0.2};
static double gl4w[] = {-0.56250000, 0.520833333333, 0.520833333333, 0.520833333333};

static double gl7x[] = {1./3., 0.05971587, 0.47014206, 0.47014206, 0.79742698, 0.10128650, 0.10128650};
static double gl7y[] = {1./3., 0.47014206, 0.05971587, 0.47014206, 0.10128650, 0.79742698, 0.10128650};
static double gl7w[] = {0.225, 0.132394153333333, 0.132394153333333, 0.132394153333333, 0.125939186666666, 0.125939186666666, 0.125939186666666};

static double gl16x[] = {0.333333333333333, 0.081414823414554, 0.459292588292723, 0.459292588292723, 0.65886138449648, 0.17056930775176, 0.17056930775176, 0.898905543365938, 0.050547228317031, 0.050547228317031, 0.008394777409958001, 0.263112829634638, 0.728492392955404, 0.263112829634638, 0.728492392955404, 0.008394777409958001};
static double gl16y[] = {0.333333333333333, 0.459292588292723, 0.459292588292723, 0.081414823414554, 0.17056930775176, 0.17056930775176, 0.65886138449648, 0.050547228317031, 0.050547228317031, 0.898905543365938, 0.263112829634638, 0.728492392955404, 0.008394777409958001, 0.008394777409958001, 0.263112829634638, 0.728492392955404};
static double gl16w[] = {0.144315607677787,0.09509163426728499,0.09509163426728499,0.09509163426728499,0.103217370534718,0.103217370534718,0.103217370534718,0.032458497623198,0.032458497623198,0.032458497623198,0.027230314174435,0.027230314174435,0.027230314174435,0.027230314174435,0.027230314174435,0.027230314174435};

class BasicTriangle
{
public:

    BasicTriangle(const Vector& _v1,
                  const Vector& _v2,
                  const Vector& _v3) : v1(_v1), v2(_v2), v3(_v3)
    {
        planar_area = calc_planar_area();
        planar_centroid = calc_planar_centroid();
        planar_normal = calc_planar_normal();
    }

    BasicTriangle(const BasicTriangle& other) : v1(other.v1), v2(other.v2), v3(other.v3) {
        planar_area = other.planar_area;
        planar_centroid = other.planar_centroid;
        planar_normal = other.planar_normal;
    }

    virtual ~BasicTriangle() {}

    virtual std::string type() const { return "BasicTriangle"; }

    // Quadrature Rules -- return static objects on first use 
    // (avoids the so-called static instantiation fiasco)
    static const QuadratureRule& triangle_order_6() {
        static QuadratureRule *t6 = new QuadratureRule(6, tri6x, tri6y, tri6w);
        return *t6;
    }
    static const QuadratureRule& gauss_legendre_pts_1() {
        static QuadratureRule *gl1 = new QuadratureRule(1, gl1x, gl1y, gl1w);
        return *gl1;
    }
    static const QuadratureRule& gauss_legendre_pts_4() {
        static QuadratureRule *gl4 = new QuadratureRule(4, gl4x, gl4y, gl4w);
        return *gl4;
    }
    static const QuadratureRule& gauss_legendre_pts_7() {
        static QuadratureRule *gl7 = new QuadratureRule(7, gl7x, gl7y, gl7w);
        return *gl7;
    }
    static const QuadratureRule& gauss_legendre_pts_16() {
        static QuadratureRule *gl16 = new QuadratureRule(16, gl16x, gl16y, gl16w);
        return *gl16;
    }

    inline const Vector& get_planar_centroid() const { return planar_centroid; }
    inline double get_planar_area() const { return planar_area; }
    inline const Vector& get_planar_normal() const { return planar_normal; }

    inline Vector get_parametric_centre() const { return get_point(ONE_THIRD, ONE_THIRD); }

    double solid_angle(const Vector& origin) const {

        Vector a = v1 - origin;
        Vector b = v2 - origin;
        Vector c = v3 - origin;

        double det = fabs(a.dot(b.cross(c)));

        double al = a.length();
        double bl = b.length();
        double cl = c.length();

        double div = al*bl*cl + a.dot(b)*cl + a.dot(c)*bl + b.dot(c)*al;
        double at = atan2(det, div);
        if(at < 0) at += pi; // If det>0 && div<0 atan2 returns < 0, so add pi.
        double omega = 2.0 * at;

        return omega;
    }

    inline const Vector& get_v1() const { return v1; }
    inline const Vector& get_v2() const { return v2; }
    inline const Vector& get_v3() const { return v3; }

    virtual double get_area() const {
        return planar_area;
    }

    std::string kinemage() const {
        std::ostringstream buf;
        buf << "X " << v1.x << " " << v1.y << " " << v1.z << " " << 
                       v2.x << " " << v2.y << " " << v2.z << " " << 
                       v3.x << " " << v3.y << " " << v3.z << "\n";
        return buf.str();
    }

    virtual double get_jacobian(double u, double v) const 
    {
        return 2.0*planar_area;
    }

    virtual Vector get_point(double u, double v) const {
        //std::cout << "BasicTriangle::get_point" << std::endl;
        double w = 1.0 - u - v;
        return w*v1 + u*v2 + v*v3;
    }

    virtual Vector get_normal(double u, double v) const {
        return planar_normal;
    }

    inline PointNormal get_point_normal(double u, double v) const 
    {
        return PointNormal(get_point(u,v), get_normal(u,v));
    }

    virtual Vector get_dU(double u, double v, double w) const
    {
        Vector dU = (v2 -v3*0.5 -v1*0.5);
        return dU;
    }

    virtual Vector get_dV(double u, double v, double w) const
    {
        Vector dV = (v3 -v1*0.5 -v2*0.5);
        return dV;
    }

    template <typename QuadPointType>
    void get_quad_points(std::vector<QuadPointType>& quad_points, const QuadratureRule& rule) const
    {
        // this function can throw exceptions if the triangle is so slender as to produce
        // NaNs in Quadrature Points.  So don't commit changes to the quad_points until
        // all succeed.
        std::vector<QuadPointType> tmp_quad_points;
    
        for (int qp_ctr=0; qp_ctr < rule.num_points; ++qp_ctr)
        {
            double u = rule.quad_points_x[qp_ctr];
            double v = rule.quad_points_y[qp_ctr];
            
            // our rules are already scaled such that the total weight is 1.0, so normally 
            // multiple by the area to get the quadrature weight for a given triangle.
            // However if we are transforming onto a curved surface, need to use the
            // jacobian (which for a planar triangle is 2.0*area). Factor of two needed
            // because we are using quadrature weights which sum to 1.0 not 0.5 (which
            // would be more technically correct).
            double area_weighting = get_jacobian(u,v) / 2.0;
            Vector qp = get_point(u,v);
            Vector norm = get_normal(u,v);
            tmp_quad_points.push_back(QuadPointType(qp, norm, rule.quad_weights[qp_ctr]*area_weighting));
        }
        
        // commit changes
        quad_points.insert(quad_points.end(), tmp_quad_points.begin(), tmp_quad_points.end());
        
        return;
    }

    template <typename QuadPointType>
    void get_quad_points(std::vector<QuadPointType>& quad_points, const QuadratureRule& rule, unsigned int subdivides) const
    {
        // this function can throw exceptions if the triangle is so slender as to produce
        // NaNs in Quadrature Points.  So don't commit changes to the quad_points until
        // all succeed.
        std::vector<QuadPointType> tmp_quad_points;
        std::vector< QuadPointT<double> > parametric_quad_points;

        double inc = 1.0 / (1 + subdivides);
        for(unsigned int vctr=0; vctr <= subdivides; ++vctr)
        {
            double bary_u1 = 0.0;
            double bary_v1 = vctr*inc;
            double bary_w1 = 1.0 - bary_v1 - bary_u1;

            double bary_u3 = bary_u1;
            double bary_v3 = bary_v1 + inc;
            double bary_w3 = 1.0 - bary_u3 - bary_v3;

            double bary_w2 = bary_w3;
            double bary_v2 = bary_v1;
            double bary_u2 = 1.0 - bary_w2 - bary_v2;

            unsigned int row_count = 2*(subdivides - vctr) + 1;

            while (row_count)
            {
                if (bary_v3 <= bary_v1 || bary_u2 <= bary_u1) { throw std::exception(); }
                
                BasicTriangle parametric_triangle( Vector(bary_u1,bary_v1,bary_w1), 
                                                   Vector(bary_u2,bary_v2,bary_w2),
                                                   Vector(bary_u3,bary_v3,bary_w3) );
                parametric_triangle.get_quad_points(parametric_quad_points, rule);

                row_count--;

                double bary_v4=bary_v3;
                double bary_u4=bary_u2;
                double bary_w4 = 1.0 - bary_v4 - bary_u4;

                if (row_count)
                {
                    row_count--;

                    BasicTriangle parametric_triangle2( Vector(bary_u2,bary_v2,bary_w2), 
                                                        Vector(bary_u4,bary_v4,bary_w4),
                                                        Vector(bary_u3,bary_v3,bary_w3) );
                    parametric_triangle2.get_quad_points(parametric_quad_points, rule);

                    bary_u1 = bary_u2;
                    bary_v1 = bary_v2;
                    bary_w1 = bary_w2;
                    bary_u3 = bary_u4;
                    bary_v3 = bary_v4;
                    bary_w3 = 1.0 - bary_u3 - bary_v3;
                    bary_w2 = bary_w3;
                    bary_v2 = bary_v1;
                    bary_u2 = 1.0 - bary_w2 - bary_v2;
                }

            }

        }
        
        // convert parametric coords into real coords
        for (typename std::vector< QuadPointT<double> >::const_iterator pit=parametric_quad_points.begin(), pend=parametric_quad_points.end();
             pit != pend; ++pit)
        {
            const Vector& qp_bary = pit->pt();
            const double& u = qp_bary.x;
            const double& v = qp_bary.y;
            const double& w = qp_bary.z;
            assert(fabs(u+v+w - 1.0) < 1e-6);

            // renormalise weights to unity (parametric triangle has total area= sqrt(3.0)/2.0)
            double wt = pit->weight() * 2.0 / sqrt(3.0);
            
            // use jacobean to convert from barycentric coordinates to real.
            // note the factor of 2.0 because our quadrature weights add to 1.0
            // rather than 0.5.
            wt *= get_jacobian(u,v) / 2.0;

            Vector qp = get_point(u, v);
            Vector norm = get_normal(u, v);
            tmp_quad_points.push_back(QuadPointType(qp, norm, wt));
        }
        
        // commit changes
        quad_points.insert(quad_points.end(), tmp_quad_points.begin(), tmp_quad_points.end());
        
        return;
    }


    void get_constant_basis_edge_points(const Vector& ref_pt, std::vector<PointNormal>& edge_points) const
    {
        const double centre_uvw = 1./3.;
        const Vector centre(centre_uvw, centre_uvw, centre_uvw);

        // klumsy way to figure out which barycentric quantity to use as the linear
        // weighting function
        Vector origin;
        Vector mid_e1;
        Vector mid_e2;
        if (ref_pt == v1)
        {
            origin = Vector(0,0,1);
            mid_e1 = Vector(0.5,0,0.5);
            mid_e2 = Vector(0,0.5,0.5);
            
        }
        else if (ref_pt == v2)
        {
            origin = Vector(1,0,0);
            mid_e1 = Vector(0.5,0.5,0);
            mid_e2 = Vector(0.5,0,0.5);
            
        }
        else if (ref_pt == v3)
        {
            origin = Vector(0,1,0);
            mid_e1 = Vector(0,0.5,0.5);
            mid_e2 = Vector(0.5,0.5,0);
            
        }
        else
        {
            // shouldn't get here -- ref_pt should be one of the
            // vertices of the original triangle upon which the
            // bezier patch is based
            throw std::exception();
        }
        
        edge_points.push_back(get_point_normal(mid_e1.x, mid_e1.y));
        edge_points.push_back(get_point_normal(centre.x, centre.y));
        edge_points.push_back(get_point_normal(mid_e2.x, mid_e2.y));
        
        return;
    }

    template<typename QuadPointType>
    void get_constant_basis_quad_points(const Vector& ref_pt, 
                                        std::vector<QuadPointType>& quad_points, 
                                        const QuadratureRule& rule,
                                        int subdivides) const
    {
        const double centre_uvw = 1./3.;
        const Vector centre(centre_uvw, centre_uvw, centre_uvw);

        std::vector<QuadPointType> tmp_quad_points;
        
        // klumsy way to figure out which barycentric quantity to use as the linear
        // weighting function
        Vector origin;
        Vector mid_e1;
        Vector mid_e2;
        if (ref_pt == v1)
        {
            origin = Vector(0.0, 0.0, 1.0);
            mid_e1 = Vector(0.5, 0.0, 0.5);
            mid_e2 = Vector(0.0, 0.5, 0.5);
        }
        else if (ref_pt == v2)
        {
            origin = Vector(1.0, 0.0, 0.0);
            mid_e1 = Vector(0.5, 0.5, 0.0);
            mid_e2 = Vector(0.5, 0.0, 0.5);
        }
        else if (ref_pt == v3)
        {
            origin = Vector(0.0, 1.0, 0.0);
            mid_e1 = Vector(0.0, 0.5, 0.5);
            mid_e2 = Vector(0.5, 0.5, 0.0);
        }
        else
        {
            // shouldn't get here -- ref_pt should be one of the
            // vertices of the original triangle
            throw std::exception();
        }
        
        // these are two triangles in parametric space, which make up 
        // one third of the total triangle we are operating on
        BasicTriangle t1(centre, origin, mid_e1);
        BasicTriangle t2(mid_e2, origin, centre);
        
        // get the quadrature points in parametric space -- these are always 
        // handled in double precision, regardless of the precision used by
        // the calling function to hold the final quadrature points
        std::vector< QuadPointT<double> > parametric_quad_points;
        t1.get_quad_points(parametric_quad_points, rule, subdivides);
        t2.get_quad_points(parametric_quad_points, rule, subdivides);
 
        double total=0;
        for (typename std::vector< QuadPointT<double> >::const_iterator pit=parametric_quad_points.begin(), pend=parametric_quad_points.end();
             pit != pend; ++pit)
        {
            total += pit->weight();
        }
        //std::cout << "Total parametric weight (const elements): " << total << std::endl;

        // the parametric quad points are u,v,w coordinates mapping
        // onto points on this triangle.  They already have the correct
        // weighting function with respect to the parametric triangle, 
        // so just need to multiply by the jacobian of the transformation
        // from parameteric coordinates to this triangle.
        for (typename std::vector< QuadPointT<double> >::const_iterator pit=parametric_quad_points.begin(), pend=parametric_quad_points.end();
             pit != pend; ++pit)
        {
            const Vector& qp_bary = pit->pt();
            const double& u = qp_bary.x;
            const double& v = qp_bary.y;
            const double& w = qp_bary.z;
            assert(fabs(u+v+w - 1.0) < 1e-6);
            
            // renormalise weights to unity (parametric triangle is area sqrt(3.0)/2.0
            double wt = pit->weight() * 2.0 / sqrt(3.0); //t1.get_area();
            
            // use jacobean to convert from barycentric coordinates to real.
            // note the factor of 2.0 because our quadrature weights add to 1.0
            // rather than 0.5.
            wt *= get_jacobian(u,v) / 2.0;

            Vector qp = get_point(u,v);
            Vector norm = get_normal(u,v);
            tmp_quad_points.push_back(QuadPointType(qp, norm, wt));
        }

        // commit changes
        quad_points.insert(quad_points.end(), tmp_quad_points.begin(), tmp_quad_points.end());
        
        return;
    }

    template <typename QuadPointType>
    void get_weighted_quad_points(const Vector& ref_pt, std::vector<QuadPointType>& quad_points, const QuadratureRule& rule, unsigned int subdivides=0) const
    {
        assert(false); // TODO: need to fix this function
        
        // temporary store of new quadpoints - just in case we throw exceptions halfway through
        std::vector<QuadPointType> tmp_quad_points;
        
        const double area_weighting = get_area() / ((subdivides+1)*(subdivides+1));

        // klumsy way to figure out which barycentric quantity to use as the linear
        // weighting function
        double u_quadpt;
        double v_quadpt;
        double w_quadpt;
        double* wt;
        if (ref_pt == v1)
        {
            wt = &w_quadpt;
        }
        else if (ref_pt == v2)
        {
            wt = &u_quadpt;
        }
        else if (ref_pt == v3)
        {
            wt = &v_quadpt;
        }
        else
        {
            // shouldn't get here -- ref_pt should be one of the
            // vertices of the original triangle upon which the
            // bezier patch is based
            throw std::exception();
        }

        double inc = 1.0 / (1 + subdivides);
        for(unsigned int vctr=0; vctr <= subdivides; ++vctr)
        {
            double bary_u1 = 0.0;
            double bary_v1 = vctr*inc;
            double bary_w1 = 1.0 - bary_v1 - bary_u1;

            double bary_u3 = bary_u1;
            double bary_v3 = bary_v1 + inc;
            double bary_w3 = 1.0 - bary_u3 - bary_v3;

            double bary_w2 = bary_w3;
            double bary_v2 = bary_v1;
            double bary_u2 = 1.0 - bary_w2 - bary_v2;

            unsigned int row_count = 2*(subdivides - vctr) + 1;

            while (row_count)
            {
                if (bary_v3 <= bary_v1 || bary_u2 <= bary_u1) { throw std::exception(); }

                for (int qp_ctr=0; qp_ctr < rule.num_points; ++qp_ctr)
                {
                    u_quadpt = bary_u1 + (bary_u2-bary_u1)*rule.quad_points_x[qp_ctr] + (bary_u3-bary_u1)*rule.quad_points_y[qp_ctr];
                    v_quadpt = bary_v1 + (bary_v2-bary_v1)*rule.quad_points_x[qp_ctr] + (bary_v3-bary_v1)*rule.quad_points_y[qp_ctr];
                    w_quadpt = 1.0 - u_quadpt - v_quadpt;
                    Vector qp = get_point(u_quadpt, v_quadpt);
                    Vector norm = get_normal(u_quadpt, v_quadpt);
                    tmp_quad_points.push_back(QuadPointType(qp, norm, (*wt)*rule.quad_weights[qp_ctr]*area_weighting, rule.quad_weights[qp_ctr]*area_weighting));
                }

                row_count--;

                double bary_v4=bary_v3;
                double bary_u4=bary_u2;
                double bary_w4 = 1.0 - bary_v4 - bary_u4;

                if (row_count)
                {
                    row_count--;

                    for (int qp_ctr=0; qp_ctr < rule.num_points; ++qp_ctr)
                    {
                        u_quadpt = bary_u2 + (bary_u4-bary_u2)*rule.quad_points_x[qp_ctr] + (bary_u3-bary_u2)*rule.quad_points_y[qp_ctr];
                        v_quadpt = bary_v2 + (bary_v4-bary_v2)*rule.quad_points_x[qp_ctr] + (bary_v3-bary_v2)*rule.quad_points_y[qp_ctr];
                        w_quadpt = 1.0 - u_quadpt - v_quadpt;
                        Vector qp = get_point(u_quadpt, v_quadpt);
                        Vector norm = get_normal(u_quadpt, v_quadpt);
                        tmp_quad_points.push_back(QuadPointType(qp, norm, (*wt)*rule.quad_weights[qp_ctr]*area_weighting, rule.quad_weights[qp_ctr]*area_weighting));
                    }

                    bary_u1 = bary_u2;
                    bary_v1 = bary_v2;
                    bary_w1 = bary_w2;
                    bary_u3 = bary_u4;
                    bary_v3 = bary_v4;
                    bary_w3 = 1.0 - bary_u3 - bary_v3;
                    bary_w2 = bary_w3;
                    bary_v2 = bary_v1;
                    bary_u2 = 1.0 - bary_w2 - bary_v2;
                }

            }

        }
        
        // commit
        quad_points.insert(quad_points.end(), tmp_quad_points.begin(), tmp_quad_points.end());
        
        return;
    }

    static double area_from_three_vectors(const Vector& a, const Vector& b, const Vector& c)
    {
        Vector cross_product = (b-a).cross(c-a);
        return cross_product.length() * 0.5;
    }

    static Vector normal_from_three_vectors(const Vector& a, const Vector& b, const Vector& c)
    {
        Vector cross_product = (b-a).cross(c-a);
        return cross_product.normalised();
    }

    Vector calculate_force(Vector fvals,
                            Vector hvals,
                            double epsilon_int,
                            double epsilon_ext,
                            double kappa,
                            unsigned int subdivides) const
    {
        KahanVector MST_external;
        KahanVector MST_internal;
        KahanVector dbf;
        KahanVector ionic;
        ::calculate_force_components(*this, fvals, hvals, epsilon_int, epsilon_ext, kappa, BasicTriangle::gauss_legendre_pts_7(), subdivides, MST_external, MST_internal, dbf, ionic);
        return *MST_external;
    }
 
    virtual const Vector& get_n1() const { return planar_normal; }
    virtual const Vector& get_n2() const { return planar_normal; }
    virtual const Vector& get_n3() const { return planar_normal; }
  
protected:

    Vector v1;
    Vector v2;
    Vector v3;

private:

    double planar_area;
    Vector planar_centroid;
    Vector planar_normal;

    virtual Vector calc_planar_centroid() const
    {
        return (v1+v2+v3)/3.0;
    }

    Vector calc_planar_normal() const
    {
        Vector product( Vector(v2 - v1).cross(Vector(v3 - v1)) );
        product.normalise();
        return product;
    }

    double calc_planar_area() const
    {
        Vector product( (v2 - v1).cross(v3 - v1) );
        double a = 0.5 * product.length();
        assert(a >= 0.0);
        return a;
    }

};

class Triangle : public BasicTriangle
{
    typedef BasicTriangle super;

protected:

    static const unsigned int MAX_SUBDIVIDES = 5;
    static constexpr double AREA_TOLERANCE = 1e-6;

public:

    friend class Mesh;

    virtual std::string type() const { return "Triangle"; }
    
    Triangle(const PointNormal& a, const PointNormal& b, const PointNormal& c) : 
        BasicTriangle(a.pt(), b.pt(), c.pt()),
        n1(a.n()),
        n2(b.n()),
        n3(c.n())
    {
        n12 = BasicTriangle::get_planar_normal();
        n23 = BasicTriangle::get_planar_normal();
        n31 = BasicTriangle::get_planar_normal();
    }
    
    Triangle(const std::vector<Vertex>& vertex_list,
            const unsigned int _v1_idx,
            const unsigned int _v2_idx,
            const unsigned int _v3_idx)
            :   BasicTriangle(vertex_list[_v1_idx], vertex_list[_v2_idx], vertex_list[_v3_idx]),
                v1_idx(_v1_idx),
                v2_idx(_v2_idx),
                v3_idx(_v3_idx)
    {
        n1 = vertex_list[_v1_idx].normal;
        n2 = vertex_list[_v2_idx].normal;
        n3 = vertex_list[_v3_idx].normal;

        n12 = Vector(0,0,0);
        n23 = Vector(0,0,0);
        n31 = Vector(0,0,0);
    }

    virtual ~Triangle() {}

    Triangle(const Triangle& other) : BasicTriangle(other) {
        v1_idx = other.v1_idx;
        v2_idx = other.v2_idx;
        v3_idx = other.v3_idx;
        n1 = other.n1;
        n2 = other.n2;
        n3 = other.n3;
        n12 = other.n12;
        n23 = other.n23;
        n31 = other.n31;
    }

protected:

    unsigned int v1_idx;
    unsigned int v2_idx;
    unsigned int v3_idx;
    Vector n1;
    Vector n2;
    Vector n3;

    // normal vectors of adjacent triangles -- used for
    // curved meshes in order to get G1 continuity over
    // the mesh (see PNG1.h).
    Vector n12;
    Vector n23;
    Vector n31;

public:

    const Vector& get_n1() const { return n1; }
    const Vector& get_n2() const { return n2; }
    const Vector& get_n3() const { return n3; }

    inline void reinit_vertex_normals(const std::vector<Vertex>& vertex_list) {
        n1 = vertex_list[v1_idx].get_normal();
        n2 = vertex_list[v2_idx].get_normal();
        n3 = vertex_list[v3_idx].get_normal();
        return;
    }

    void set_adjacent_normal(const Triangle& other)
    {
        const unsigned int other_v1_idx = other.get_v1_idx();
        const unsigned int other_v2_idx = other.get_v2_idx();
        const unsigned int other_v3_idx = other.get_v3_idx();

        if (v1_idx != other_v1_idx && v1_idx != other_v2_idx && v1_idx != other_v3_idx)
        {
            n23 = other.get_planar_normal();
        }
        else if (v2_idx != other_v1_idx && v2_idx != other_v2_idx && v2_idx != other_v3_idx)
        {
            n31 = other.get_planar_normal();
        }
        else if (v3_idx != other_v1_idx && v3_idx != other_v2_idx && v3_idx != other_v3_idx)
        {
            n12 = other.get_planar_normal();
        }
    }

    inline unsigned int get_v1_idx() const { return v1_idx; }
    inline unsigned int get_v2_idx() const { return v2_idx; }
    inline unsigned int get_v3_idx() const { return v3_idx; }

    // function to return the other two vertices of the triangle in anticlockwise
    // order from the one passed in as an argument.
    // This is used when constructing NodePatch objects.
    void get_anti_clockwise_vertices_from(const Vector& v_123, Vector &a, Vector &b) const
    {
        if (vectors_are_equal(v_123, v1))
        {
            a = v2;
            b = v3;
        }
        else if(vectors_are_equal(v_123, v2))
        {
            a = v3;
            b = v1;
        }
        else if(vectors_are_equal(v_123, v3))
        {
            a = v1;
            b = v2;
        }
        else
        {
            throw Vertex::BadVertex();
        }
        return;

    }

    void get_anti_clockwise_vertices_from(const Vector& v_123, unsigned int &a, unsigned int &b) const
    {
        if (vectors_are_equal(v_123, v1))
        {
            a = v2_idx;
            b = v3_idx;
        }
        else if(vectors_are_equal(v_123, v2))
        {
            a = v3_idx;
            b = v1_idx;
        }
        else if(vectors_are_equal(v_123, v3))
        {
            a = v1_idx;
            b = v2_idx;
        }
        else
        {
            throw Vertex::BadVertex();
        }
        return;

    }

    Vector circumcentre() const
    {

        // circumcentre
        double denom = 2.0f*(super::v1-super::v2).cross(super::v2-super::v3).length2();
        double alpha = (super::v2-super::v3).length2() * ((super::v1-super::v2).dot(super::v1-super::v3)) / denom;
        double beta =  (super::v1-super::v3).length2() * ((super::v2-super::v1).dot(super::v2-super::v3)) / denom;
        double gamma = (super::v1-super::v2).length2() * ((super::v3-super::v1).dot(super::v3-super::v2)) / denom;

        return Vector(super::v1*alpha + super::v2*beta + super::v3*gamma);

    }

    Vector curveycentre() const;

    void get_vertex_indices(unsigned int &a, unsigned int &b, unsigned int &c) const
    {
        a = v1_idx;
        b = v2_idx;
        c = v3_idx;
        return;
    }

    friend std::ostream& operator<< (std::ostream &os, BasicTriangle t);
 
};

// stream operator for Triangle class
inline std::ostream& operator<<(std::ostream& os, BasicTriangle t)
{
    os << "[" << t.get_v1() << "," << t.get_v2() << "," << t.get_v3() << "](a=" << t.get_planar_area() << ")";
    return os;
}

#endif /* TRIANGLE_H_ */
