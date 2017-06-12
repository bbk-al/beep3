/*
* png1.h
*
*  Created on: 12 Jan 2010
*      Author: david
*/

#ifndef PNG1_H_
#define PNG1_H_

#include "triangle.h"
#include "../common/math_vector.h"

class PNG1_Triangle : public Triangle
{

public:

    virtual std::string type() const { return "PNG1_Triangle"; }
    
    // complicated function! transforms the third vertex of each triangle into
    // the plane defined by the other triangle, such that two new triangles can
    // be formed (tr1.v1 -> tr1.v2 -> tr2_v3_in_tr1 = affine transformation of
    // tr2 into tr1) and (tr2.v1 -> tr2.v2 -> tr1_v3_in_tr2 = affine transformation
    // of tr1 into tr2.
    static void affine_transformation_between_planes(const BasicTriangle& tr1,
                                                     const BasicTriangle& tr2,
                                                     Vector& tr1_v3_in_tr2,
                                                     Vector& tr2_v3_in_tr1)
    {
        // this could be simplified

        // get location of the third vertex of tr1
        const Vector& origin_plane_1 = tr1.get_v1();
        Vector x_axis_p1 = (tr1.get_v2() - origin_plane_1);
        double scale1 = x_axis_p1.length();
        x_axis_p1.normalise();
        Vector y_axis_p1 = tr1.get_planar_normal().cross(x_axis_p1).normalised();
        double dx1 = (tr1.get_v3() - origin_plane_1).dot(x_axis_p1);
        double dy1 = (tr1.get_v3() - origin_plane_1).dot(y_axis_p1);

        // re-express tr1.v3() in the coordinate frame of tr2, scaling the
        // coordinates such that the base lengths match
        const Vector& origin_plane_2 = tr2.get_v1();
        Vector x_axis_p2 = (tr2.get_v2() - origin_plane_2);
        double scale2 = x_axis_p2.length();
        x_axis_p2.normalise();
        Vector y_axis_p2 = tr2.get_planar_normal().cross(x_axis_p2).normalised();
        double dx2 = (tr2.get_v3() - origin_plane_2).dot(x_axis_p2);
        double dy2 = (tr2.get_v3() - origin_plane_2).dot(y_axis_p2);

        tr1_v3_in_tr2 = origin_plane_2 + (x_axis_p2*dx1 + y_axis_p2*dy1)*scale2/scale1;
        tr2_v3_in_tr1 = origin_plane_1 + (x_axis_p1*dx2 + y_axis_p1*dy2)*scale1/scale2;

        return;

    }

    static Vector affine_transformation_between_planes(const BasicTriangle& tr1,
                                                       const Vector& a,
                                                       const Vector& b,
                                                       const Vector& norm)
    {
        // get location of the third vertex of tr1
        Vector x_axis_p1 = (tr1.get_v2() - tr1.get_v1());
        double scale1 = x_axis_p1.length();
        x_axis_p1.normalise();
        Vector y_axis_p1 = tr1.get_planar_normal().cross(x_axis_p1).normalised();
        Vector v3rel = tr1.get_v3() - tr1.get_v1();
        double dx = (v3rel).dot(x_axis_p1);
        double dy = (v3rel).dot(y_axis_p1);

        // re-express tr1.v3() in the coordinate frame of tr2, scaling the
        // coordinates such that the base lengths match
        Vector x_axis_p2 = b - a;
        double scale2 = x_axis_p2.length();
        x_axis_p2.normalise();
        Vector y_axis_p2 = norm.cross(x_axis_p2).normalised();

        return a + (x_axis_p2*dx + y_axis_p2*dy)*scale2/scale1;

    }

    virtual ~PNG1_Triangle() {}

    PNG1_Triangle(const Triangle& planar_triangle)
        : Triangle(planar_triangle)
    {
        // init the control points

        // corners
        b300 = v1;
        b030 = v2;
        b003 = v3;

        // control points immediately adjacent to corners
        b_p1_210 = project_vector_to_plane_and_rescale(v2, v1,n1);
        b_p1_201 = project_vector_to_plane_and_rescale(v3, v1,n1);
        b_p2_120 = project_vector_to_plane_and_rescale(v1, v2,n2);
        b_p2_021 = project_vector_to_plane_and_rescale(v3, v2,n2);
        b_p3_102 = project_vector_to_plane_and_rescale(v1, v3,n3);
        b_p3_012 = project_vector_to_plane_and_rescale(v2, v3,n3);

        // project p201 in tangent plane of p1 into tangent plane of p2 --> find b_p1_021
        affine_transformation_between_planes(BasicTriangle(b300, b_p1_210, b_p1_201), BasicTriangle(b_p2_120, b030, b_p2_021), b_p1_021, b_p2_201);
        affine_transformation_between_planes(BasicTriangle(b030, b_p2_021, b_p2_120), BasicTriangle(b_p3_012, b003, b_p3_102), b_p2_102, b_p3_120);
        affine_transformation_between_planes(BasicTriangle(b003, b_p3_102, b_p3_012), BasicTriangle(b_p1_201, b300, b_p1_210), b_p3_210, b_p1_012);

        const Vector& Nt = get_planar_normal(); // call to (grand)parent class function
        Vector average_norm_12 = (Nt + n12).normalised(); // n12 is adjacent normal across side 1-2
        Vector average_norm_23 = (Nt + n23).normalised(); // n23 is adjacent normal across side 2-3
        Vector average_norm_31 = (Nt + n31).normalised(); // n31 is adjacent normal across side 3-1
        b_12_p1_111 = affine_transformation_between_planes(BasicTriangle(b300, b_p1_210, b_p1_201), b_p1_210, b_p2_120, average_norm_12);
        b_12_p2_111 = affine_transformation_between_planes(BasicTriangle(b_p2_120, b030, b_p2_021), b_p1_210, b_p2_120, average_norm_12);
        b_23_p2_111 = affine_transformation_between_planes(BasicTriangle(b030, b_p2_021, b_p2_120), b_p2_021, b_p3_012, average_norm_23);
        b_23_p3_111 = affine_transformation_between_planes(BasicTriangle(b_p3_012, b003, b_p3_102), b_p2_021, b_p3_012, average_norm_23);
        b_31_p3_111 = affine_transformation_between_planes(BasicTriangle(b003, b_p3_102, b_p3_012), b_p3_102, b_p1_201, average_norm_31);
        b_31_p1_111 = affine_transformation_between_planes(BasicTriangle(b_p1_201, b300, b_p1_210), b_p3_102, b_p1_201, average_norm_31);

        // find the area of the curved patch by recursive subdivision
        set_area_and_normal();

    }

    static Vector project_vector_to_plane_and_rescale(const Vector& pt, const Vector& point_in_plane, const Vector& normal)
    {
        Vector projected_point = pt + normal * normal.dot(point_in_plane - pt);
        Vector relative_direction = (projected_point - point_in_plane).normalised();
        return point_in_plane + relative_direction*((point_in_plane - pt).length() / 3.0);
    }

    virtual double get_area() const {
        return curved_area;
    }

    // returns the point with parametric coordinates (u,v,w) (u+v+w = 1.0)
    virtual Vector get_point(double u, double v) const 
    {

        const double w = 1.0 - u - v;
        Vector S =  b300*w*w*w + b030*u*u*u + b003*v*v*v
                + b210(u,v,w)*3*w*w*u + b120(u,v,w)*3*w*u*u + b201(u,v,w)*3*w*w*v
                + b021(u,v,w)*3*u*u*v + b102(u,v,w)*3*w*v*v + b012(u,v,w)*3*u*v*v
                + b111(u,v,w)*6*w*u*v;
        return S;
    }

    virtual double get_jacobian(double u, double v) const 
    {
        // return the jacobean of the projection from planar to whatever
        // type of triangle this is.  For planar, should return 2.0*area.
        // For curved elements, will depend on local curvature.
        double w = 1.0 - u - v;
        Vector du = get_dU(u,v,w);
        Vector dv = get_dV(u,v,w);
        double E = du.length2();
        double F = du.dot(dv);
        double G = dv.length2();
        double jac = sqrt(E*G - F*F);
        return jac;
    }

    // returns the derivative w.r.t parametric coordinate u @ (u,v,w)
    virtual Vector get_dU(double u, double v, double w) const
    {
        Vector dU = b030*3*u*u - b300*3*w*w + b120(u,v,w)*(6*u*w - 3*u*u)
                    + b021(u,v,w)*6*u*v - b102(u,v,w)*3*v*v + b012(u,v,w)*3*v*v
                    + b210(u,v,w)*(3*w*w - 6*w*u) - b201(u,v,w)*6*w*v
                    + b111(u,v,w)*(6*w*v - 6*u*v);
        return dU;
    }

    // returns the derivative w.r.t parametric coordinate v @ (u,v,w)
    virtual Vector get_dV(double u, double v, double w) const
    {
        Vector dV = b003*3*v*v - b300*3*w*w - b120(u,v,w)*3*u*u
                    + b021(u,v,w)*3*u*u + b102(u,v,w)*(6*v*w - 3*v*v) + b012(u,v,w)*6*v*u
                    - b210(u,v,w)*6*w*u + b201(u,v,w)*(3*w*w - 6*w*v)
                    + b111(u,v,w)*(6*w*u - 6*u*v);
        return dV;
    }

    virtual Vector get_normal(double u, double v) const
    {
//         static const double eps=1e-6;
//         Vector pt = get_point(u,v);
//         Vector dU = get_point(u+eps,v) - pt;
//         Vector dV = get_point(u,v+eps) - pt;
//         Vector n_numerical = dU.cross(dV).normalised();
//         return n_numerical;
        double w = 1.0 - u - v;
        Vector dBu = b030*3*u*u - b300*3*w*w + b120(u,v,w)*(6*u*w - 3*u*u)
                    + b021(u,v,w)*6*u*v - b102(u,v,w)*3*v*v + b012(u,v,w)*3*v*v
                    + b210(u,v,w)*(3*w*w - 6*w*u) - b201(u,v,w)*6*w*v
                    + b111(u,v,w)*(6*w*v - 6*u*v);

        Vector dBv = b003*3*v*v - b300*3*w*w - b120(u,v,w)*3*u*u
                    + b021(u,v,w)*3*u*u + b102(u,v,w)*(6*v*w - 3*v*v) + b012(u,v,w)*6*v*u
                    - b210(u,v,w)*6*w*u + b201(u,v,w)*(3*w*w - 6*w*v)
                    + b111(u,v,w)*(6*w*u - 6*u*v);

        return dBu.cross(dBv).normalised();
    }

private:

    // Blending functions
    inline Vector b201(double u, double v, double w) const {
        double x = 1.0 - u;
        double xx = x*x;
        double uu = u*u;
        return (xx*b_p1_201 + uu*b_p2_201) / (uu+xx);
    }
    inline Vector b102(double u, double v, double w) const {
        double x = 1.0 - u;
        double xx = x*x;
        double uu = u*u;
        return (xx*b_p3_102 + uu*b_p2_102) / (uu+xx);
    }
    inline Vector b012(double u, double v, double w) const {
        double x = 1.0 - w;
        double xx = x*x;
        double ww = w*w;
        return (xx*b_p3_012 + ww*b_p1_012) / (ww+xx);
    }
    inline Vector b021(double u, double v, double w) const {
        double x = 1.0 - w;
        double xx = x*x;
        double ww = w*w;
        return (xx*b_p2_021 + ww*b_p1_021) / (ww+xx);
    }
    inline Vector b120(double u, double v, double w) const {
        double x = 1.0 - v;
        double xx = x*x;
        double vv = v*v;
        return (xx*b_p2_120 + vv*b_p3_120) / (vv+xx);
    }
    inline Vector b210(double u, double v, double w) const {
        double x = 1.0 - v;
        double xx = x*x;
        double vv = v*v;
        return (xx*b_p1_210 + vv*b_p3_210) / (vv+xx);
    }

    inline Vector b111(double u, double v, double w) const
    {
        static const double eps=1e-6;
        // need this for stability - otherwise the corner values
        // (i.e. u=1 or v=1 or w=1) will lead to 0/0 ratios
        u += eps;
        v += eps;
        w += eps;
        double tot = u+v+w;
        u /= tot;
        v /= tot;
        w /= tot;

        double pre1 = u*w / (u*v + u*w + v*w);
        Vector term1 = (1.-u)*w*b_12_p1_111 + (1.-w)*u*b_12_p2_111;
        double denom1 = w+u-(2.*u*w);

        double pre2 = u*v / (u*v + u*w + v*w);
        Vector term2 = (1.-v)*u*b_23_p2_111 + (1.-u)*v*b_23_p3_111;
        double denom2 = u+v-(2.*u*v);

        double pre3 = v*w / (u*v + u*w + v*w);
        Vector term3 = (1.-w)*v*b_31_p3_111 + (1.-v)*w*b_31_p1_111;
        double denom3 = w+v-(2.*v*w);

        return term1*(pre1/denom1) + term2*(pre2/denom2) + term3*(pre3/denom3);
    }

    void set_area_and_normal()
    {
        const double planar_area = get_planar_area();
        const double tol = AREA_TOLERANCE * planar_area;

        Vector weighted_area_normal(0,0,0);
        double weighted_area=0;

        // integrate over u and v parametric coords
        unsigned int n=1;
        double a1=0;
        subdivided_area(n, a1, weighted_area, mean_normal, weighted_area_normal);
        double a2=0;
        while (true)
        {
            subdivided_area(++n, a2, weighted_area, mean_normal, weighted_area_normal);
            if ( (a2 - a1) <= tol ) { break; }

            if (n > MAX_SUBDIVIDES)
            { break; }
            a1 = a2;
        }

        curved_area = a2;
        area_times_normal = mean_normal;
        mean_normal.normalise();
        num_subdivides = n;

        return;

    }

    void subdivided_area(unsigned int num_subdivisions,
                        double& area,
                        double& weighted_area,
                        Vector& normal,
                        Vector& weighted_area_normal) const
    {

        area = 0;
        weighted_area = 0;
        normal = Vector(0,0,0);
        weighted_area_normal = Vector(0,0,0);

        double inc = 1.0 / (1 + num_subdivisions);

        for(unsigned int vctr=0; vctr <= num_subdivisions; ++vctr)
        {
            double u1 = 0.0;
            double v1 = vctr*inc;
            double w1 = 1.0 - v1 - u1;

            double u3 = u1;
            double v3 = v1 + inc;
            double w3 = 1.0 - u3 - v3;

            double w2 = w3;
            double v2 = v1;
            double u2 = 1.0 - w2 - v2;

            unsigned int row_count = 2*(num_subdivisions - vctr) + 1;

            while (row_count)
            {
                if (v3 <= v1 || u2 <= u1) { throw std::exception(); }

                // create triangle
                Vector a = get_point(u1,v1);
                Vector b = get_point(u2,v2);
                Vector c = get_point(u3,v3);

                double wmean = (w1 + w2 + w3) / 3.0;
                double area_tmp = area_from_three_vectors(a,b,c);
                area += area_tmp;
                weighted_area += area_tmp*wmean;
                normal += area_tmp * normal_from_three_vectors(a,b,c);
                weighted_area_normal += area_tmp*wmean*normal_from_three_vectors(a,b,c);
                row_count--;

                double v4=v3;
                double u4=u2;
                double w4 = 1.0 - v4 - u4;

                if (row_count)
                {
                    Vector d = get_point(u4,v4);
                    wmean = (w2+w3+w4) / 3.0;
                    area_tmp = area_from_three_vectors(b,d,c);
                    area += area_tmp;
                    weighted_area += area_tmp*wmean;
                    normal += area_tmp * normal_from_three_vectors(b,d,c);
                    weighted_area_normal += area_tmp*wmean*normal_from_three_vectors(b,d,c);
                    row_count--;

                    u1 = u2;
                    v1 = v2;
                    w1 = w2;
                    u3 = u4;
                    v3 = v4;
                    w3 = 1.0 - u3 - v3;
                    w2 = w3;
                    v2 = v1;
                    u2 = 1.0 - w2 - v2;
                }

            }

        }

        return;
    }

    unsigned int num_subdivides;
    double curved_area;
    Vector area_times_normal;
    Vector mean_normal;

    Vector b300;
    Vector b030;
    Vector b003;

    Vector b_p1_201;
    Vector b_p2_201;
    Vector b_p3_102;
    Vector b_p2_102;
    Vector b_p3_012;
    Vector b_p1_012;
    Vector b_p2_021;
    Vector b_p1_021;
    Vector b_p2_120;
    Vector b_p3_120;
    Vector b_p1_210;
    Vector b_p3_210;

    Vector b_12_p1_111;
    Vector b_12_p2_111;
    Vector b_23_p2_111;
    Vector b_23_p3_111;
    Vector b_31_p1_111;
    Vector b_31_p3_111;

};

#endif /* PNG1_H_ */
