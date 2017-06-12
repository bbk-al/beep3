/*
* bezier.h
*
*  Created on: 5 Aug 2010
*      Author: david
*/

#ifndef BEZIER_H_
#define BEZIER_H_

#include "triangle.h"
#include "../common/math_vector.h"
#include "quad_point.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <exception>

#define FORCE_PLANAR false

// interpolates a point normal triangle into a curved triangular bezier patch
class HybridBezierTriangle : public Triangle
{

public:

    HybridBezierTriangle(const PointNormal& a, const PointNormal& b, const PointNormal& c) : Triangle(a,b,c)
    {
        init();
    }
    
    HybridBezierTriangle(const Triangle& planar_triangle) : Triangle(planar_triangle)
    {
        init();
    }

    virtual std::string type() const { return "HybridBezierTriangle"; } 

    inline void init()
    {
        // failsafe for dodgy triangles which do strange things like wraparound
        if (FORCE_PLANAR || n1.dot(n2) < 0 || n2.dot(n3) < 0 || n3.dot(n1) < 0)
        {

            std::cerr << "Resetting normal vectors to v1-v2-v3 anticlockwise definition because you passed some weird normals..." << std::endl;
            Vector n = normal_from_three_vectors(v1,v2,v3);
            n1 = n;
            n2 = n;
            n3 = n;
        }

        // cubic bezier triangle geometry fitting
        b300 = v1;
        b030 = v2;
        b003 = v3;
        double w12 = (v2 - v1).dot(n1);
        double w21 = (v1 - v2).dot(n2);
        double w13 = (v3 - v1).dot(n1);
        double w31 = (v1 - v3).dot(n3);
        double w23 = (v3 - v2).dot(n2);
        double w32 = (v2 - v3).dot(n3);
        b210 = (v1*2 + v2 - n1*w12)/3.0;
        b120 = (v2*2 + v1 - n2*w21)/3.0;
        b021 = (v2*2 + v3 - n2*w23)/3.0;
        b012 = (v3*2 + v2 - n3*w32)/3.0;
        b102 = (v3*2 + v1 - n3*w31)/3.0;
        b201 = (v1*2 + v3 - n1*w13)/3.0;
        E = (b210 + b120 + b021 + b012 + b102 + b201) / 6.0;
        V = (v1 + v2 + v3) / 3.0;
        b111 = E + (E-V)/2.0;

        // calculate area (assume 5 subdivisions is enough)
        curved_area = subdivided_area(15);
    }

    // derived classes need virtual destructors
    virtual ~HybridBezierTriangle() {}

    virtual double get_area() const {
        return curved_area;
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
        return sqrt(E*G - F*F);
    }

    virtual Vector get_point(double u, double v) const
    {

        // get position on patch at u,v,w coordinates
        double w = 1.0 - u - v;

        Vector b =  b300*w*w*w   + b030*u*u*u   + b003*v*v*v
                  + b210*3*w*w*u + b120*3*w*u*u + b201*3*w*w*v
                  + b021*3*u*u*v + b102*3*w*v*v + b012*3*u*v*v
                  + b111*6*w*u*v;

        //std::cout << "HybridBezierTriangle::get_point @ ( " << u << " " << v << ") = " << b << "\n";
        return b;
    }

    virtual Vector get_dU(double u, double v, double w) const
    {
        //std::cout << "HybridBezierTriangle::get_dU" << std::endl;
        Vector dU = b030*3*u*u - b300*3*w*w + b120*(6*u*w - 3*u*u)
                    + b021*6*u*v - b102*3*v*v + b012*3*v*v
                    + b210*(3*w*w - 6*w*u) - b201*6*w*v
                    + b111*(6*w*v - 6*u*v);
        return dU;
    }

    virtual Vector get_dV(double u, double v, double w) const
    {
        //std::cout << "HybridBezierTriangle::get_dV" << std::endl;
        Vector dV = b003*3*v*v - b300*3*w*w - b120*3*u*u
                    + b021*3*u*u + b102*(6*v*w - 3*v*v) + b012*6*v*u
                    - b210*6*w*u + b201*(3*w*w - 6*w*v)
                    + b111*(6*w*u - 6*u*v);
        return dV;
    }

    virtual Vector get_normal(double u, double v) const
    {
        double w = 1.0 - u - v;
        return get_dU(u,v,w).cross(get_dV(u,v,w)).normalised();
        
    }

private:

    double subdivided_area(unsigned int num_subdivisions) const
    {

        double ar = 0;

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

                ar += area_from_three_vectors(a,b,c);
                row_count--;

                double v4=v3;
                double u4=u2;

                if (row_count)
                {
                    Vector d = get_point(u4,v4);
                    ar += area_from_three_vectors(b,d,c);
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

        return ar;
    }

    // geometry control points
    Vector b300;
    Vector b030;
    Vector b003;
    Vector b210;
    Vector b120;
    Vector b021;
    Vector b012;
    Vector b102;
    Vector b201;
    Vector E;
    Vector V;
    Vector b111;

    double curved_area;

};

#endif /* BEZIER_H_ */
