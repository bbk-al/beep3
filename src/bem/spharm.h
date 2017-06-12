/*
 * spharm.h
 *
 *  Created on: 21 Jul 2010
 *      Author: david
 */

#ifndef SPHARM_H_
#define SPHARM_H_

#include "../common/math_vector.h"
#include <boost/shared_array.hpp>
#include <memory.h>
#include <vector>

class SpharmPt
{
public:

    SpharmPt(const Vector& r, const Vector& centre, double fval=0) : f(fval)
    {
        double rho, theta;

        // NB: I have reversed the usual (maths) convention of phi and theta, since
        // spherical harmnonics more usually follow the physics convention.  So, to clarify:
        // theta is the polar angle, ranging from 0 at 'north pole' to pi at 'south pole'
        // phi is the rotational angle (which, confusingly, would be denoted theta in
        // simple polar coordinates) and ranges from 0 to 2pi.
        (r - centre).toSpherical(rho, phi, theta); // phi and theta reversed from usual calling convention

        cos_theta = cos(theta);
        sin_theta = sin(theta);
    }

    // copy ctor
    SpharmPt(const SpharmPt& other) : cos_theta(other.cos_theta), sin_theta(other.sin_theta), phi(other.phi), f(other.f) {}

    double cos_theta;
    double sin_theta;
    double phi;
    double f;
};

typedef std::vector<SpharmPt> SpharmList;

class SpharmHolder
{
public:

    SpharmHolder(unsigned int n) : max_n(n)
    {
        __size = (n+1)*(n+1);
        data = boost::shared_array<double>(new double[__size]);
    }

    SpharmHolder(const SpharmHolder& other) : __size(other.__size),  max_n(other.max_n)
    {
        data = boost::shared_array<double>(new double[__size]);
        memcpy(data.get(), other.data.get(), __size * sizeof(double));
    }

    inline double& operator() (int n, int m)
    {
        assert(n >= 0);
        assert(n <= static_cast<int>(max_n));
        assert(fabs(m) <= n);
        return data[n*n + m+n];
    }

    inline const double& operator() (int n, int m) const
    {
        assert(n >= 0);
        assert(n <= static_cast<int>(max_n));
        assert(fabs(m) <= n);
        return data[n*n + m+n];
    }

    inline double& operator() (unsigned int i)
    {
        assert(i < __size);
        return data[i];
    }

    inline const double& operator() (unsigned int i) const
    {
        assert(i < __size);
        return data[i];
    }

    inline unsigned int size() const { return __size; }
    inline unsigned int max_degree() const { return max_n; }

    double evaluate(const Vector& pt, const Vector& origin) const
    {
        SpharmPt spt(pt, origin);
        double val = evaluate_spt(spt);
        //double dtheta = eval_dtheta(spt);
        //double dphi = eval_dphi(spt);

        return val;
    }

    Vector dArea(const Vector& pt, const Vector& origin) const;
    Vector dphi(const Vector& pt, const Vector& origin) const;
    Vector dtheta(const Vector& pt, const Vector& origin) const;

    double evaluate_spt(const SpharmPt& sp) const
    {
        double sp_val = 0.0;

        for (int nctr=0; nctr <= static_cast<int>(max_n); ++nctr)
        {
            for (int mctr=-nctr; mctr <= nctr; ++mctr)
            {
                double yval = SpharmHolder::Ynm(nctr, mctr, sp.cos_theta, sp.phi);
                //std::cout << "(n,m)=(" << nctr << "," << mctr << ") = " << (*this)(nctr,mctr) << std::endl;
                sp_val += (*this)(nctr,mctr) * yval;
            }
        }
        return sp_val;
    }

    static double Ynm(int n, int m, double cos_theta, double phi);
    static double Ynm_dphi(int n, int m, double cos_theta, double phi);
    static void Ynm_dtheta(int nmax, double cos_theta, double sin_theta, double phi, SpharmHolder& holder);

private:

    unsigned int __size;
    unsigned int max_n;
    boost::shared_array<double> data;
};

SpharmHolder least_squares_spherical_harmonic(const SpharmList& sample_pts, unsigned int n);
SpharmHolder least_squares_spherical_harmonic(const std::vector<Vector>& sample_pts, const Vector& origin, const std::vector<double>& vals, unsigned int n);



#endif /* SPHARM_H_ */
