/*
 * spharm.cpp
 *
 *  Created on: 21 Jul 2010
 *      Author: david
 */

#ifdef PYTHON_MODULE
#include <boost/python.hpp>
using namespace boost::python;
#endif

#include "spharm.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

const double root2 = 1.4142135623730951;

double SpharmHolder::Ynm(int n, int m, double cos_theta, double phi)
{
    assert(n >= 0);
    assert(fabs(m) <= n);
    if (m >= 0)
    {
        return root2 * gsl_sf_legendre_sphPlm(n, m, cos_theta) * cos(phi*fabs(m));
    }
    else
    {
        return root2 * gsl_sf_legendre_sphPlm(n, -m, cos_theta) * sin(phi*fabs(m));
    }
}

double SpharmHolder::Ynm_dphi(int n, int m, double cos_theta, double phi)
{
    assert(n >= 0);
    assert(fabs(m) <= n);
    if (m >= 0)
    {
        return root2 * gsl_sf_legendre_sphPlm(n, m, cos_theta) * (fabs(m)*-sin(phi*fabs(m)));
    }
    else
    {
        return root2 * gsl_sf_legendre_sphPlm(n, -m, cos_theta) * (fabs(m)*cos(phi*fabs(m)));
    }
}

void SpharmHolder::Ynm_dtheta(int nmax, double cos_theta, double sin_theta, double phi, SpharmHolder& holder)
{
    assert(nmax >= 0);
    double *result_array = new double[(nmax+1)*(nmax+1)];
    double *result_deriv_array = new double[(nmax+1)*(nmax+1)];

    for (int m=-nmax; m <= nmax; ++m)
    {
        if (m >= 0)
        {
            gsl_sf_legendre_sphPlm_deriv_array(nmax, m, cos_theta, result_array, result_deriv_array);
            int res_ctr=0;
            for (int nctr=fabs(m); nctr <= nmax; ++nctr)
            {
                holder(nctr, m) = result_deriv_array[res_ctr++] * root2 * cos(phi*fabs(m)) * -sin_theta;
            }
        }
        else
        {
            gsl_sf_legendre_sphPlm_deriv_array(nmax, -m, cos_theta, result_array, result_deriv_array);
            int res_ctr=0;
            for (int nctr=fabs(m); nctr <= nmax; ++nctr)
            {
                holder(nctr, m) = result_deriv_array[res_ctr++] * root2 * sin(phi*fabs(m)) * -sin_theta;
            }
        }
    }

    delete[] result_array;
    delete[] result_deriv_array;

    return;
}

Vector SpharmHolder::dphi(const Vector& pt, const Vector& origin) const
{
    SpharmPt sp(pt, origin);
    double sp_fval = 0.0;

    for (int nctr=0; nctr <= static_cast<int>(max_n); ++nctr)
    {
        for (int mctr=-nctr; mctr <= nctr; ++mctr)
        {
            double yval = SpharmHolder::Ynm_dphi(nctr, mctr, sp.cos_theta, sp.phi);
            sp_fval += (*this)(nctr,mctr) * yval;
        }
    }

    Vector dphi_dir(-sp.sin_theta*sin(sp.phi),
                     sp.sin_theta*cos(sp.phi),
                     0);

    if (sp.sin_theta == 0.0)
    {
        return Vector(0,0,0);
    }

    return dphi_dir.normalised() * sp_fval;

}

Vector SpharmHolder::dArea(const Vector& pt, const Vector& origin) const
{
    SpharmPt sp(pt, origin);
    Vector dtheta(sp.cos_theta * cos(sp.phi),
                  sp.cos_theta * sin(sp.phi),
                  -sp.sin_theta);
    Vector dphi(-sp.sin_theta*sin(sp.phi),
                 sp.sin_theta*cos(sp.phi),
                 0);
    Vector dS = dtheta.cross(dphi);
    return dS;
}

Vector SpharmHolder::dtheta(const Vector& pt, const Vector& origin) const
{
    SpharmPt sp(pt, origin);
    double val = 0.0;
    std::auto_ptr<SpharmHolder> tmp_ptr(new SpharmHolder(max_n));
    SpharmHolder& tmp = *tmp_ptr;
    Ynm_dtheta(static_cast<int>(max_n), sp.cos_theta, sp.sin_theta, sp.phi, tmp);

    for (int nctr=0; nctr <= static_cast<int>(max_n); ++nctr)
    {
        for (int mctr=-nctr; mctr <= nctr; ++mctr)
        {
            val += (*this)(nctr,mctr) * tmp(nctr, mctr);
        }
    }

    Vector dtheta_dir(sp.cos_theta * cos(sp.phi),
                      sp.cos_theta * sin(sp.phi),
                      -sp.sin_theta);

    return dtheta_dir.normalised() * val;
}

SpharmHolder least_squares_spherical_harmonic(const std::vector<Vector>& sample_pts, const Vector& origin, const std::vector<double>& vals, unsigned int n)
{
    SpharmList sph_pts;
    unsigned int ctr=0;
    for (std::vector<Vector>::const_iterator it=sample_pts.begin(), end=sample_pts.end();
         it != end;
         ++it)
    {
        sph_pts.push_back(SpharmPt(*it, origin, vals[ctr++]));
    }

    return least_squares_spherical_harmonic(sph_pts, n);
}

SpharmHolder least_squares_spherical_harmonic(const SpharmList& sample_pts, unsigned int n)
{
    // n is max order of spherical harmonic to combine to approximate the function at the sample_pts
    SpharmHolder spharm(n);

    // Allocate memory for matrix
    double *a_data = new double[spharm.size() * sample_pts.size()];
    double *f_data = new double[sample_pts.size()];
    gsl_matrix_view m = gsl_matrix_view_array(a_data, sample_pts.size(), spharm.size());
    gsl_vector_view f = gsl_vector_view_array(f_data, sample_pts.size());

    // for each spherical harmonic point put an element in the matrix
    int row_ctr=0;
    for (SpharmList::const_iterator it=sample_pts.begin(), end=sample_pts.end(); it != end; ++it)
    {
        const SpharmPt& sp = *it;

        // set vector element in this row
        gsl_vector_set(&f.vector, row_ctr, sp.f);

        // set matrix elements in this row
        int col_ctr=0;
        for (int nctr=0; nctr <= static_cast<int>(n); ++nctr)
        {
            for (int mctr=-nctr; mctr <= nctr; ++mctr)
            {
                double yval = SpharmHolder::Ynm(nctr, mctr, sp.cos_theta, sp.phi);

                // set matrix val
                gsl_matrix_set(&m.matrix, row_ctr, col_ctr, yval);

                col_ctr++;
            }
        }

        row_ctr++;
    }

    // do the SVD
    gsl_vector *solution = gsl_vector_alloc(spharm.size());
    gsl_vector *work = gsl_vector_alloc(spharm.size());
    gsl_matrix *V = gsl_matrix_alloc(spharm.size(), spharm.size());
    gsl_vector *S = gsl_vector_alloc(spharm.size());

    int retval = gsl_linalg_SV_decomp(&(m.matrix), V, S, work);
    if (retval != 0) { std::cerr << "SVD_decomp retval: " << retval << std::endl; }
    assert(retval == 0);
    retval = gsl_linalg_SV_solve(&m.matrix, V, S, &f.vector, solution);
    if (retval != 0) { std::cerr << "SVD_solve retval: " << retval << std::endl; }
    assert(retval == 0);

    // get the results and stick them in the spherical harmonic holder
    for (unsigned int i=0; i < spharm.size(); ++i)
    {
        spharm(i) = gsl_vector_get(solution, i);
    }

    // delete allocated memory
    gsl_vector_free(solution);
    gsl_vector_free(work);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    delete[] a_data;
    delete[] f_data;

    return spharm;
}

#ifdef PYTHON_MODULE
SpharmHolder least_squares_spherical_harmonic_pylist(const boost::python::list& pts, const boost::python::list& vals, const Vector& origin, unsigned int n)
{
    // check same number of vals as points
    assert(len(pts) == len(vals));

    // n is max order of spherical harmonic to combine to approximate the function at the sample_pts
    std::auto_ptr<SpharmHolder> spharm_ptr(new SpharmHolder(n));
    SpharmHolder& spharm = *spharm_ptr;

    // check that the number of parameters to fit is smaller than the number of constraints
    // (i.e. is overspecified problem for SVD/least squares)
    std::cout << "SPHARM SIZE: " << spharm.size() << " CONSTRAINTS: " << len(pts) << std::endl;
    assert(spharm.size() <= len(pts));

    // Allocate memory for matrix
    int num_pts = len(pts);
    double *a_data = new double[num_pts*spharm.size()];
    double *f_data = new double[num_pts];
    gsl_matrix_view m = gsl_matrix_view_array(a_data, num_pts, spharm.size());
    gsl_vector_view f = gsl_vector_view_array(f_data, num_pts);

    // for each spherical harmonic point put an element in the matrix
    for (int row_ctr=0; row_ctr < num_pts; ++row_ctr)
    {
        Vector pt = static_cast<Vector>(extract<Vector>(pts[row_ctr]));
        double fval = static_cast<double>(extract<double>(vals[row_ctr]));
        const SpharmPt sp(pt, origin);

        // set vector element in this row
        gsl_vector_set(&f.vector, row_ctr, fval);

        // set matrix elements in this row
        int col_ctr=0;
        for (int nctr=0; nctr <= static_cast<int>(n); ++nctr)
        {
            for (int mctr=-nctr; mctr <= nctr; ++mctr)
            {
                double yval = SpharmHolder::Ynm(nctr, mctr, sp.cos_theta, sp.phi);

                // set matrix val
                //std::cout << "Setting row " << row_ctr+1 << " (of " << num_pts << ") col " << col_ctr+1 << " (of " <<  spharm.size() << ") with val " << yval;
                gsl_matrix_set(&m.matrix, row_ctr, col_ctr, yval);
                //std::cout << " (done)" << std::endl;

                col_ctr++;
            }
        }
    }

    // do the SVD
    gsl_vector *solution = gsl_vector_alloc(spharm.size());
    gsl_vector *work = gsl_vector_alloc(spharm.size());
    gsl_matrix *V = gsl_matrix_alloc(spharm.size(), spharm.size());
    gsl_vector *S = gsl_vector_alloc(spharm.size());

    int retval = gsl_linalg_SV_decomp(&(m.matrix), V, S, work);
    //int retval = gsl_linalg_SV_decomp_jacobi(&(m.matrix), V, S);
    assert(retval == 0);
    // get the results and stick them in the spherical harmonic holder
    retval = gsl_linalg_SV_solve(&(m.matrix), V, S, &(f.vector), solution);
    assert(retval == 0);

    // get the results and stick them in the spherical harmonic holder
    for (unsigned int i=0; i < spharm.size(); ++i)
    {
        spharm(i) = gsl_vector_get(solution, i);
        //std::cout << spharm(i) << std::endl;
    }

    // delete allocated memory
    gsl_vector_free(solution);
    gsl_vector_free(work);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    delete[] a_data;
    delete[] f_data;

    return spharm;
}

BOOST_PYTHON_MODULE(_SPHARM)
{

    class_<SpharmHolder>("SpharmHolder", init<const SpharmHolder&>())
        .def("evaluate", &SpharmHolder::evaluate)
        .def("dArea", &SpharmHolder::dArea)
        .def("dtheta", &SpharmHolder::dtheta)
        .def("dphi", &SpharmHolder::dphi)
        .def("max_degree", &SpharmHolder::max_degree);

    def("least_squares_spherical_harmonic", least_squares_spherical_harmonic_pylist, "Fit spherical harmonics");

}

#endif

