/*
* fmm.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef FMM_H_
#define FMM_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ltl/marray.h>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "../common/matrix_types.h"
#include "../common/multipole_holder.h"

#include <cassert>
#include <iostream>
#include <vector>
#include <complex>
#include <limits>

#include "../common/octree_indexer.h"
#include "../common/charge.h"
#include "eval_pt.h"
#include "fmm_globals.h"
#include "fmm_math_funcs.h"

#define PRECIS std::numeric_limits<double>::epsilon()

namespace fmm
{

template<int n, int m, unsigned short nmax>
inline void qnm(double beta,
                const DblMatrix1D & bi,
                const TriangularMemory<nmax, double> &p,
                const CmplxMatrix1D& ephi,
                CmplxMatrix2D &q)
{
    assert(n >= 0);
    assert(n <= nmax);
    const int absm = m < 0 ? -m : m;
    assert(absm <= n);

    if (n==0)
    {
        q(n,m) = bi(0);
    }
    else if (m < 0)
    {
        q(n,m) = ((absm % 2) ? -1.0 : 1.0) * conj(bi(n)*p(n,-m)*ephi(-m));
    }
    else
    {
        q(n,m) = bi(n)*p(n,m)*ephi(m);
    }

    return;
}

template< int i, int j, unsigned short nmax>
class Q_LOOP {
public:
    static inline void EXEC(double beta,
                            const DblMatrix1D& bi,
                            const TriangularMemory<nmax, double> &workspace_p,
                            const CmplxMatrix1D &ephi,
                            CmplxMatrix2D &q)
    {
        qnm<i,j>(beta,bi,workspace_p,ephi,q);
        Q_LOOP< i, j-1, nmax >::EXEC(beta,bi,workspace_p,ephi,q);
    }

};

template<int i, unsigned short nmax>
class Q_LOOP<i,0,nmax>{
public:
    static inline void EXEC(double beta,
                            const DblMatrix1D& bi,
                            const TriangularMemory<nmax, double> &workspace_p,
                            const CmplxMatrix1D &ephi,
                            CmplxMatrix2D &q)
    {
        qnm<i,0>(beta,bi,workspace_p,ephi,q);
        Q_LOOP<i-1,i-1,nmax>::EXEC(beta,bi,workspace_p,ephi,q);

    }
};

template<unsigned short nmax>
class Q_LOOP<0,0,nmax>{
public:
    static inline void EXEC(double beta,
                            const DblMatrix1D& bi,
                            const TriangularMemory<nmax, double> &workspace_p,
                            const CmplxMatrix1D &ephi,
                            CmplxMatrix2D &q)
    {
        qnm<0,0>(beta,bi,workspace_p,ephi,q);
    }
};

template<int nt>
class QNM_LOOP {
public:
    static inline void EXEC(double beta,
                            const DblMatrix1D& bi,
                            const TriangularMemory<nt, double> &workspace_p,
                            const CmplxMatrix1D &ephi,
                            CmplxMatrix2D &q)
    {
        Q_LOOP< nt, nt, nt >::EXEC(beta,bi,workspace_p,ephi,q);
    }
};

template<int n, int m>
inline void grad_q(const double beta,
                    const double scale,
                    const double inv_scale,
                    const CmplxMatrix2D &q,
                    CmplxMatrix3D &gradq)
{
    assert(n >= 0);
    assert(m >= 0);
    assert(m <= n);

    // clear gradq
    std::complex<double> left[3];
    std::complex<double> right[3];

    if (n==0 && m==0)
    {
        left[0] = 0.0;
        left[1] = 0.0;
        left[2] = 0.0;
        right[0] = -conj(q(1,1))*scale;
        right[1] =      -q(1,0)*scale;
        right[2] =       q(1,1)*scale;
    }
    else if (n==1 && m==0)
    {
        left[0] = 0.0;
        left[1] = q(0,0)*inv_scale;
        left[2] = 0.0;
        right[0] = -conj(q(2,1))*scale;
        right[1] = -2.0*q(2,0)*scale;
        right[2] = q(2,1)*scale;
    }
    else if (m==0) //&& n > 1
    {
        left[0] = -conj(q(n-1,1))*inv_scale;
        left[1] = q(n-1,0)*static_cast<double>(n)*inv_scale;
        left[2] = q(n-1,1)*inv_scale;
        right[0] = -conj(q(n+1,1))*scale;
        right[1] = -static_cast<double>(n+1)*q(n+1,0)*scale;
        right[2] = q(n+1,1)*scale;
    }
    else if (n==1 && m==1)
    {
        left[0] = 2.0*q(0,0)*inv_scale;
        left[1] = 0.0;
        left[2] = 0.0;
        right[0] = 2.0*q(2,0)*scale;
        right[1] = -q(2,1)*scale;
        right[2] = q(2,2)*scale;
    }
    else if (m==n)
    {
        left[0]= q(n-1,n-1)*(static_cast<double>((n+n-1)*(n+n))*inv_scale);
        left[1] = 0.0;
        left[2] = 0.0;
        right[0] =  q(n+1,n-1)*2.0*scale;
        right[1] = -q(n+1,n)*scale;
        right[2] =  q(n+1,n+1)*scale;
    }
    else if (m==(n-1))
    {
        left[0] = q(n-1,n-2)*(static_cast<double>((n+m-1)*(n+m))*inv_scale);
        left[1] = q(n-1,n-1)*(static_cast<double>(n+m)*inv_scale);
        left[2] = 0.0;
        right[0] =  q(n+1,m-1)*(6.0*scale);
        right[1] = -q(n+1,m)*(2.0*scale);
        right[2] =  q(n+1,m+1)*scale;
    }
    else
    {
        left[0] = q(n-1,m-1)*(static_cast<double>((n+m-1)*(n+m))*inv_scale);
        left[1] = q(n-1,m)*(static_cast<double>(n+m)*inv_scale);
        left[2] = q(n-1,m+1)*inv_scale;
        right[0] = q(n+1,m-1)*((static_cast<double>(n-m+1)*(n-m+2))*scale);
        right[1] = -q(n+1,m)*(static_cast<double>(n-m+1)*scale);
        right[2] = q(n+1,m+1)*scale;
    }

    const double premult = -0.5 / static_cast<double>(2*n + 1);
    gradq(0,n,m) = premult * beta *        ((left[0]-right[0]) - (left[2]-right[2]));
    gradq(1,n,m) = std::complex<double>(0,premult * beta) * ((left[0]-right[0]) + (left[2]-right[2]));
    gradq(2,n,m) = premult * beta * -2.0 *  (left[1]-right[1]);

    return;
}

template< int i, int j>
class G_LOOP {
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix2D &q,
                            CmplxMatrix3D &gqk)
    {
        grad_q<i,j>(beta,scale,inv_scale,q,gqk);
        G_LOOP< i, j-1 >::EXEC(beta,scale,inv_scale,q,gqk);
    }

};

template<int i>
class G_LOOP<i,0>{
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix2D &q,
                            CmplxMatrix3D &gqk)

    {
        grad_q<i,0>(beta,scale,inv_scale,q,gqk);
        G_LOOP<i-1,i-1>::EXEC(beta,scale,inv_scale,q,gqk);

    }
};

template<>
class G_LOOP<0,0>{
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix2D &q,
                            CmplxMatrix3D &gqk)

    {
        grad_q<0,0>(beta,scale,inv_scale,q,gqk);
    }
};

template<int nt>
class GRADQ_LOOP {
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix2D &q,
                            CmplxMatrix3D &gqk)

    {
        G_LOOP< nt, nt >::EXEC(beta,scale,inv_scale,q,gqk);
    }
};

template<int n, int m>
inline void grad2_q(double beta,
                    const double scale,
                    const double inv_scale,
                    const CmplxMatrix3D &gradq,
                    CmplxMatrix4D &grad2q)
{

    assert(n >= 0);
    assert(m >= 0);
    assert(m <= n);

    std::complex<double> left[3];
    std::complex<double> right[3];
    const double premult = -0.5 / static_cast<double>(2*n + 1);

    for (int i=0; i < 3; ++i)
    {
        if (n==0 && m==0)
        {
            left[0] = 0.0;
            left[1] = 0.0;
            left[2] = 0.0;
            right[0] = -conj(gradq(i,1,1))*scale;
            right[1] =      -gradq(i,1,0)*scale;
            right[2] =       gradq(i,1,1)*scale;
        }
        else if (n==1 && m==0)
        {
            left[0] = 0.0;
            left[1] = gradq(i,0,0)*inv_scale;
            left[2] = 0.0;
            right[0] = -conj(gradq(i,2,1))*scale;
            right[1] = -2.0*gradq(i,2,0)*scale;
            right[2] = gradq(i,2,1)*scale;
        }
        else if (m==0) //&& n > 1
        {
            left[0] = -conj(gradq(i,n-1,1))*inv_scale;
            left[1] = gradq(i,n-1,0)*static_cast<double>(n)*inv_scale;
            left[2] = gradq(i,n-1,1)*inv_scale;
            right[0] = -conj(gradq(i,n+1,1))*scale;
            right[1] = -static_cast<double>(n+1)*gradq(i,n+1,0)*scale;
            right[2] = gradq(i,n+1,1)*scale;
        }
        else if (n==1 && m==1)
        {
            left[0] = 2.0*gradq(i,0,0)*inv_scale;
            left[1] = 0.0;
            left[2] = 0.0;
            right[0] = 2.0*gradq(i,2,0)*scale;
            right[1] = -gradq(i,2,1)*scale;
            right[2] = gradq(i,2,2)*scale;
        }
        else if (m==n)
        {
            left[0]= gradq(i,n-1,n-1)*(static_cast<double>((n+n-1)*(n+n))*inv_scale);
            left[1] = 0.0;
            left[2] = 0.0;
            right[0] =  gradq(i,n+1,n-1)*2.0*scale;
            right[1] = -gradq(i,n+1,n)*scale;
            right[2] =  gradq(i,n+1,n+1)*scale;
        }
        else if (m==(n-1))
        {
            left[0] = gradq(i,n-1,n-2)*(static_cast<double>((n+m-1)*(n+m))*inv_scale);
            left[1] = gradq(i,n-1,n-1)*(static_cast<double>(n+m)*inv_scale);
            left[2] = 0.0;
            right[0] =  gradq(i,n+1,m-1)*(6.0*scale);
            right[1] = -gradq(i,n+1,m)*(2.0*scale);
            right[2] =  gradq(i,n+1,m+1)*scale;
        }
        else
        {
            left[0] = gradq(i,n-1,m-1)*(static_cast<double>((n+m-1)*(n+m))*inv_scale);
            left[1] = gradq(i,n-1,m)*(static_cast<double>(n+m)*inv_scale);
            left[2] = gradq(i,n-1,m+1)*inv_scale;
            right[0] = gradq(i,n+1,m-1)*((static_cast<double>(n-m+1)*(n-m+2))*scale);
            right[1] = -gradq(i,n+1,m)*(static_cast<double>(n-m+1)*scale);
            right[2] = gradq(i,n+1,m+1)*scale;
        }

        grad2q(0,i,n,m) = premult * beta *        ((left[0]-right[0]) - (left[2]-right[2]));
        grad2q(1,i,n,m) = std::complex<double>(0,premult * beta)  * ((left[0]-right[0]) + (left[2]-right[2]));
        grad2q(2,i,n,m) = premult * beta * -2.0 *  (left[1]-right[1]);

    }

    return;
}

template< int i, int j>
class G2_LOOP {
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix3D &gqk,
                            CmplxMatrix4D &g2qk)
    {
        grad2_q<i,j>(beta,scale,inv_scale,gqk,g2qk);
        G2_LOOP< i, j-1 >::EXEC(beta,scale,inv_scale,gqk,g2qk);
    }

};

template<int i>
class G2_LOOP<i,0>{
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix3D &gqk,
                            CmplxMatrix4D &g2qk)
    {
        grad2_q<i,0>(beta,scale,inv_scale,gqk,g2qk);
        G2_LOOP<i-1,i-1>::EXEC(beta,scale,inv_scale,gqk,g2qk);

    }
};

template<>
class G2_LOOP<0,0>{
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix3D &gqk,
                            CmplxMatrix4D &g2qk)
    {
        grad2_q<0,0>(beta,scale,inv_scale,gqk,g2qk);
    }
};

template<int nt>
class GRAD2Q_LOOP {
public:
    static inline void EXEC(double beta,
                            double scale,
                            double inv_scale,
                            const CmplxMatrix3D &gqk,
                            CmplxMatrix4D &g2qk)
    {
        G2_LOOP< nt, nt >::EXEC(beta,scale,inv_scale,gqk,g2qk);
    }
};


template<int NTERMS>
inline void evaluate_local_expansion_at_xyz(const double beta,
                                            double scale,
                                            EvalPt_2ndDerivs& eval_pt,
                                            const Vector& x0y0z0,
                                            const BaseMultipoleHolder<NTERMS> &locals)
{

    double pot=0;
    Vector field(0,0,0);
    GradField3x3 field2(0.0);
    
    const Vector& xyz = eval_pt.pt();

    // workspace for legendre polynomials
    TriangularMemory<NTERMS+2,double> workspace_p;
    std::complex<double> another_workspace[6];

    // need this for the derivative terms
    const double inv_scale = 1.0 / scale;

    const Vector r = (xyz - x0y0z0);

    double proj = r.x*r.x+r.y*r.y;
    const double rr = proj+r.z*r.z;
    proj = sqrt(proj);
    const double d = sqrt(rr);
    double ctheta;

    if (d <= PRECIS) {
        ctheta = 0.0;
    }
    else {
        ctheta = r.z/d;
    }

    CmplxMatrix1D ephi(Range(0,NTERMS+4));
    ephi(0) = std::complex<double>(1.0,0.0);
    if ( proj <= PRECIS*d ) {
        ephi(1) = 1.0;
    }
    else {
        ephi(1) = std::complex<double>(r.x/proj,r.y/proj);
    }

    for (int i=1; i <= NTERMS+2; ++i)
    {
        ephi(i+1) = ephi(i)*ephi(1);
    }

    const double rk=d*beta;

    DblMatrix1D bi(Range(0,NTERMS+3));
    i_n(scale,rk,NTERMS+2,bi);

    // could use GSL instead?
    lgndr(ctheta, workspace_p);

    pot = 0.0;
    field.x = 0.0;
    field.y = 0.0;
    field.z = 0.0;
    for (int i=0; i < 3; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
            field2(i,j)=0;
        }
    }

    // value of qk for all combos of n and m
    CmplxMatrix2D qk(Range(0,NTERMS+2),Range(0,NTERMS+2));
    CmplxMatrix3D gqk(Range(0,3),Range(0,NTERMS+1),Range(0,NTERMS+1));
    CmplxMatrix4D g2qk(Range(0,3),Range(0,3),Range(0,NTERMS),Range(0,NTERMS));

    QNM_LOOP<NTERMS+2>::EXEC(beta,bi,workspace_p,ephi,qk);
    GRADQ_LOOP<NTERMS+1>::EXEC(beta,scale,inv_scale,qk,gqk);
    GRAD2Q_LOOP<NTERMS>::EXEC(beta,scale,inv_scale,gqk,g2qk);

    for (int n=0; n <=NTERMS; ++n)
    {
        for (int m=1; m <= n; ++m)
        {
            pot += std::real(qk(n,m) * locals(n,m));
            field.x += std::real(gqk(0,n,m)*locals(n,m));
            field.y += std::real(gqk(1,n,m)*locals(n,m));
            field.z += std::real(gqk(2,n,m)*locals(n,m));
            for (int i=0; i < 3; ++i)
            {
                for (int j=0; j < 3; ++j)
                {
                    field2(i,j) += std::real(g2qk(i,j,n,m) * locals(n,m));
                }
            }

        }
    }

    pot *= 2.0;
    field.x *= 2.0;
    field.y *= 2.0;
    field.z *= 2.0;
    for (int i=0; i < 3; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
            field2(j,i) *= 2.0;
        }
    }

    for (int n=0; n <=NTERMS; ++n)
    {
        pot += std::real(qk(n,0) * locals(n,0));
        field.x += std::real(gqk(0,n,0)*locals(n,0));
        field.y += std::real(gqk(1,n,0)*locals(n,0));
        field.z += std::real(gqk(2,n,0)*locals(n,0));
        for (int i=0; i < 3; ++i)
        {
            for (int j=0; j < 3; ++j)
            {
                field2(i,j) += std::real(g2qk(i,j,n,0) * locals(n,0));
            }
        }
    }

    const double two_beta_over_pi = beta*TWO_OVER_PI;
    pot *= two_beta_over_pi;
    field.x *= two_beta_over_pi;
    field.y *= two_beta_over_pi;
    field.z *= two_beta_over_pi;
    for (int i=0; i < 3; ++i)
    {
        for (int j=0; j < 3; ++j)
        {
            field2(i,j) *= two_beta_over_pi;
        }
    }

    eval_pt.add_potential(pot);
    eval_pt.add_field(field);
    eval_pt.add_field2(field2);

    return;

}

template<int NTERMS>
inline void evaluate_local_expansion_at_xyz(const double beta,
                                            double scale,
                                            EvalPtBase& eval_pt,
                                            const Vector& x0y0z0,
                                            const BaseMultipoleHolder<NTERMS> &locals)
{

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    evaluates local expansion at arbitrary point.
    //
    //  on input:
    //
    //    beta : the frequency of the equation.
    //    local: coefficients of local expansion(scaled)
    //    x0y0z0: the center of the expansion
    //    point: point of evaluation
    //    NTERMS: order of expansion
    //    p: work arrays to hold legendre polynomials
    //       and associated legendre functions, respectively.
    //
    //  on output:
    //    rpot: computed potential
    //
    //  note: the current version will compute the potential only,
    //    and the computation of the field will be added later.
    //
    //   --------------------------------------------------
    //  subroutine called :
    //
    //  called from : brfrc()
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    double pot=0;
    Vector field(0,0,0);
    
    const Vector& xyz = eval_pt.pt();

    // workspace for legendre polynomials
    TriangularMemory<NTERMS+1,double> workspace_p;

    const double inv_scale = 1.0 / scale;
    const Vector r = (xyz - x0y0z0);

    double proj = r.x*r.x+r.y*r.y;
    const double rr = proj+r.z*r.z;
    proj = sqrt(proj);
    const double d = sqrt(rr);
    double ctheta;

    if (d <= PRECIS) {
        ctheta = 0.0;
    }
    else {
        ctheta = r.z/d;
    }

    CmplxMatrix1D ephi(Range(0,NTERMS+4));
    ephi(0) = std::complex<double>(1.0,0.0);
    if ( proj <= PRECIS*d ) {
        ephi(1) = 1.0;
    }
    else {
        ephi(1) = std::complex<double>(r.x/proj,r.y/proj);
    }

    for (int i=1; i <= NTERMS+2; ++i)
    {
        ephi(i+1) = ephi(i)*ephi(1);
    }

    const double rk=d*beta;

    DblMatrix1D bi(Range(0,NTERMS+2));
    i_n(scale,rk,NTERMS+1,bi);

    // could use GSL instead?
    lgndr(ctheta, workspace_p);

    pot = 0.0;
    field.x = 0.0;
    field.y = 0.0;
    field.z = 0.0;


    // value of qk for all combos of n and m
    CmplxMatrix2D qk(Range(0,NTERMS+1),Range(0,NTERMS+1));
    CmplxMatrix3D gqk(Range(0,3),Range(0,NTERMS),Range(0,NTERMS));

    QNM_LOOP<NTERMS+1>::EXEC(beta,bi,workspace_p,ephi,qk);
    GRADQ_LOOP<NTERMS>::EXEC(beta,scale,inv_scale,qk,gqk);

    for (int n=0; n <=NTERMS; ++n)
    {
        for (int m=1; m <= n; ++m)
        {
            pot += std::real(qk(n,m) * locals(n,m));
            field.x += std::real(gqk(0,n,m)*locals(n,m));
            field.y += std::real(gqk(1,n,m)*locals(n,m));
            field.z += std::real(gqk(2,n,m)*locals(n,m));
        }
    }

    pot *= 2.0;
    field.x *= 2.0;
    field.y *= 2.0;
    field.z *= 2.0;

    for (int n=0; n <=NTERMS; ++n)
    {
        pot += std::real(qk(n,0) * locals(n,0));
        field.x += std::real(gqk(0,n,0)*locals(n,0));
        field.y += std::real(gqk(1,n,0)*locals(n,0));
        field.z += std::real(gqk(2,n,0)*locals(n,0));
    }

    const double two_beta_over_pi = beta*TWO_OVER_PI;
    pot *= two_beta_over_pi;
    field.x *= two_beta_over_pi;
    field.y *= two_beta_over_pi;
    field.z *= two_beta_over_pi;

    eval_pt.add_potential(pot);
    eval_pt.add_field(field);
    
    return;
}

template<int NTERMS>
inline void evaluate_local_expansion_at_xyz(double beta_kappa,
                                    double beta_kappa0,
                                    double scale_kappa,
                                    double scale_kappa0,
                                    const Vector& xyz,
                                    const Vector& x0y0z0,
                                    const MultiHolder<12,BaseMultipoleHolder<NTERMS> > &locals_12mer,
                                    DblMatrix1D& pot_twelvelet,
                                    DblMatrix2D& field_twelvelet,
                                    DblMatrix3D& field2_twelvelet
                                    )
{

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    evaluates local expansion at arbitrary point.
    //
    //  on input:
    //
    //    beta : the frequency of the equation.
    //    local: coefficients of local expansion(scaled)
    //    x0y0z0: the center of the expansion
    //    point: point of evaluation
    //    NTERMS: order of expansion
    //    p: work arrays to hold legendre polynomials
    //       and associated legendre functions, respectively.
    //
    //  on output:
    //    rpot: computed potential
    //
    //  note: the current version will compute the potential only,
    //    and the computation of the field will be added later.
    //
    //   --------------------------------------------------
    //  subroutine called :
    //
    //  called from : brfrc()
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    // workspace for legendre polynomials
    TriangularMemory<NTERMS+2,double> workspace_p;
    Vector r = (xyz - x0y0z0);

    double proj = r.x*r.x+r.y*r.y;
    double rr = proj+r.z*r.z;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta;

    if (d <= PRECIS) {
        ctheta = 0.0;
    }
    else {
        ctheta = r.z/d;
    }

    assert(fabs(ctheta) <= 1.0);

    CmplxMatrix1D ephi(Range(0,NTERMS+4));
    ephi(0) = std::complex<double>(1.0,0.0);
    if ( proj <= PRECIS*d ) {
        ephi(1) = 1.0;
    }
    else {
        //double phi = atan2(r.y,r.x);
        //ephi(1) = std::complex<double>(cos(phi),sin(phi));
        ephi(1) = std::complex<double>(r.x/proj,r.y/proj);
    }

    for (int i=1; i <= NTERMS+2; ++i)
    {
        ephi(i+1) = ephi(i)*ephi(1);
    }

    DblMatrix1D bi(Range(0,NTERMS+3));

    // value of qk for all combos of n and m
    CmplxMatrix2D qk(Range(0,NTERMS+2),Range(0,NTERMS+2));
    CmplxMatrix3D gqk(Range(0,3),Range(0,NTERMS+1),Range(0,NTERMS+1));
    CmplxMatrix4D g2qk(Range(0,3),Range(0,3),Range(0,NTERMS),Range(0,NTERMS));

    std::complex<double> tmp[3];

    // beta_use depends on which charge set we're working on
    for (unsigned short BEM_12way_idx=0; BEM_12way_idx < 12; ++BEM_12way_idx)
    {
        const double beta = (BEM_12way_idx >= 8) ? beta_kappa0 : beta_kappa;
        const double scale = (BEM_12way_idx >= 8) ? scale_kappa0 : scale_kappa;

        double& pot = pot_twelvelet(BEM_12way_idx);
        DblMatrix1D field_mtrx = field_twelvelet(Range(0,3),BEM_12way_idx);
        field_mtrx.setBase(0);
        DblMatrix2D field2 = field2_twelvelet(Range(0,3),Range(0,3),BEM_12way_idx);
        field2.setBase(0,0);

        pot=0.0;
        Vector field(0,0,0);
        for (int i=0; i < 3; ++i)
        {
            for (int j=0; j < 3; ++j)
            {
                field2(j,i) = 0.0;
            }
        }

        // generate the qnm/grad_qnm/grad2_qnm terms only twice -- once for
        // each type of beta
        if (BEM_12way_idx == 0 || (BEM_12way_idx == 8 && beta_kappa0 != beta_kappa) )
        {

            const double inv_scale = 1.0 / scale;
            const double rk=d*beta;

            DblMatrix1D bi(Range(0,NTERMS+3));
            i_n(scale,rk,NTERMS+2,bi);

            // could use GSL instead?
            lgndr(ctheta, workspace_p);

            QNM_LOOP<NTERMS+2>::EXEC(beta,bi,workspace_p,ephi,qk);
            GRADQ_LOOP<NTERMS+1>::EXEC(beta,scale,inv_scale,qk,gqk);
            GRAD2Q_LOOP<NTERMS>::EXEC(beta,scale,inv_scale,gqk,g2qk);

        }

        const BaseMultipoleHolder<NTERMS> &locals = locals_12mer[BEM_12way_idx];
        for (int n=0; n <=NTERMS; ++n)
        {
            for (int m=1; m <= n; ++m)
            {
                pot += std::real(qk(n,m) * locals(n,m));
                field.x += std::real(gqk(0,n,m)*locals(n,m));
                field.y += std::real(gqk(1,n,m)*locals(n,m));
                field.z += std::real(gqk(2,n,m)*locals(n,m));
                for (int i=0; i < 3; ++i)
                {
                    for (int j=0; j < 3; ++j)
                    {
                        field2(i,j) += std::real(g2qk(i,j,n,m) * locals(n,m));
                    }
                }

            }
        }

        pot *= 2.0;
        field.x *= 2.0;
        field.y *= 2.0;
        field.z *= 2.0;
        for (int i=0; i < 3; ++i)
        {
            for (int j=0; j < 3; ++j)
            {
                field2(j,i) *= 2.0;
            }
        }

        for (int n=0; n <=NTERMS; ++n)
        {
            pot += std::real(qk(n,0) * locals(n,0));
            field.x += std::real(gqk(0,n,0)*locals(n,0));
            field.y += std::real(gqk(1,n,0)*locals(n,0));
            field.z += std::real(gqk(2,n,0)*locals(n,0));
            for (int i=0; i < 3; ++i)
            {
                for (int j=0; j < 3; ++j)
                {
                    field2(i,j) += std::real(g2qk(i,j,n,0) * locals(n,0));
                }
            }
        }

        const double two_beta_over_pi = beta*TWO_OVER_PI;
        pot *= two_beta_over_pi;
        field.x *= two_beta_over_pi;
        field.y *= two_beta_over_pi;
        field.z *= two_beta_over_pi;
        for (int i=0; i < 3; ++i)
        {
            for (int j=0; j < 3; ++j)
            {
                field2(i,j) *= two_beta_over_pi;
            }
        }

        field_mtrx(0) = field.x;
        field_mtrx(1) = field.y;
        field_mtrx(2) = field.z;

    }

    return;
}

template<int NLAMBS, int NWAVES>
inline void yphystof(const int nexptot,
            const int nexptotp,
            const int fsize,
            const DblMatrix1D& quad_nodes,
            const UShrtMatrix1D &numfour,
            const UShrtMatrix1D &numphys,
            const CmplxMatrix1D& fexpback,
            const BasePlaneWaveHolder<NLAMBS,NWAVES>& mexpphys,
            CmplxMatrix1D& mexpf)
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    this subroutine converts the discretized exponential moment function
    //    into its fourier expansion.
    //    it calculates the inner sum of the exp->local expansion.
    //      (/sum_{j=1}^{m(k)} w(k,j)*e^{-im*alpha_j})/m(k) for k=1, NLAMBS,
    //      and m=0, numfour.
    //      numfour is the total number of the fourier modes.
    //      or in other words, those l_n^m <>0.
    //      the summation is over the numphys.
    //
    //  on input:
    //
    //    mexpphys(*):  discrete values of the moment function
    //                  m(\lambda,\alpha), ordered as follows.
    //
    //        mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
    //             m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
    //        mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
    //             m(\lambda_2,0),...,
    //                 m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
    //        etc.
    //
    //    NLAMBS:        number of discretization pts. in lambda integral
    //    quad_nodes(NLAMBS): discretization points in lambda integral.
    //    numfour(j):   number of fourier modes in the expansion
    //                      of the function m(\lambda_j,\alpha)
    //    fexpback : contains the precomputed e^{-im *alpha_j}
    //
    //  on output:
    //
    //    mexpf(*):     fourier coefficients of the function
    //                  mexp(lambda,m) for discrete lambda values.
    //                  they are ordered as follows:
    //
    //               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
    //               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
    //                                              for lambda_2
    //               etc.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    int next = 0;
    int nftot = 0;
    int nptot = 0;

    for (int i=0; i < NLAMBS; ++i)
    {
        int nalpha = numphys(i);
        int nalpha2 = nalpha/2;

        //
        //-------first mm=0 case.
        //
        mexpf(nftot) = std::complex<double>(0.0,0.0);
        for (int ival=0; ival < nalpha2; ++ival)
        {
            mexpf(nftot) += 2.0*std::real(mexpphys(nptot+ival));
        }
        mexpf(nftot) /= double(nalpha);

        //
        //-------even mm.  (w(k,j)+conj(w(k,j))*e^{-im alpha_j}
        //
        for (int mm=3; mm <= numfour(i); mm+=2)
        {
            mexpf(nftot+mm-1) = std::complex<double>(0.0,0.0);
            for (int ival=0; ival < nalpha2; ++ival)
            {
                double rtmp = 2.0 * std::real(mexpphys(nptot+ival));
                mexpf(nftot+mm-1) += fexpback(next)*rtmp;
                next++;
            }
            mexpf(nftot+mm-1) /= double(nalpha);
        }

        //
        //-------odd mm. (w(k,j)-conj(w(k,j))*e^{-im alpha_j}
        //
        for (int mm=2; mm <= numfour(i); mm+=2)
        {
            mexpf(nftot+mm-1) = std::complex<double>(0.0,0.0);
            for (int ival=0; ival < nalpha2; ++ival)
            {
                std::complex<double> ztmp = std::complex<double>(0.0, 2.0 * std::imag(mexpphys(nptot+ival)));
                mexpf(nftot+mm-1) += fexpback(next)*ztmp;
                next++;
            }
            mexpf(nftot+mm-1) /= double(nalpha);
        }

        nftot += numfour(i);
        nptot += numphys(i)/2;
    }

    return;
}

template<int NLAMBS, int NTERMS>
inline void yexptolocal(const int nexptot,
                double beta,
                const DblMatrix3D& rlsc,
                BaseMultipoleHolder<NTERMS>& local,
                const DblMatrix1D& quad_nodes,
                const DblMatrix1D& quad_weights,
                const DblMatrix2D& scale_factors,
                const UShrtMatrix1D &numfour,
                const CmplxMatrix1D& mexpup,
                const CmplxMatrix1D& mexpdown)
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //    this subroutine converts the fourier representation of two
    //    exponential moment functions into a local multipole expansion
    //    (with respect to the same box center).
    //      l_n^m= (see reference). and is scaled.
    //
    //    u(x,y,z) = \int_0^\infty e^{-\lambda z}
    //                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
    //                mexpup(lambda,alpha) dalpha dlambda
    //            +
    //                \int_0^\infty e^{\lambda z}
    //                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
    //                mexpdown(lambda,alpha) dalpha dlambda
    //
    //             = \sum_{n=0}^{NTERMS} \sum_{m=-n,n}
    //                local(n,m) y_n^m(cos theta) e^{i m \phi} r^{n}
    //
    //  on input:
    //    beta : the scaled beta used in the calculation of l_n^m.
    //    rlsc : the precomputed and scaled p_n^m *sc^n
    //    NTERMS : the total number of expansions.
    //    mexpup(nexptot): fourier coefficients of the function
    //                    mexpup for discrete lambda
    //                    values. they are ordered as follows:
    //
    //                 mexpup(1,...,numtets(1)) = fourier modes
    //                             for lambda_1
    //                 mexpup(numtets(1)+1,...,numtets(2)) = fourier modes
    //                             for lambda_2
    //                 etc.
    //    mexpdown(nexptot): as above for down expansion
    //    quad_nodes(NLAMBS): discretization points in lambda integral
    //    quad_weights(NLAMBS): quadrature weights in lambda integral
    //    NLAMBS:      number of discretization pts. in lambda integral
    //    numtets(j): number of fourier modes in expansion of alpha
    //                variable for lambda_j.
    //    nexptot:    sum_j numtets(j)
    //
    //  on output:
    //    local(0:NTERMS,0:NTERMS): output multipole expansion of order
    //                              NTERMS.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //

    //
    //-----compute necessary powers of -i
    //
    CmplxMatrix1D zeye(Range(0,NTERMS+1));
    CmplxMatrix1D mexpplus(Range(0,nexptot));
    CmplxMatrix1D mexpminus(Range(0,nexptot));

    zeye(0) = std::complex<double>(1.0, 0.0);
    for (int i=1; i <= NTERMS; ++i)
    {
        zeye(i) = zeye(i-1) * ima;
    }

    //
    //-----initialize local expansion
    //
    for (int nm=0; nm <= NTERMS; ++nm)
    {
        for (int mth=0; mth <= nm; ++mth)
        {
            local(nm,mth) = std::complex<double>(0.0,0.0);
        }
    }

    //
    //-----compute sum and difference of mexpup and mexpdown
    //
    for (int nm=0; nm < nexptot; ++nm)
    {
        mexpplus(nm) = mexpdown(nm) + mexpup(nm);
        mexpminus(nm) = mexpdown(nm) - mexpup(nm);
    }

    //
    //-----loop over multipole order to generate mexp values.
    //
    int ntot = 0;
    for (int nl=0; nl < NLAMBS; ++nl)
    {

        //
        //-------add contributions to local expansion. first compute
        //       p_n^m*w_k*mexplus/minus.
        //
        for (int nm=0; nm <= NTERMS; nm+=2)
        {
            const int mmax = numfour(nl)-1 > nm ? nm : numfour(nl)-1;
            for (int mth=0; mth <= mmax; ++mth)
            {
                local(nm,mth) += rlsc(nm,mth,nl)*quad_weights(nl)*mexpplus(ntot+mth);
            }
        }

        for (int nm=1; nm <= NTERMS; nm+=2)
        {
            const int mmax = numfour(nl)-1 > nm ? nm : numfour(nl)-1;
            for (int mth=0; mth <= mmax; ++mth)
            {
                local(nm,mth) += rlsc(nm,mth,nl)*quad_weights(nl)*mexpminus(ntot+mth);
            }
        }
        ntot += numfour(nl);
    }

    //
    //-----scale the expansions according to formula
    //
    double rscale = pi / (beta*2.0);
    for (int nm=0; nm <=NTERMS; ++nm)
    {
        for (int mth=0; mth <= nm; ++mth)
        {
            local(nm,mth) *= zeye(mth)*rscale*scale_factors(nm,mth);
        }
    }

    return;
}

template<int NLAMBS, int NTERMS>
inline void ympoletoexp(const int nexptot,
                const UShrtMatrix1D &numfour,
                const DblMatrix3D& rlsc,
                const BaseMultipoleHolder<NTERMS>& mpole,
                CmplxMatrix1D& mexpup,
                CmplxMatrix1D& mexpdn)
{

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //     this subroutine converts a multipole expansion mpole into the
    //     corresponding exponential moment function mexp for the
    //     both the +z direction and the -z direction.
    //
    //     u(x,y,z) = \sum_{n=0}^{NTERMS} \sum_{m=-n,n}
    //                mpole(n,m) y_n^m(cos theta) e^{i m \phi}/r^{n+1}
    //
    //              = (1/2pi) \int_0^\infty e^{-\lambda z}
    //                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
    //                mexpup(lambda,alpha) dalpha dlambda
    //
    //     for +z direction and
    //
    //              = (1/2pi) \int_0^\infty e^{\lambda z}
    //                \int_0^{2\pi} e^{-i\lambda(xcos(alpha)+ysin(alpha))}
    //                mexpdown(lambda,alpha) dalpha dlambda
    //
    //     for -z direction.
    //
    //     note: the expression for the -z direction corresponds to the
    //     mapping (x,y,z) -> (-x,-y,-z), i.e. reflection through the origin.
    //     one could also use rotation about the y axis, for which
    //     (x,y,z) -> (-x,y,-z) but we stick to the reflected convention.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //     note: the multipole expansion is assumed to have been rescaled
    //           so that the box containing sources has unit dimension.
    //
    //     note: we only store mpole(n,m) for n,m >= 0, since mpole(n,-m)=
    //           conj(mpole(n,m)). since we store the exponential
    //           moment function in the fourier domain (w.r.t. the alpha
    //           variable), we compute
    //
    //       m_lambda(m) = (i)**m \sum_{n=m}^n c(n,m) mpole(n,m) lambda^n
    //
    //           for m >= 0 only, where c(n,m) = 1/sqrt((n+m)!(n-m)!).
    //
    //       for possible future reference, it should be noted that
    //       it is not true that m_lamb(-m) = conj(m_lamb(m)).
    //       inspection of the integral formula for y_n^{-m} shows that
    //       m_lamb(-m) = conj(m_lamb(m)) * (-1)**m.
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //     on input:
    //
    //     mpole(0:NTERMS,0:NTERMS): the multipole expansion
    //
    //     quad_nodes(NLAMBS):  discretization points in lambda integral
    //
    //     NLAMBS:         number of discretization pts.
    //
    //     numfour(NLAMBS): number of fourier modes needed in expansion
    //                    of alpha variable for each lambda value.
    //                    note : the numfour is given by numthehalf().
    //
    //     nexptot =      sum_j numfour(j)
    //
    //     rlsc() : p_n^m for different lambda_k
    //
    //     on output:
    //
    //     mexpf(nexptot): fourier coefficients of the function
    //                     mexp(lambda,alpha) for successive discrete
    //                     lambda values. they are ordered as follows:
    //
    //                 mexpf(1,...,numfour(1)) = fourier modes
    //                             for lambda_1
    //                 mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
    //                             for lambda_2
    //                 etc.
    //     note by huangjf : in return, we will output
    //       in mexpup, sum_{n=m}^{NTERMS} m_n^m*p_n^m. (all are scaled)
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //     note by jingfang :
    //     1.this subroutine will compute the inner sum.
    //       instead of compute all the NTERMS modes, only
    //       the necessary modes are calculated. the number of mode
    //       needed are provided by the subroutine numthehalf().
    //       note that the number is always less than NTERMS needed,
    //       and the worst case is NTERMS+1?.
    //
    //     2.
    //       subroutine called :
    //       called from : mkudexp(), mknsexp(), mkewexp()
    //
    //     3. the down-list will have the same fourier modes as the
    //        up-list if we only change the sign of the z. so we don't need to
    //        compute them separately.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    //
    //-----loop over multipole order to generate mexpup and mexpdown values.
    //
    int ntot = 0;
    int ncurrent = 0;

    for (int nl = 0; nl < NLAMBS; ++nl)
    {
        double sgn = -1.0;
        for (int mth=0; mth < numfour(nl); ++mth)
        {
            ncurrent = ntot+mth;
            std::complex<double> ztmp1(0.0, 0.0);
            std::complex<double> ztmp2(0.0, 0.0);
            sgn = -sgn;
            for (int nm=mth; nm <= NTERMS; nm+=2)
            {
                ztmp1 += rlsc(nm,mth,nl)*mpole(nm,mth);
            }
            for (int nm=mth+1; nm <= NTERMS; nm+=2)
            {
                ztmp2 += rlsc(nm,mth,nl)*mpole(nm,mth);
            }
            mexpup(ncurrent) = ztmp1 + ztmp2;
            mexpdn(ncurrent) = sgn*(ztmp1 - ztmp2);
        }
        ntot += numfour(nl);
    }

    return;
}

template<int NLAMBS, int NWAVES>
inline void yftophys(const int nexptot,
            const int nexptotp,
            const int fsize,
            const DblMatrix1D& quad_nodes,
            const UShrtMatrix1D &numfour,
            const UShrtMatrix1D &numphys,
            const CmplxMatrix1D& fexpe,
            const CmplxMatrix1D& fexpo,
            const CmplxMatrix1D& mexpf,    // input
            BasePlaneWaveHolder<NLAMBS,NWAVES>& mexpphys        // output
            )
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    this subroutine evaluates the fourier expansion of the
    //    exponential moment function m(\lambda,\alpha) at equispaced
    //    nodes.
    //
    //  on input:
    //
    //    mexpf(*):     fourier coefficients of the function
    //                  mexp(lambda,alpha) for discrete lambda values.
    //                  they are ordered as follows:
    //
    //               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
    //               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
    //                                              for lambda_2
    //               etc.
    //
    //    NLAMBS:        number of discretization pts. in lambda integral
    //    quad_nodes(NLAMBS): discretization points in lambda integral.
    //    numfour(j):   number of fourier modes in the expansion
    //                      of the function m(\lambda_j,\alpha)
    //    numphys : number of fourier modes in the plane wave expansion.
    //    fexpe =      precomputed array of exponentials needed for
    //                 fourier series evaluation. even terms.
    //    fexpo =      precomputed array of exponentials needed for
    //                 fourier series evaluation. odd terms.
    //  note : we will keep these two terms because in
    //         the helmholtz equation, we will need these.
    //         however, in yukawa, it is not necessary to have
    //         them separated.--huangjf
    //
    //  on output:
    //    mexpphys(*):  discrete values of the moment function
    //                  m(\lambda,\alpha), ordered as follows.
    //
    //        mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
    //             m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
    //        mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
    //             m(\lambda_2,0),...,
    //                 m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
    //        etc.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //      this subroutine computes the outer sum, it is possible
    //      to do this using fft. but the current version will not do that.
    //
    //  subroutine called :
    //
    //  called from : mkudexp(), mknsexp(), mkewexp()
    //
    //  note :
    //    the current subroutine computes sum_{m=-numfour, numfour}
    //      e^{im*alpha} * i^|m| *inner(m)
    //
    //    the constant will left to the pw_local.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    int nftot = 0;
    int nptot = 0;
    int nexte = 0;
    int nexto = 0;
    int sgn;

    for(int i=0; i < NLAMBS; ++i)
    {
        for (int ival=0; ival < numphys(i)/2; ++ival)
        {
            mexpphys(nptot+ival) = mexpf(nftot);
            sgn=-2;
            for (int mm=1; mm < numfour(i); mm+=2)
            {
                sgn = -sgn;
                double rtmp = sgn*std::real(fexpe(nexte)*mexpf(nftot+mm));
                nexte++;
#ifdef DELETED
                double& cmplx_part = mexpphys(nptot+ival).imag();
                cmplx_part += rtmp;
#else
                std::complex<double>& data = mexpphys(nptot+ival);
				data.imag(data.imag() + rtmp);
#endif
            }

            sgn=2;
            for (int mm=2; mm < numfour(i); mm+=2)
            {
                sgn = -sgn;
                double rtmp = sgn*std::real( fexpo(nexto)*mexpf(nftot+mm) );
                nexto++;
#ifdef DELETED
                double& real_part = mexpphys(nptot+ival).real();
                real_part += rtmp;
#else
                std::complex<double>& data = mexpphys(nptot+ival);
				data.real(data.real() + rtmp);
#endif
            }
        }

        nftot = nftot+numfour(i);
        nptot = nptot+numphys(i)/2;
    }
    return;
}

#ifdef DELETED
template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void fmm::convert_mpole_to_six_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                const MultiHolder<multiplic,MultipoleHolderT>& mpole,
                MultiHolder<multiplic,PlaneWaveHolderT>& up,
                MultiHolder<multiplic,PlaneWaveHolderT>& down,
                MultiHolder<multiplic,PlaneWaveHolderT>& north,
                MultiHolder<multiplic,PlaneWaveHolderT>& south,
                MultiHolder<multiplic,PlaneWaveHolderT>& east,
                MultiHolder<multiplic,PlaneWaveHolderT>& west)
{

    MultipoleHolderT mpole_rotation;
    #pragma omp parallel for private(mpole_rotation)
    for (int m=begin; m < end; ++m)
    {
        // convert multipole expansion to planewave expansions
        convert_mp_to_exp(fmm_globs, fmm_level_globs, mpole[m], up[m], down[m]);

        rotate_ztoy(mpole[m], mpole_rotation, fmm_globs.rdmpi2);
        convert_mp_to_exp(fmm_globs, fmm_level_globs, mpole_rotation, north[m], south[m]);

        rotate_ztox(mpole[m], mpole_rotation, fmm_globs.rdpi2);
        convert_mp_to_exp( fmm_globs, fmm_level_globs, mpole_rotation, east[m], west[m]);
    }
    return;

}
#else // ! DELETED
template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void convert_mpole_to_six_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                const MultiHolder<multiplic,MultipoleHolderT>& mpole,
                MultiHolder<multiplic,PlaneWaveHolderT>& up,
                MultiHolder<multiplic,PlaneWaveHolderT>& down,
                MultiHolder<multiplic,PlaneWaveHolderT>& north,
                MultiHolder<multiplic,PlaneWaveHolderT>& south,
                MultiHolder<multiplic,PlaneWaveHolderT>& east,
                MultiHolder<multiplic,PlaneWaveHolderT>& west);
#endif // ! DELETED

template <int multiplic, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void convert_mpole_to_six_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                const MultiHolder<multiplic,MultipoleHolderT>& mpole,
                MultiHolder<multiplic,PlaneWaveHolderT>& up,
                MultiHolder<multiplic,PlaneWaveHolderT>& down,
                MultiHolder<multiplic,PlaneWaveHolderT>& north,
                MultiHolder<multiplic,PlaneWaveHolderT>& south,
                MultiHolder<multiplic,PlaneWaveHolderT>& east,
                MultiHolder<multiplic,PlaneWaveHolderT>& west)
{
    convert_mpole_to_six_planewaves<multiplic,0,multiplic>(fmm_globs,fmm_level_globs,mpole,up,down,north,south,east,west);
}

#ifdef DELETED
template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void convert_and_add_accumulated_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                                            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& fmm_level_globs,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& up,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& down,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& north,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& south,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& east,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& west,
                                            MultiHolder<multiplic,MultipoleHolderT>& child_lexp)
{
    // workspace
    MultipoleHolderT local_out;
    MultipoleHolderT rotated_mpole;

    #pragma omp parallel for private(local_out, rotated_mpole)
    for (int m=begin; m < end; ++m)
    {
        // UP/DOWN
        pw_to_local(fmm_globs, fmm_level_globs, up[m], down[m], local_out);
        child_lexp[m] += local_out;

        // NORTH/SOUTH
        pw_to_local(fmm_globs, fmm_level_globs, north[m], south[m], local_out);
        rotate_ytoz(local_out, rotated_mpole, fmm_globs.rdpi2);
        child_lexp[m] += rotated_mpole;

        // EAST/WEST
        pw_to_local(fmm_globs, fmm_level_globs, east[m], west[m], local_out);
        rotate_ztox(local_out, rotated_mpole, fmm_globs.rdmpi2);
        child_lexp[m] += rotated_mpole;
    }

    return;
}
#else // ! DELETED
template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void convert_and_add_accumulated_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                                            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& fmm_level_globs,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& up,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& down,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& north,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& south,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& east,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& west,
                                            MultiHolder<multiplic,MultipoleHolderT>& child_lexp);
#endif // ! DELETED

template<typename ContentType, int n, int begin, int end, int NTERMS>
inline void yformmp(double beta,
            const Vector& x0y0z0,
            const std::vector<ContentType*>& charges,
            MultiHolder<n,BaseMultipoleHolder<NTERMS> > &mpoles,
            double scale,
            LegendreHolder<NTERMS>,
            const DblMatrix2D &c) // c initialised in frmini
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    this subroutine forms the multipole expansion cause by the nparts
    //      particles in the box.
    //
    //  on input:
    //
    //    beta : the frequency.
    //    x0y0z0: center of the expansion
    //    nparts: number of sources
    //    zparts(3,nparts): array of coordinates of sources
    //    charge(nparts): array of strengths of sources
    //    nparts: the total number of particles.
    //    NTERMS: order of desired expansion
    //    scale: the scaling factor.
    //
    //  on output:
    //
    //    mpole: coefficients of multipole expansion
    //
    //  working space :
    //
    //    p: used for storing the associate legendre polynomials.
    //
    //  subroutine called : dsqrt(), in(), lgndr()
    //  called from : brfrc()
    //
    //  note 1: this subroutine needs the precomputed variables c(,)
    //          derived from entry frmini()
    //
    //       2: the multipole expansion is scaled to avoid over- and
    //          under-flow.
    //
    //       3: only the n=0, ... ,NTERMS, m=0, ..., NTERMS coefficients
    //          are calculated.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

       boost::scoped_ptr< LegendreHolder<NTERMS> > p_ptr;
       boost::scoped_ptr<CmplxMatrix1D> ephi_ptr;
       BaseMultipoleHolder<NTERMS> mpole; // working copy

       for (int charge_ctr=0; charge_ctr < charges.size(); ++charge_ctr)
       {
            if (p_ptr.get() == NULL)
            {
                p_ptr.reset(new LegendreHolder<NTERMS>);
            }
            LegendreHolder<NTERMS>& p = *p_ptr;
            
            if (ephi_ptr.get() == NULL)
            {
                ephi_ptr.reset(new CmplxMatrix1D(Range(0,NTERMS)));
            }
            CmplxMatrix1D& ephi = *ephi_ptr;
            
            const ContentType& thing = *(charges[charge_ctr]);
            
            double d, costheta;
    
            // scope this bit off to avoid having too many local variables floating about
            {
                // get location of charge from centre of cube in spherical coordinates.

                Vector r = (thing - x0y0z0);

                double proj = r.x*r.x + r.y*r.y;
                double rr = proj + r.z*r.z;
                proj = sqrt(proj);
                d = sqrt(rr);

                if ( d <= PRECIS ) {
                    costheta = 1.0;
                }
                else {
                    costheta = r.z/d;
                }

                if ( proj <= PRECIS*d ) {
                    ephi(0) = std::complex<double>(1.0,0.0);
                }
                else {
                    ephi(0) = std::complex<double>(r.x/proj,-r.y/proj);
                }
            }

            for (int i=1; i < NTERMS; ++i)
            {
                ephi(i) = ephi(i-1)*ephi(0);
            }

            double rk=d*beta;
            DblMatrix1D bi(Range(0,NTERMS+1));
            i_n(scale,rk,NTERMS,bi);

    //
    //-------compute legendre polynomials of argument cos(theta) = ctheta
    //         and add contributions from legendre polynomials
    //
            lgndr(costheta, p);

            mpole.reset();
            mpole(0,0) = bi(0);
            for (int l=1; l <= NTERMS; ++l)
            {
                mpole(l,0) = p(l,0)*bi(l)*c(l,0);
                for (int m=1; m <= l; ++m)
                {
                    double cp = bi(l)*c(l,m)*p(l,m);
                    mpole(l,m) = cp*ephi(m-1);
                }
            }

            #pragma omp parallel for
            for (int nctr=begin; nctr < end; ++nctr)
            {
                double charge = thing.get_charge(nctr);
                {
                    BaseMultipoleHolder<NTERMS>& actual_mpole = mpoles[nctr];
                    for (int l=0; l <= NTERMS; ++l)
                    {
                        for (int m=0; m <= l; ++m)
                        {
                            {
                                actual_mpole(l,m) += mpole(l,m)*charge;
                            }
                        }
                    }
                }
            }
    }

    return;
}

template<typename ContentType, int n, int NTERMS>
inline void yformmp(double beta,
        const Vector& x0y0z0,
        const std::vector<ContentType*>& charges,
        MultiHolder<n,BaseMultipoleHolder<NTERMS> > &mpoles,
        double scale,
        LegendreHolder<NTERMS> &p,
        const DblMatrix2D &c)
{
    yformmp<ContentType,n,0,n,NTERMS>(beta,x0y0z0,charges,mpoles,scale,p,c);
}

template<typename ContentType, int NTERMS>
inline void yformmp(double beta,
            const Vector& x0y0z0,
            const std::vector<ContentType*>& charges,
            BaseMultipoleHolder<NTERMS> &mpoles,
            double scale,
            LegendreHolder<NTERMS> &p,
            const DblMatrix2D &c) // c initialised in frmini
{
    yformmp<ContentType,1,0,1,NTERMS>(beta,x0y0z0,charges,mpoles,scale,p,c);
}

template <int multiplic, int begin, int end, int NTERMS, int NLAMBS, int NWAVES>
inline void translate_up_wave(const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                    const MultiHolder<multiplic, BasePlaneWaveHolder<NLAMBS,NWAVES> >& waves,
                    MultiHolder<multiplic, BasePlaneWaveHolder<NLAMBS,NWAVES> >& targets,
                    short dx,
                    short dy,
                    short dz)
{
    unsigned short nexptotp = fmm_level_globs.nexptotp;
    const CmplxMatrix2D& xs = fmm_level_globs.xs;
    const CmplxMatrix2D& ys = fmm_level_globs.ys;
    const DblMatrix2D& zs = fmm_level_globs.zs;

    for (int multiplicity=begin; multiplicity < end; ++multiplicity)
    {
        const BasePlaneWaveHolder<NLAMBS,NWAVES>& wave = waves[multiplicity];
        BasePlaneWaveHolder<NLAMBS,NWAVES>& target = targets[multiplicity];

        assert(dz >= 0);
        
        #pragma omp for
        for (int i=0; i < nexptotp; ++i)
        {
            double tmp_real = 1.0;
            if (dz > 0)
            {
                tmp_real = zs(dz-1,i);
            }

            std::complex<double> tmp(tmp_real, 0.0);

            if (dy > 0)
            {
                tmp *= ys(dy-1,i);
            }
            else if (dy < 0)
            {
                tmp *= conj(ys((-dy)-1,i));
            }

            if (dx > 0)
            {
                tmp *= xs(dx-1,i);
            }
            else if (dx < 0)
            {
                tmp *= conj(xs((-dx)-1,i));
            }

            if (dz < 0)
            {
                tmp /= zs((-dz)-1,i);
            }

            target(i) += wave(i)*tmp;
        }
    }

    return;
}

template <int multiplic, int begin, int end, int NTERMS, int NLAMBS, int NWAVES>
inline void translate_down_wave(const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                        const MultiHolder<multiplic, BasePlaneWaveHolder<NLAMBS,NWAVES> >& waves,
                        MultiHolder<multiplic, BasePlaneWaveHolder<NLAMBS,NWAVES> >& targets,
                        short dx,
                        short dy,
                        short dz)
{
    unsigned short nexptotp = fmm_level_globs.nexptotp;
    const CmplxMatrix2D& xs = fmm_level_globs.xs;
    const CmplxMatrix2D& ys = fmm_level_globs.ys;
    const DblMatrix2D& zs = fmm_level_globs.zs;

    for (int multiplicity=begin; multiplicity < end; ++multiplicity)
    {
        const BasePlaneWaveHolder<NLAMBS,NWAVES>& wave = waves[multiplicity];
        BasePlaneWaveHolder<NLAMBS,NWAVES>& target = targets[multiplicity];

        assert(dz <= 0);
        #pragma omp for
        for (int i=0; i < nexptotp; ++i)
        {
            double tmp_real = 1.0;
            if (dz < 0)
            {
                tmp_real = zs((-dz)-1,i);
            }

            std::complex<double> tmp(tmp_real, 0.0);

            if (dy > 0)
            {
                tmp *= conj(ys(dy-1,i));
            }
            else if (dy < 0)
            {
                tmp *= ys((-dy)-1,i);
            }

            if (dx > 0)
            {
                tmp *= conj(xs(dx-1,i));
            }
            else if (dx < 0)
            {
                tmp *= xs((-dx)-1,i);
            }

            if (dz > 0)
            {
                tmp_real /= zs(dz-1,i);
            }

            target(i) += wave(i)*tmp;
        }
    }

    return;
}


template <int multiplic, int begin, int end, int NTERMS>
inline void ympshift(int ifl, // id within parent node
                    const MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &mpoles_in, // input multipoles
                    MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &mpoles_out, // output multipoles
                    BaseMultipoleHolder<NTERMS>&,
                    BaseMultipoleHolder<NTERMS>&,
                    const DblMatrix3D& dc, // precomputed shift coefficients
                    const DblMatrix3D& rd) // rotation matrix
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    this subroutine shifts the center of a child box multipole
    //      expansion to the parent center, via the rotation scheme.
    //      we rotate the coordinate system, shift along the z-axis,
    //      and then rotate back.
    //      there are eight possible child locations, defined in the
    //      calling sequence of this routine by the parameters ifl
    //      and rd (see below).
    //
    //  on input:
    //
    //     integer *4  ifl    = flag which describes the quadrant in which
    //                          the child box lies (1,2,3 or 4).
    //     complex *16 mpole  = coefficients of original multipole exp.
    //     real *8     dc     = precomputed array containing
    //                          the shifting coefficients
    //                          along the z-axis. this is precomputed by
    //                          the subroutine mpshftcoef() at the beginning
    //                          of different levels.
    //     real *8     rd     = precomputed array containing rotation matrix
    //                          about y-axis.
    //                 there are two possible y rotations, depending on
    //                 whether the child box lies in +z half space or the
    //                 -z half space. they are referred to in the calling
    //                 program as rdp and rdm, respectively.
    //                 this is precomputed in the subroutine rotgen().
    //
    //     complex *16 marray = work array
    //
    //  on output:
    //
    //     complex *16 mpolen = coefficients of shifted multipole exp.
    //
    //     note 1 : the rotation part is the same as the old subroutine
    //              of the laplace equation. the shifting along the z-axis
    //              is changed to the new version. the rotation matrix
    //              is precomputed at the very beginning and the shifting
    //              matrix is computed at different levels.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    CmplxMatrix1D ephi(Range(0,NTERMS+3));
    ephi(0) = 1.0;
    const double arg = sqrt(2.0)/2.0;
    if (ifl == 4 || ifl == 8) {
        ephi(1) = std::complex<double>(-arg,arg);
    }
    else if (ifl == 3 || ifl == 7) {
        ephi(1) = std::complex<double>(arg,arg);
    }
    else if (ifl == 1 || ifl == 5) {
        ephi(1) = std::complex<double>(arg,-arg);
    }
    else if (ifl == 2 || ifl ==6) {
        ephi(1) = std::complex<double>(-arg,-arg);
    }
    else {
        std::cerr << "Bad octant-id passed into ympshift function." << std::endl;
        throw;
    }

    //
    //-----create array of powers of e^(i*m*phi).
    //
    for (int ell=1; ell <= NTERMS+1; ++ell) {
        ephi(ell+1) = ephi(ell)*ephi(1);
    }

    BaseMultipoleHolder<NTERMS> mwrk1;
    BaseMultipoleHolder<NTERMS> mwrk2;

    #pragma omp parallel for private(mwrk1, mwrk2)
    for (int multiplicity=begin; multiplicity < end; ++multiplicity)
    {
        const BaseMultipoleHolder<NTERMS>& mpole_in = mpoles_in[multiplicity];
        BaseMultipoleHolder<NTERMS>& mpole_out = mpoles_out[multiplicity];

    //-----a rotation of phi radians about the z-axis in the
    //       original coordinate system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk1(ell,m)=conj(ephi(m))*mpole_in(ell,m);
            }
        }
    //
    //-----a rotation about the y'-axis  in the rotated system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk2(ell,m)=mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; ++mp) {
                    mwrk2(ell,m)=mwrk2(ell,m)+mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m);
                }
            }
        }
    //
    //-----shift along z-axis.
    //       note that everything is scaled.
    //
        for (int n=0; n <= NTERMS; ++n) {
            for (int m=0; m <= n; ++m) {
                mwrk1(n,m)=0.0;
                for (int n_dash=m; n_dash <= NTERMS; ++n_dash) {
                    mwrk1(n,m) += mwrk2(n_dash,m)*dc(m,n,n_dash);
                }
            }
        }

    //
    //-----reverse rotation about the y'-axis.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {

            for (int m=0; m<=ell; m+=2) {
                mwrk2(ell,m)=mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)-(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
                for (int mp=2; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)+(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
            }

            for (int m=1; m<=ell; m+=2) {
                mwrk2(ell,m) = -mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)+(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
                for (int mp=2; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)-(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
            }
        }

    //
    //-----rotate back phi radians about the z-axis in the above system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk1(ell,m)=ephi(m)*mwrk2(ell,m);
            }
        }

        // store result
        mpole_out += mwrk1;
    }

    return;

}

template <int multiplic, int NTERMS>
inline void ympshift(int ifl, // id within parent node
                const MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &mpoles_in, // input multipoles
                MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &mpoles_out, // output multipoles
                BaseMultipoleHolder<NTERMS>& mwrk1,
                BaseMultipoleHolder<NTERMS>& mwrk2,
                const DblMatrix3D& dc, // precomputed shift coefficients
                const DblMatrix3D& rd) // rotation matrix
{
    ympshift<multiplic,0,multiplic>(ifl, mpoles_in, mpoles_out, mwrk1, mwrk2, dc, rd);
}

template <int multiplic, int begin, int end, int NTERMS>
inline void ylcshift(int ifl, // id within parent node
                const MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &locals_in, // input multipoles
                MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &locals_out, // output multipoles
                BaseMultipoleHolder<NTERMS>& mwrk1,
                BaseMultipoleHolder<NTERMS>& mwrk2,
                const DblMatrix3D& dc, // precomputed shift coefficients
                const DblMatrix3D& rd) // rotation matrix
{

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //     this subroutine shifts the local expansion of a parent cell
    //     to the center of one of its children, via the rotation scheme.
    //     that is, we rotate the coordinate system, shift along the z-axis,
    //     and then rotate back.
    //     there are eight possible child locations, defined in the
    //     calling sequence of this routine by the parameters ifl
    //     and rd (see below).
    //
    // on input:
    //
    //     integer *4  ifl    = flag which describes the quadrant in which
    //                          the child box lies (1,2,3 or 4).
    //     complex *16 local  = coefficients of original multipole exp.
    //     integer *4  NTERMS = integer indicates the terms retained in the
    //                          expansion.
    //     real *8     dc     = precomputed array containing coefficients
    //                          for the local translation along the z axis.
    //                          this is precomputed by the subroutine
    //                          lcshftcoef()
    //     real *8     rd     = precomputed array containing rotation matrix
    //                          about y-axis.
    //                 there are two possible y rotations, depending on
    //                 whether the child box lies in +z half space or the
    //                 -z half space. they are referred to in the calling
    //                 program as rdp and rdm, respectively.
    //                          this is precomputed by the subroutine
    //                          rotgen().
    //     complex *16 marray = work array
    //
    // on output:
    //
    //     complex *16 localn = coefficients of shifted multipole exp.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    CmplxMatrix1D ephi(Range(0,NTERMS+3));
    ephi(0) = 1.0;
    double arg = sqrt(2.0)/2.0;

    if (ifl == 4 || ifl == 8) {
        ephi(1) = std::complex<double>(arg,-arg);
    }
    else if (ifl == 3 || ifl == 7) {
        ephi(1) = std::complex<double>(-arg,-arg);
    }
    else if (ifl == 1 || ifl == 5) {
        ephi(1) = std::complex<double>(-arg,arg);
    }
    else if (ifl == 2 || ifl == 6) {
        ephi(1) = std::complex<double>(arg,arg);
    }
    else {
        std::cerr << "Bad octant-id passed into ympshift function." << std::endl;
        throw;
    }

//
//-----create array of powers of e^(i*m*phi).
//
    for (int ell=1; ell <= NTERMS+1; ++ell) {
        ephi(ell+1) = ephi(ell)*ephi(1);
    }

    for (int multiplicity=begin; multiplicity < end; ++multiplicity)
        {
            const BaseMultipoleHolder<NTERMS>& local_in = locals_in[multiplicity];
            BaseMultipoleHolder<NTERMS>& local_out = locals_out[multiplicity];

    //
    //-----a rotation of phi radians about the z-axis in the
    //       original coordinate system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk1(ell,m)=conj(ephi(m))*local_in(ell,m);
            }
        }

    //
    //-----a rotation about the y'-axis  in the rotated system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk2(ell,m)=mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; ++mp) {
                    mwrk2(ell,m)=mwrk2(ell,m)+mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m);
                }
            }
        }

    //
    //-----shift along z-axis.
    //       note that everything is scaled.
    //
        for (int n=0; n <= NTERMS; ++n) {
            for (int m=0; m <= n; ++m) {
                mwrk1(n,m)=0.0;
                for (int n_dash=m; n_dash <= NTERMS; ++n_dash) {
                    mwrk1(n,m)=mwrk1(n,m)+mwrk2(n_dash,m)*dc(m,n,n_dash);
                }
            }
        }

    //
    //-----reverse rotation about the y'-axis.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {

            for (int m=0; m<=ell; m+=2) {
                mwrk2(ell,m)=mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)-(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
                for (int mp=2; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)+(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
            }

            for (int m=1; m<=ell; m+=2) {
                mwrk2(ell,m) = -mwrk1(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)+(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
                for (int mp=2; mp <= ell; mp += 2) {
                    mwrk2(ell,m)=mwrk2(ell,m)-(mwrk1(ell,mp)*rd(ell,mp,m)+conj(mwrk1(ell,mp))*rd(ell,mp,-m));
                }
            }
        }

    //
    //-----rotate back phi radians about the z-axis in the above system.
    //
        for (int ell=0; ell <= NTERMS; ++ell) {
            for (int m=0; m <= ell; ++m) {
                mwrk1(ell,m)=ephi(m)*mwrk2(ell,m);
            }
        }

        local_out += mwrk1;
    }

    return;

}

template <int multiplic, int NTERMS>
inline void ylcshift(int ifl, // id within parent node
                const MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &locals_in, // input multipoles
                MultiHolder<multiplic, BaseMultipoleHolder<NTERMS> > &locals_out, // output multipoles
                BaseMultipoleHolder<NTERMS>& mwrk1,
                BaseMultipoleHolder<NTERMS>& mwrk2,
                const DblMatrix3D& dc, // precomputed shift coefficients
                const DblMatrix3D& rd) // rotation matrix
{
    ylcshift<multiplic,0,multiplic>(ifl, locals_in, locals_out, mwrk1, mwrk2, dc, rd);
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void convert_mp_to_exp(const FMM_Globals<NTERMS>& fmm_globs,
                        const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                        const BaseMultipoleHolder<NTERMS>& mpole_in,     // the multipole expansions (input)
                        BasePlaneWaveHolder<NLAMBS,NWAVES>& mexpuphys,           // +ve plane wave (output)
                        BasePlaneWaveHolder<NLAMBS,NWAVES>& mexpdphys            // -ve plane wave (output)
                        )
{
    const DblMatrix2D& scale_factors = fmm_globs.scale_factors; // square roots of binomial coeffs

    unsigned short nexptot = fmm_level_globs.nexptot;
    unsigned short nexptotp = fmm_level_globs.nexptotp;
    const CmplxMatrix1D& fexpe = fmm_level_globs.fexpe;
    const CmplxMatrix1D& fexpo = fmm_level_globs.fexpo;
    const unsigned int fsize = fmm_level_globs.fsize;
    const DblMatrix3D& rlsc = fmm_level_globs.rlsc;      // scaled legendre functions
    const UShrtMatrix1D &numfour = fmm_level_globs.numfour;
    const UShrtMatrix1D &numphys = fmm_level_globs.numphys;
    const DblMatrix1D& quad_nodes = fmm_level_globs.quad_nodes;

    CmplxMatrix1D mexpup(Range(0,nexptot));
    CmplxMatrix1D mexpdn(Range(0,nexptot));

    // convert multipole to exponential (Fourier)
    ympoletoexp<NLAMBS>(nexptot, numfour, rlsc, mpole_in, mexpup, mexpdn);

    // up and down components (Fourier to phys)
    yftophys(nexptot, nexptotp, fsize, quad_nodes, numfour, numphys, fexpe, fexpo, mexpup, mexpuphys);
    yftophys(nexptot, nexptotp, fsize, quad_nodes, numfour, numphys, fexpe, fexpo, mexpdn, mexpdphys);

    return;
}

template<int NTERMS, int NLAMBS, int NWAVES>
inline void pw_to_local(const FMM_Globals<NTERMS>& fmm_globs,
                const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                const BasePlaneWaveHolder<NLAMBS,NWAVES>& wave_up ,      // +ve plane wave (input)
                const BasePlaneWaveHolder<NLAMBS,NWAVES>& wave_down,    // -ve plane wave (input)
                BaseMultipoleHolder<NTERMS>& local_out                          // the local expansions
                )
{

    double betascal = fmm_level_globs.betascal;
    unsigned short nexptot = fmm_level_globs.nexptot;
    unsigned short nexptotp = fmm_level_globs.nexptotp;
    int fsize = fmm_level_globs.fsize;
    const DblMatrix2D& scale_factors = fmm_globs.scale_factors;
    const DblMatrix3D& rlsc = fmm_level_globs.rlsc;
    const UShrtMatrix1D &numfour = fmm_level_globs.numfour;
    const UShrtMatrix1D &numphys = fmm_level_globs.numphys;
    const DblMatrix1D& quad_nodes = fmm_level_globs.quad_nodes;
    const DblMatrix1D& quad_weights = fmm_level_globs.quad_weights;
    const CmplxMatrix1D& fexpback = fmm_level_globs.fexpback;

    CmplxMatrix1D fmode_up(Range(0,nexptot));
    CmplxMatrix1D fmode_down(Range(0,nexptot));

    // convert physical waves to fourier representation
    yphystof(nexptot,
            nexptotp,
            fsize,
            quad_nodes,
            numfour,
            numphys,
            fexpback,
            wave_up, // input
            fmode_up // output
            );

    yphystof(nexptot,
            nexptotp,
            fsize,
            quad_nodes,
            numfour,
            numphys,
            fexpback,
            wave_down, // input
            fmode_down // output
            );

    yexptolocal<NLAMBS>(nexptot,
                betascal,
                rlsc,
                local_out,  // output
                quad_nodes,
                quad_weights,
                scale_factors,
                numfour,
                fmode_up,  // input
                fmode_down); // input

    return;
}

template<int NTERMS>
inline void rotate_ztoy(const BaseMultipoleHolder<NTERMS>& mpole, // input unrotated multipoles
                BaseMultipoleHolder<NTERMS>& mrotate,     // output
                const DblMatrix3D& rdminus   // rotation matrix
                )
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    the rotation matrix used in the subroutine mknsexp()
    //    so the north_south expansions are made the same way
    //    as the up-down expansions.
    //
    //  on input :
    //
    //    NTERMS : number of terms in the multipole expansion.
    //    mpole : the multipole expansion.
    //    rdminus : the rotation matrix generated in subroutine
    //      rotgen<- fstrtn()
    //
    //  output :
    //    mrotate : the rotated multiple expansion coefficients.
    //
    //  working space : mwork().
    //
    //  called from :  mknsexp()
    //
    //  subroutine called : none.
    //
    //  note : this can be further simplified?
    //
    //     end result         z_new <- y_old
    //                        y_new <- x_old
    //                        x_new <- z_old
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    BaseMultipoleHolder<NTERMS> mwork; // workspace

    //
    //-----a rotation of -pi/2 radians about the z-axis in the
    //       original coordinate system.
    //
    CmplxMatrix1D ephi(Range(0,NTERMS+1));
    ephi(0) = std::complex<double>(1.0,0.0);
    for (int m=1; m <= NTERMS; ++m)
    {
        ephi(m)=ephi(m-1)*std::complex<double>(0.0,-1.0);
    }

        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
                mwork(ell,m)=ephi(m)*mpole(ell,m);
            }
        }

        //
        //-----a rotation of -pi/2 radians about the y'-axis in the
        //       new coordinate system, bringing the +z-axis in line
        //       with the original +y axis.
        //
        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
                mrotate(ell,m)=mwork(ell,0)*rdminus(ell,0,m);
                for (int mp=1; mp <= ell; ++mp)
                {
                    mrotate(ell,m) += mwork(ell,mp)*rdminus(ell,mp,m) + conj(mwork(ell,mp))*rdminus(ell,mp,-m);
                }
            }
        }

    return;

}

template<int NTERMS>
inline void rotate_ztox(const BaseMultipoleHolder<NTERMS>& mpole, // input unrotated multipoles
                BaseMultipoleHolder<NTERMS>& mrotate,     // output
                const DblMatrix3D& rd)
{

    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    the rotation matrix used in the subroutine ladapfmm(), yadapfmm(),
    //    and mkewexp(). rotate the east-west expansions from or to the
    //    normal direction, depending on rd.
    //
    //  on input :
    //
    //    NTERMS : number of terms in the multipole/local expansion.
    //    mpole : the multipole/local expansion.
    //    rd : the rotation matrix generated in subroutine
    //         rotgen<- fstrtn(), it can either be rdplus or rdminus
    //
    //  on output :
    //    mrotate : the rotated local expansion coefficients.
    //
    //  called from :  mkewexp(), ladapfmm(), yadapfmm()
    //
    //  subroutine called : none.
    //
    //  note : this can be further simplified?
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
                mrotate(ell,m) = mpole(ell,0)*rd(ell,0,m);
                for (int mp=1; mp <= ell; ++mp)
                {
                    mrotate(ell,m) += mpole(ell,mp)*rd(ell,mp,m) + conj(mpole(ell,mp))*rd(ell,mp,-m);
                }
            }
        }

    return;
}

template<int NTERMS>
inline void rotate_ytoz(const BaseMultipoleHolder<NTERMS>& mpole, // input unrotated multipoles
                BaseMultipoleHolder<NTERMS>& mrotate,     // output
                const DblMatrix3D& rdplus   // rotation matrix
                )
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    the rotation matrix used in the subroutine ladapfmm(), yadapfmm()
    //      rotate the totated north_south expansions back to the
    //      normal direction.
    //
    //  on input :
    //    NTERMS : number of terms in the multipole expansion.
    //    mpole : the local expansion.
    //    rdplus : the rotation matrix generated in subroutine
    //      rotgen<- fstrtn()
    //
    //  on output :
    //    mrotate : the rotated local expansion coefficients.
    //
    //  working space : mwork().
    //
    //  called from :  ladapfmm(), yadapfmm()
    //  subroutine called : none.
    //
    //  note : this can be further simplified?
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    BaseMultipoleHolder<NTERMS> mwork;



        //
        //----a rotation of pi/2 radians about the y'-axis in the
        //      new coordinate system, bringing the +z-axis in line
        //      with the original +z axis.
        //
        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
                mwork(ell,m)=mpole(ell,0)*rdplus(ell,0,m);
                for (int mp=1; mp <= ell; ++mp)
                {
                    mwork(ell,m) += mpole(ell,mp)*rdplus(ell,mp,m) + conj(mpole(ell,mp))*rdplus(ell,mp,-m);
                }
            }
        }

        //
        //-----a rotation of pi/2 radians about the z-axis in the
        //       original coordinate system.
        //
        CmplxMatrix1D ephi(Range(0,NTERMS+1));
        ephi(0) = std::complex<double>(1.0,0.0);
        for (int m=1; m <= NTERMS; ++m)
        {
            ephi(m)=ephi(m-1)*std::complex<double>(0.0,1.0);
        }

        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
                mrotate(ell,m)=ephi(m)*mwork(ell,m);
            }
        }

    return;
}

}// end namespace

#ifndef __DELETED__
template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void fmm::convert_mpole_to_six_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                const Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& fmm_level_globs,
                const MultiHolder<multiplic,MultipoleHolderT>& mpole,
                MultiHolder<multiplic,PlaneWaveHolderT>& up,
                MultiHolder<multiplic,PlaneWaveHolderT>& down,
                MultiHolder<multiplic,PlaneWaveHolderT>& north,
                MultiHolder<multiplic,PlaneWaveHolderT>& south,
                MultiHolder<multiplic,PlaneWaveHolderT>& east,
                MultiHolder<multiplic,PlaneWaveHolderT>& west)
{

    MultipoleHolderT mpole_rotation;
    #pragma omp parallel for private(mpole_rotation)
    for (int m=begin; m < end; ++m)
    {
        // convert multipole expansion to planewave expansions
        convert_mp_to_exp(fmm_globs, fmm_level_globs, mpole[m], up[m], down[m]);

        rotate_ztoy(mpole[m], mpole_rotation, fmm_globs.rdmpi2);
        convert_mp_to_exp(fmm_globs, fmm_level_globs, mpole_rotation, north[m], south[m]);

        rotate_ztox(mpole[m], mpole_rotation, fmm_globs.rdpi2);
        convert_mp_to_exp( fmm_globs, fmm_level_globs, mpole_rotation, east[m], west[m]);
    }
    return;

}

template <int multiplic, int begin, int end, typename MultipoleHolderT, typename PlaneWaveHolderT, int NTERMS, int NLAMBS>
inline void fmm::convert_and_add_accumulated_planewaves(const FMM_Globals<NTERMS>& fmm_globs,
                                            const Level_Dependent_FMM_Globals<NTERMS,NLAMBS>& fmm_level_globs,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& up,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& down,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& north,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& south,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& east,
                                            const MultiHolder<multiplic,PlaneWaveHolderT>& west,
                                            MultiHolder<multiplic,MultipoleHolderT>& child_lexp)
{
    // workspace
    MultipoleHolderT local_out;
    MultipoleHolderT rotated_mpole;

    #pragma omp parallel for private(local_out, rotated_mpole)
    for (int m=begin; m < end; ++m)
    {
        // UP/DOWN
        pw_to_local(fmm_globs, fmm_level_globs, up[m], down[m], local_out);
        child_lexp[m] += local_out;

        // NORTH/SOUTH
        pw_to_local(fmm_globs, fmm_level_globs, north[m], south[m], local_out);
        rotate_ytoz(local_out, rotated_mpole, fmm_globs.rdpi2);
        child_lexp[m] += rotated_mpole;

        // EAST/WEST
        pw_to_local(fmm_globs, fmm_level_globs, east[m], west[m], local_out);
        rotate_ztox(local_out, rotated_mpole, fmm_globs.rdmpi2);
        child_lexp[m] += rotated_mpole;
    }

    return;
}

#endif // ! __DELETED__

#endif /* FMM_H_ */
