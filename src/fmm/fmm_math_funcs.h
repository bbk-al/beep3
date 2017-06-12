/*
 * fmm_math_funcs.h
 *
 *  Created on: 21 Oct 2010
 *      Author: david
 */

#ifndef FMM_MATH_FUNCS_H_
#define FMM_MATH_FUNCS_H_

// Bessel and Gamma functions ported by John Burkardt
// from the original gangster fortran to C++.
#include "bessel_gamma.h"

namespace fmm
{

// some useful mathsy-ish constants
const std::complex<double> ima(0,1);
const double halfpi=atan(1.0)*2.0;
const double pi=atan(1.0)*4.0;
const double TWO_OVER_PI = 2.0 / pi;

inline void i_n(double scal, double x, unsigned int nb, DblMatrix1D &b)
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    calculates the i_n(z)=sqrt(pi/2/z)*i_(n+1/2,z).
    //
    //  on input:
    //    scal: the scaling factor to avoid underflow.
    //    x: the parameter for i_n(x)
    //    nb: number of terms for the subindex n, n=0:nb.
    //
    //  on output :
    //    b: contains the i_n(x), n=0:nb.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    //integer *4 i,ize
        //real *8 const,ensig,enmten,halfpi
        //real *8 xscal,term1,term2,alpha
        //data ensig, enmten/1.0d-4, 1.0d-300/
//
//-----function called
//

    // x gotta be >0
    if (x < 0) {
        std::cerr << "Bad x argument to i_n function" << std::endl;
        throw;
    }

    // TODO: use <limits> for these nasty looking variables
    static const double ensig=1e-4, enmten=1e-300;
    const double inv_scal = 1.0 / scal;
//
//-----if x.le. 10e-4, then use the 2 terms taylor expansion.
//

    if (x <= ensig) {
        double xscal = x * inv_scal;
        double term1=1.0;
        double term2=0.5*x*x;
        b(0) = term1*(1.0 + term2 / 3.0);
        for(unsigned int i=1; i <= nb; ++i)
        {
            term1 *= xscal/double(2*i+1);
            if (term1 <= enmten)
            {
                term1 = 0.0;
            }
            b(i)=term1*(1.0 + term2/double(2*i+3) );
        }

//
//-------usually, scal should be less than one, however, in the
//         calculation of h_n, scal will be greater than one,
//         in this case, we should modify the program and replace every
//         small item with zero.
//
    }
    else if (x > 100.0) {
        for(unsigned int i=0; i<=nb; ++i)
        {
            b(i) = 0.0;
        }
    }
    else {
//
//-------otherwise, we will call besselj() and then scale the
//         jn by the scaling factor.
//

        double constant = sqrt(halfpi/x);
        int ncalc = ::ribesl(x, 0.5, nb+1, 1, b.data());
        
        for(unsigned int i=0; i <= nb; ++i)
        {
            
            b(i) *= constant;
            
            constant *= inv_scal;
            if ( fabs(b(i)) <= enmten ) {
                constant = 0.0;
            }
        }
    }

    return;
}

template<unsigned short nmax>
void lgndr(double x, // the argument
			TriangularMemory<nmax,double>& y // where to put the legendre associated polynomials
          )
{
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
    //
    //  purpose:
    //
    //    this subroutine computes the legendre polynomial expansion using
    //    a recursive expansion.
    //
    //  on input:
    //    nmax: the max number of terms in the expansion.
    //    x: where we want to evaluate the expansion.
    //
    //  on output:
    //    y: the function value at x.
    //
    //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    // TODO Throw exceptions on math errors (e.g. domain error)

    double u = - sqrt(1.0 - x*x);
    y(0,0) = 1.0;
    for (unsigned short m=0; m <=nmax; ++m)
    {
        if (m > 0) {
            y(m,m) = y(m-1,m-1) * u * double(2*m-1);
        }
        if (m < nmax) {
            y(m+1,m) = double(2*m+1)*x*y(m,m);
        }
        for (unsigned short n=m+2; n <= nmax; ++n)
        {
            y(n,m) = ((2.0*double(n)-1.0)*x*y(n-1,m)-double(n+m-1)*y(n-2,m)) / double(n-m);
        }
    }
    return;
}

template<unsigned short nmax>
void lgndrgt1(double scal, double x, TriangularMemory<nmax,double>& y)
{
    //**********************************************************************
    //
    //-----this subroutine calculates the legendre function for x >1.
    //     the scaled version. scal^n* p_n^m
    //     scal should be less than 1.
    //
    //     input :
    //       scal : the scaling factor.
    //       nmax : the number of terms we want to compute.
    //       x : the value at which we want to compute.
    //
    //     output :
    //       y() : the scal^n* p_n^m.
    //
    //***********************************************************************
/*      implicit real *8 (a-h,o-z)
      integer *4 nmax, m
      real *8  scal, x, y(0:nmax,0:nmax), u, v, w*/
//

    // fairly arbitary sanity test
    assert(nmax < 50);

    double v = scal*x;
    double w = scal*scal;
    double u = sqrt(x*x - 1.0) * scal;
    y(0,0) = 1.0;
    for (int m=0; m <= nmax; ++m)
    {
        if (m > 0) { y(m,m) = y(m-1,m-1)*u*double(2*m-1); }
        if (m < nmax) { y(m+1,m)=double(2*m+1)*v*y(m,m); }
        for (int n=m+2; n <= nmax; ++n)
        {
            y(n,m)= ((2.0*double(n)-1.0)*v*y(n-1,m)-double(n+m-1)*w*y(n-2,m)) / double(n-m);
        }
    }

    for (int n=0; n <=nmax; ++n)
    {
        for (int m=0; m <= n; ++m)
        {
            //double sgn = (abs(m) % 2) ? -1.0 : 1.0;
            double sgn = 1.0;
            y(n,m) *= sgn;
        }
    }

    return;
}

} // end fmm namespace

#endif /* FMM_MATH_FUNCS_H_ */
