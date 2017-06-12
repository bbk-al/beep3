#ifndef BESSEL_GAMMA_H_
#define BESSEL_GAMMA_H_

//
// This code is taken from:
//
// http://people.sc.fsu.edu/~jburkardt/cpp_src/prob/prob.C
//
// Which is a GNU LGPL licensed bunch of very very useful maths functions.
// Thanks John Burkardt for porting these functions from the original
// fortran codes.

# include <cmath>

#define MIN(i1,i2) ((i1 < i2) ? i1 : i2)
#define MAX(r1,r2) ((r1 > r2) ? r1 : r2)

//****************************************************************************80
inline double r8_gamma ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03,
    8.4171387781295E-04,
   -5.952379913043012E-04,
    7.93650793500350248E-04,
   -2.777777777777681622553E-03,
    8.333333333333333331554247E-02,
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01,
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02,
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04,
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02,
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03,
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03,
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}

//****************************************************************************80
static int ribesl ( double x, double alpha, int nb, int ize, double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    RIBESL calculates I Bessel function with non-integer orders.
//
//  Discussion:
//
//    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
//    for non-negative argument X, and non-negative order N+ALPHA,
//    with or without exponential scaling.
//
//    This program is based on a program written by David
//    Sookne that computes values of the Bessel functions J or
//    I of real argument and integer order.  Modifications include
//    the restriction of the computation to the I Bessel function
//    of non-negative real argument, the extension of the computation
//    to arbitrary positive order, the inclusion of optional
//    exponential scaling, and the elimination of most underflow.
//
//    In case of an error, NCALC will not equal NB, and not all I's are
//    calculated to the desired accuracy.
//
//    If NCALC < 0:  An argument is out of range. For example,
//    NB <= 0, IZE is not 1 or 2, or IZE = 1 and EXPARG <= ABS(X)
//    In this case, the B-vector is not calculated, and NCALC is
//    set to MIN(NB,0)-1 so that NCALC /= NB.
//
//    If 0 < NCALC < NB, then not all requested function values could
//    be calculated accurately.  This usually occurs because NB is
//    much larger than ABS(X).  In this case, B(N) is calculated
//    to the desired accuracy for N <= NCALC, but precision
//    is lost for NCALC < N <= NB.  If B(N) does not vanish
//    for NCALC < N (because it is too small to be represented),
//    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
//    significant figures of B(N) can be trusted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Frank Olver, David Sookne,
//    A Note on Backward Recurrence Algorithms,
//    Mathematics of Computation,
//    Volume 26, 1972, pages 941-947.
//
//    David Sookne,
//    Bessel Functions of Real Argument and Integer Order,
//    NBS Journal of Research B,
//    Volume 77B, 1973, pages 125-132.
//
//    William Cody,
//    Algorithm 597:
//    Sequence of Modified Bessel Functions of the First Kind,
//    ACM Transactions of Mathematical Software,
//    Volume 9, Number 2, June 1983, pages 242-245.
//
//  Parameters:
//
//    Input, double X, the argument for which the functions
//    are to be calculated.
//
//    Input, double ALPHA,the fractional part of the order
//    for which the functions are to be calculated.
//    0 <= ALPHA < 1.0.
//
//    Input, int NB, the number of functions to be calculated.
//    The first function calculated is of order ALPHA, and the
//    last is of order (NB - 1 + ALPHA).  1 <= NB.
//
//    Input, int IZE, scaling option.
//    1, unscaled I's are to calculated,
//    2, exponentially scaled I's are to be calculated.
//
//    Output, double B[NB], the values of the functions
//    I(ALPHA,X) through I(NB-1+ALPHA,X), with scaling if requested.
//
//    Output, int RIBESL, the value of NCALC, the error indicator.
//    If NCALC = NB, then all the requested values were calculated
//    to the desired accuracy.
//
//  Local Parameeters:
//
//    BETA, the radix for the floating-point system.
//
//    MINEXP, smallest representable power of BETA.
//
//    MAXEXP, smallest power of BETA that overflows
//
//    IT, number of bits in the mantissa of a working precision variable.
//
//    NSIG, decimal significance desired.  Should be set to
//    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
//    in decreased accuracy while setting NSIG higher will
//    increase CPU time without increasing accuracy.  The
//    truncation error is limited to a relative error of
//    T=.5*10^(-NSIG).
//
//    ENTEN, 10.0^K, where K is the largest integer such that
//    ENTEN is machine-representable in working precision
//
//    ENSIG, 10.0^NSIG
//
//    RTNSIG, 10.0^(-K) for the smallest integer K such that
//    NSIG/4 <= K.
//
//    ENMTEN, smallest ABS(X) such that X/4 does not underflow
//
//    XLARGE, upper limit on the magnitude of X when IZE=2.  Bear
//    in mind that if ABS(X)=N, then at least N iterations
//    of the backward recursion will be executed.  The value
//    of 10.0^4 is used on every machine.
//
//    EXPARG, largest working precision argument that the library
//    EXP routine can handle and upper limit on the
//    magnitude of X when IZE=1; approximately log(BETA^MAXEXP).
//
//    Approximate values for some important machines are:
//
//                        beta       minexp      maxexp       it
//
//  CRAY-1        (S.P.)    2        -8193        8191        48
//  Cyber 180/855
//    under NOS   (S.P.)    2         -975        1070        48
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)    2         -126         128        24
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)    2        -1022        1024        53
//  IBM 3033      (D.P.)   16          -65          63        14
//  VAX           (S.P.)    2         -128         127        24
//  VAX D-Format  (D.P.)    2         -128         127        56
//  VAX G-Format  (D.P.)    2        -1024        1023        53
//
//
//                        NSIG       ENTEN       ENSIG      RTNSIG
//
// CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
// Cyber 180/855
//   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
// IEEE (IBM/XT,
//   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
// IEEE (IBM/XT,
//   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
// IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
// VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
// VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
// VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
//
//
//                         ENMTEN      XLARGE   EXPARG
//
// CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
// Cyber 180/855
//   under NOS   (S.P.)   1.25E-293    1.0E+4     741
// IEEE (IBM/XT,
//   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
// IEEE (IBM/XT,
//   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
// IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
// VAX           (S.P.)   1.17E-38     1.0E+4      88
// VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
// VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
//
{
  double constant = 1.585;
  double em;
  double empal;
  double emp2al;
  double en;
  double enmten = 8.9E-308;
  double ensig = 1.0E+16;
  double enten = 1.0E+308;
  double exparg = 709.0;
  bool flag;
  double half = 0.5;
  double halfx;
  int i;
  int k;
  int l;
  int magx;
  int n;
  int nbmx;
  int ncalc;
  int nend;
  int nsig = 16;
  int nstart;
  double one = 1.0;
  double p;
  double plast;
  double pold;
  double psave;
  double psavel;
  double rtnsig = 1.0E-04;
  double tempa;
  double tempb;
  double tempc;
  double test;
  double total;
  double tover;
  double two = 2.0;
  double xlarge = 1.0E+04;
  double zero = 0.0;
//
//  Check for X, NB, OR IZE out of range.
//
  if ( nb <= 0 )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( x < 0.0 )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( alpha < 0.0 )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( 1.0 <= alpha )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( ize == 1 && exparg < x )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( ize == 2 && xlarge < x )
  {
    ncalc = MIN ( nb, 0 ) - 1;
    return ncalc;
  }
//
//  Use 2-term ascending series for small X.
//
  ncalc = nb;
  magx = ( int ) ( x );
//
//  Initialize the forward sweep, the P-sequence of Olver.
//
  if ( rtnsig <= x )
  {
    nbmx = nb - magx;
    n = magx + 1;
    en = ( double ) ( n + n ) + ( alpha + alpha );
    plast = one;
    p = en / x;
//
//  Calculate general significance test.
//
    test = ensig + ensig;

    if ( 5 * nsig < 2 * magx )
    {
      test = sqrt ( test * p );
    }
    else
    {
      test = test / pow ( constant, magx );
    }
//
//  Calculate P-sequence until N = NB-1.  Check for possible overflow.
//
    flag = false;

    if ( 3 <= nbmx )
    {
      tover = enten / ensig;
      nstart = magx + 2;
      nend = nb - 1;

      for ( k = nstart; k <= nend; k++ )
      {
        n = k;
        en = en + two;
        pold = plast;
        plast = p;
        p = en * plast / x + pold;
//
//  To avoid overflow, divide P-sequence by TOVER.  Calculate
//  P-sequence until 1 < ABS(P).
//
        if ( tover < p )
        {
          tover = enten;
          p = p / tover;
          plast = plast / tover;
          psave = p;
          psavel = plast;
          nstart = n + 1;

          for ( ; ; )
          {
            n = n + 1;
            en = en + two;
            pold = plast;
            plast = p;
            p = en * plast / x + pold;

            if ( 1.0 < p )
            {
              break;
            }
          }

          tempb = en / x;
//
//  Calculate backward test, and find NCALC, the highest N
//  such that the test is passed.
//
          test = pold * plast / ensig;
          test = test * ( half - half / ( tempb * tempb ) );
          p = plast * tover;
          n = n - 1;
          en = en - two;
          nend = MIN ( nb, n );

          ncalc = nend + 1;

          for ( l = nstart; l <= nend; l++ )
          {
            pold = psavel;
            psavel = psave;
            psave = en * psavel / x + pold;

            if ( test < psave * psavel )
            {
              ncalc = l;
              break;
            }
          }
          ncalc = ncalc - 1;
          flag = true;
          break;
        }
      }

      if ( !flag )
      {
        n = nend;
        en = ( double ) ( n + n ) + ( alpha + alpha );
//
//  Calculate special significance test for 2 < NBMX.
//
        test = MAX ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) );
      }
    }
//
//  Calculate P-sequence until significance test passed.
//
    if ( !flag )
    {
      for ( ; ; )
      {
        n = n + 1;
        en = en + two;
        pold = plast;
        plast = p;
        p = en * plast / x + pold;

        if ( test <= p )
        {
          break;
        }
      }
    }
//
//  Initialize the backward recursion and the normalization sum.
//
    n = n + 1;
    en = en + two;
    tempb = zero;
    tempa = one / p;
    em = ( double ) ( n ) - one;
    empal = em + alpha;
    emp2al = ( em - one ) + ( alpha + alpha );
    total = tempa * empal * emp2al / em;
    nend = n - nb;
//
//  N < NB, so store B(N) and set higher orders to zero.
//
    if ( nend < 0 )
    {
      b[n-1] = tempa;
      nend = -nend;

      for ( l = 1; l <= nend; l++ )
      {
        b[n+l-1] = zero;
      }

      nend = n - 2;
//
//  Calculate via difference equation and store B(N), until N = 2.
//
      if ( 0 < nend )
      {
        for ( l = 1; l <= nend; l++ )
        {
          n = n - 1;
          en = en - two;
          b[n-1] = ( en * b[n] ) / x + b[n+1];
          em = em - one;
          emp2al = emp2al - one;
          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + b[n-1] * empal ) * emp2al / em;
        }
      }
//
//  Calculate B(1).
//
      b[0] = two * empal * b[1] / x + b[2];

      total = ( total + total ) + b[0];
    }
//
//  Recur backward via difference equation, calculating (but
//  not storing) B(N), until N = NB.
//
    else
    {
      if ( 0 < nend )
      {
        for ( l = 1; l <= nend; l++ )
        {
          n = n - 1;
          en = en - two;
          tempc = tempb;
          tempb = tempa;
          tempa = ( en * tempb ) / x + tempc;
          em = em - one;
          emp2al = emp2al - one;

          if ( n == 1 )
          {
            break;
          }

          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + tempa * empal ) * emp2al / em;
        }
      }
//
//  Store B(NB).
//
      b[n-1] = tempa;

      if ( nb <= 1 )
      {
        total = ( total + total ) + tempa;
      }
//
//  Calculate and Store B(NB-1).
//
      else
      {
        n = n - 1;
        en = en - two;
        b[n-1] = ( en * tempa ) / x + tempb;

        if ( 1 < n  )
        {
          em = em - one;
          emp2al = emp2al - one;

          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + b[n-1] * empal ) * emp2al / em;

          nend = n - 2;
//
//  Calculate via difference equation and store B(N), until N = 2.
//
          if ( 0 < nend )
          {
            for ( l = 1; l <= nend; l++ )
            {
              n = n - 1;
              en = en - two;
              b[n-1] = ( en * b[n] ) / x + b[n+1];
              em = em - one;
              emp2al = emp2al - one;
              if ( n == 2 )
              {
                emp2al = one;
              }
              empal = empal - one;
              total = ( total + b[n-1] * empal ) * emp2al / em;
            }
          }
//
//  Calculate B(1).
//
          b[0] = two * empal * b[1] / x + b[2];
        }
        total = ( total + total ) + b[0];
      }
    }
//
//  Normalize.  Divide all B(N) by TOTAL.
//
    if ( alpha != zero )
    {
       total = total * r8_gamma ( one + alpha ) * pow ( x * half, -alpha );
    }

    if ( ize == 1 )
    {
      total = total * exp ( -x );
    }

    tempa = enmten;

    if ( 1.0 < total )
    {
      tempa = tempa * total;
    }

    for ( n = 1; n <= nb; n++ )
    {
      if ( b[n-1] < tempa )
      {
        b[n-1] = zero;
      }
      b[n-1] = b[n-1] / total;
    }

    return ncalc;
  }
//
//  Two-term ascending series for small X.
//
  else
  {
    tempa = one;
    empal = one + alpha;
    halfx = zero;

    if ( enmten < x )
    {
      halfx = half * x;
    }

    if ( alpha != zero )
    {
      tempa = pow ( halfx, alpha ) / r8_gamma ( empal );
    }

    if ( ize == 2 )
    {
      tempa = tempa * exp ( - x );
    }

    tempb = zero;

    if ( one < x + one )
    {
      tempb = halfx * halfx;
    }

    b[0] = tempa + tempa * tempb / empal;

    if ( x != zero && b[0] == zero )
    {
      ncalc = 0;
    }

    if ( 1 < nb )
    {
      if ( x == zero )
      {
        for ( i = 1; i < nb; i++ )
        {
          b[i] = zero;
        }
      }
//
//  Calculate higher-order functions.
//
      else
      {
        tempc = halfx;
        tover = ( enmten + enmten ) / x;

        if ( tempb != zero )
        {
          tover = enmten / tempb;
        }

        for ( n = 2; n <= nb; n++ )
        {
          tempa = tempa / empal;
          empal = empal + one;
          tempa = tempa * tempc;

          if ( tempa <= tover * empal )
          {
            tempa = zero;
          }

          b[n-1] = tempa + tempa * tempb / empal;

          if ( b[n-1] == zero && n < ncalc )
          {
            ncalc = n - 1;
          }
        }
      }
    }
  }

  return ncalc;
}


#endif
