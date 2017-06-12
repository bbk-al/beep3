/*
 * fmm_globals.h
 *
 *  Created on: 16 Aug 2010
 *      Author: david
 */

#ifndef FMM_GLOBALS_H_
#define FMM_GLOBALS_H_
#include "math.h"
#include "fmm_math_funcs.h"
#include "../common/octree_indexer.h"

namespace fmm
{

// returns 2 to the power of a number
// two_powers is an array of ulongs defined
// in octree_indexer.h
template <typename T, typename IntType>
inline T two_pow(IntType exponent)
{
    assert(abs(exponent) < 64);
    if (exponent < 0)
    {
        return static_cast<T> (1.0 / two_powers[-exponent]);
    }
    else
    {
        return static_cast<T>(two_powers[exponent]);
    }
}

template <typename T>
inline T two_pow(unsigned short exponent)
{
    assert(abs(exponent) < 64);
    return static_cast<T>(two_powers[exponent]);
}

template<int NTERMS>
class FMM_Globals
{

public:

    DblMatrix3D rdpi2;   // init'd by yhrotgen
    DblMatrix3D rdmpi2;  // init'd by yhrotgen
    DblMatrix3D rdsq3;   // init'd by yhrotgen
    DblMatrix3D rdmsq3;  // init'd by yhrotgen
    DblMatrix2D scale_factors; // yhfrmini

    // Constructor
    FMM_Globals() :
        rdpi2(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(-NTERMS,NTERMS+1)),
        rdmpi2(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(-NTERMS,NTERMS+1)),
        rdsq3(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(-NTERMS,NTERMS+1)),
        rdmsq3(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(-NTERMS,NTERMS+1)),
        scale_factors(Range(0,NTERMS+1),Range(0,NTERMS+1))
    {
        // DEBUG
        //std::cout << "Creating FMM globals." << std::endl;

        yhrotgen(rdpi2, rdmpi2, rdsq3, rdmsq3);
        yhfrmini(scale_factors);
    }

    static double get_scale(double beta, double top_level_edge_length, unsigned short level)
    {
    	double scale_factor = beta*top_level_edge_length;
    	if (scale_factor >= 1.0) { scale_factor = 1.0;}// / scale_factor; }

    	unsigned int two_power = static_cast<unsigned int>(pow(2,level));
    	return scale_factor / two_power;
    }

    static void yhrotgen(DblMatrix3D& rdpi2,
    		 	 	 	 DblMatrix3D& rdmpi2,
    		 	 	 	 DblMatrix3D& rdsq3,
    		 	 	 	 DblMatrix3D& rdmsq3)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //    precomputes the rotation matrix for
        //     1. mp->mp
        //     2. local->local
        //     3. east-west expansion
        //     4. north-south expansion
        //
        //  on input :
        //    NTERMS : the number of terms in the multipole expansion.
        //
        //  on output :
        //    rdpi2, rdmpi2 : the rotation matrix for 3 and 4.
        //    rdsq3, rdmsq3 : the rotation matrix for 1 and 2.
        //
        //  workspace :
        //    carray : the square root of the binomial numbers.
        //      these numbers are only used here in this subroutine.
        //
        //  subroutine called :
        //    bnlcft(), fstrtn(), datan, dacos,
        //
        //  called from : main()
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        DblMatrix2D carray(Range(0,4*NTERMS+1),Range(0,4*NTERMS+1));

        // get binomial coefficients
        bnlcft(carray,4*NTERMS);

        // fast rotation calls
        yhfstrtn(rdpi2,carray,pi/2.0);
        yhfstrtn(rdmpi2,carray,-pi/2.0);

        double theta = acos(sqrt(3.0)/3.0);
        yhfstrtn(rdsq3,carray,theta);
        theta = acos(-sqrt(3.0)/3.0);
        yhfstrtn(rdmsq3,carray,theta);

        return;
    }

    static void yhfrmini(DblMatrix2D& factorial_scale_factors)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //----- initialization entry point - create factorial scaling factors
        //
        //     ...... not the most stable way of doing this ......
        //
        //     the vector fact() is only used here. and c(l,m) will used
        //       in subroutine form_mp() and entry brtaev()
        //
        //     note :
        //       the current program is changed a little bit from entry to
        //       a subroutine, since the result c is needed at several places.
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        // calculate some factorials -- these will get quite large quite quickly!!
        DblMatrix1D facts(Range(0,2*NTERMS+1));
        facts(0) = 1.0;
        for (int ell=1; ell<=2*NTERMS; ++ell)
        {
            facts(ell) = facts(ell-1)*double(ell);
        }

        for (int ell=0; ell <= NTERMS; ++ell)
        {
            for (int m=0; m <= ell; ++m)
            {
            	factorial_scale_factors(ell,m) = facts(ell-m)/facts(ell+m)*double(2*ell+1);
            }
        }

        return;

    }

    static void bnlcft(DblMatrix2D& c, int nmax)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //    computes the binomial coefficients c_NTERMS^n, where n=0,1,2,...,NTERMS.
        //
        //  on input:
        //
        //    NTERMS: an integer indicates the number we are going to choose from.
        //
        //  on output:
        //
        //    c:    an array consists of the squre root of the
        //                 binomial coefficients.
        //
        //  note : this is a different version from the laplace equation.
        //    since the binomial coefficients are not needed, but
        //    only the square root is needed, so we will not store
        //    the binomial coefficients. and we will use that space
        //    for some must-be computed coefficients at different levels.
        //
        //  subroutine called : dsqrt()
        //
        //  called from : rotgen()
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        // init binomial coefficients
        for (int n=0; n <= nmax; ++n)
        {
            c(n,0) = 1.0;
        }
        for (int m=1; m <= nmax; ++m)
        {
            c(m,m) = 1.0;
            for (int n=m+1; n <= nmax; ++n)
            {
                c(n,m) = c(n-1,m) + c(n-1,m-1);
            }
        }

        // compute square roots
        for (int m=1; m <= nmax; ++m)
        {
            for (int n=m+1; n <= nmax; ++n)
            {
                c(n,m) = sqrt(c(n,m));
            }
        }

        return;
    }

    static void yhfstrtn(DblMatrix3D& d, // where the rotation matrix ends up
    					 const DblMatrix2D& sqc, // square roots of binomial coeff's
    					 double theta // rotation angle
                 	 	)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //    implement the fast version of rotation matrices from
        //      the recurrences formulas.
        //
        //  on input:
        //    NTERMS: an integer indicates the dimension of d.
        //    sqc: an array contains the square root of the
        //       binormial coefficients.
        //    theta:  the rotate angle about the y-axis.
        //
        //  on output:
        //    d: an array which contains the rotation matrix.
        //
        //  note: only half of d are evaluated, the other
        //    half can be obtained by using the symmetricity.
        //
        //  called from : rotgen()
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        const double precis = 1.0e-19;
        const double ww = 0.7071067811865476;

        double ctheta=cos(theta);
        if (fabs(ctheta) <= precis) { ctheta=0.0; }
        double stheta=sin(-theta);
        if (fabs(stheta) <= precis) { stheta=0.0; }
        double hsthta=ww*stheta;
        double cthtap=ww*(1.0+ctheta);
        double cthtan=-ww*(1.0-ctheta);

    //
    //-----initial setup for some coefficient matrix.
    //
        d(0,0,0) = 1.0;
    //
        for (int ij=1; ij <= NTERMS; ++ij)
        {
    //
    //-------compute the result for m'=0 case, use formula (1).
    //
            for (int im=-ij; im <= -1; ++im)
            {
                d(ij,0,im) = -sqc(ij-im,2)*d(ij-1,0,im+1);
                if (im > (1-ij)) {
                    d(ij,0,im)=d(ij,0,im)+sqc(ij+im,2)*d(ij-1,0,im-1);
                }
                d(ij,0,im)=d(ij,0,im)*hsthta;
                if (im > -ij) {
                    d(ij,0,im)=d(ij,0,im) + d(ij-1,0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1);
                }
                d(ij,0,im) = d(ij,0,im) / ij;
            }

            d(ij,0,0)=d(ij-1,0,0)*ctheta;
            if (ij > 1) {
                d(ij,0,0)=d(ij,0,0)+hsthta*sqc(ij,2)*(d(ij-1,0,-1)+d(ij-1,0,1))/ij;
            }
    //
            for (int im=1; im <= ij; ++im)
            {
                d(ij,0,im)=-sqc(ij+im,2)*d(ij-1,0,im-1);
                if (im < (ij-1)) {
                    d(ij,0,im)=d(ij,0,im)+sqc(ij-im,2)*d(ij-1,0,im+1);
                }
                d(ij,0,im)=d(ij,0,im)*hsthta;
                if (im < ij) {
                    d(ij,0,im)=d(ij,0,im)+d(ij-1,0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1);
                }
                d(ij,0,im)=d(ij,0,im)/ij;
            }
    //
    //-------compute the result for 0<m'<=j case][ use formula (2).
    //
            for (int imp=1; imp <= ij; ++imp)
            {
                for (int im=-ij;im <= -1; ++im)
                {
                    d(ij,imp,im)=d(ij-1,imp-1,im+1)*cthtan*sqc(ij-im,2);
                    if (im > (1-ij)) {
                        d(ij,imp,im)=d(ij,imp,im)-d(ij-1,imp-1,im-1)*cthtap*sqc(ij+im,2);
                    }
                    if (im > -ij) {
                        d(ij,imp,im)=d(ij,imp,im)+  d(ij-1,imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1);
                    }
                    d(ij,imp,im)=d(ij,imp,im)*ww/sqc(ij+imp,2);
                }

                d(ij,imp,0)=ij*stheta*d(ij-1,imp-1,0);
                if (ij > 1) {
                    d(ij,imp,0)=d(ij,imp,0)-sqc(ij,2)*(d(ij-1,imp-1,-1)*cthtap+d(ij-1,imp-1,1)*cthtan);
                }
                d(ij,imp,0)=d(ij,imp,0)*ww/sqc(ij+imp,2);

                for (int im=1; im <= ij; ++im)
                {
                    d(ij,imp,im)=d(ij-1,imp-1,im-1)*cthtap*sqc(ij+im,2);
                    if (im < (ij-1)) {
                        d(ij,imp,im) = d(ij,imp,im) - d(ij-1,imp-1,im+1)*cthtan*sqc(ij-im,2);
                    }
                    if (im < ij) {
                        d(ij,imp,im)=d(ij,imp,im) + d(ij-1,imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1);
                    }
                    d(ij,imp,im)=d(ij,imp,im)*ww/sqc(ij+imp,2);
                }
    //
    //---------note: the lower part of the matrix can be computed using
    //           symmetry, i.e. formula (3.80) in biedenharn & louck's
    //           book.
    //
            }
          }
    //
    //-----now scale the rotation matrix to avoid y_n^m
    //       note : since in yukawa use p_n^m instead of
    //            y_n^m.
    //
        // create factorials.
        //double facts(2*NTERMS+1);
        DblMatrix1D facts(Range(0,2*NTERMS+1)); // this is better than raw array since we get bounds checking in debug mode
        facts(0)=1.0;
        for (int n=1; n <= 2*NTERMS; ++n) {
            facts(n) = facts(n-1)*double(n);
        }

        for (int n=0; n <= NTERMS; ++n)
        {
            for (int m=0; m <= n; ++m)
            {
                for (int mp=-n; mp <= n; ++mp)
                {
                    unsigned int mpabs = abs(mp);
                    d(n,m,mp)=d(n,m,mp)*sqrt(facts(n+m)/facts(n+mpabs)*facts(n-mpabs)/facts(n-m) );
                }
            }
        }

        return;
    }

};

template<int NTERMS, int NLAMBS>
class Level_Dependent_FMM_Globals
{

public:

    DblMatrix1D quad_nodes;    // Yarvin/Rokhlin quadrature points
    DblMatrix1D quad_weights;  // Yarvin/Rokhlin quadrature weights
    UShrtMatrix1D numfour;
    UShrtMatrix1D numphys;

    CmplxMatrix2D xs;
    CmplxMatrix2D ys;
    DblMatrix2D zs;
    DblMatrix3D rlsc; // init by yrlscini
    DblMatrix3D multipole_shift_coefficients;
    DblMatrix3D local_shift_coefficients;
    double betascal;

    unsigned short nexptot;
    unsigned short nexptotp;

    unsigned int fsize;
    CmplxMatrix1D fexpe;
    CmplxMatrix1D fexpo;
    CmplxMatrix1D fexpback;

    Level_Dependent_FMM_Globals(double beta, unsigned short level, double edge_length) :
        quad_nodes(Range(0,NLAMBS)),
        quad_weights(Range(0,NLAMBS)),
        numfour(Range(0,NLAMBS)),
        numphys(Range(0,NLAMBS)),
        rlsc(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(0,NLAMBS)),
        multipole_shift_coefficients(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(0,NTERMS+1)),
    	local_shift_coefficients(Range(0,NTERMS+1),Range(0,NTERMS+1),Range(0,NTERMS+1))
    {

        // DEBUG
        //std::cout << "Creating level-dependent FMM globals: beta=" << beta << " level=" << level << std::endl;

    	double scale = FMM_Globals<NTERMS>::get_scale(beta, edge_length, level);
    	double level_edge_length = edge_length / two_pow<double>(level);
    	betascal = beta*level_edge_length/2.0;

        // Multipole translation coefficients (for upward pass)
    	{
			double r0 = sqrt(3.0) * level_edge_length / 2.0;
			ympshiftcoef(scale,
			 		     beta,
						 r0,
						 multipole_shift_coefficients);

    	}

        // Yarvin/Rokhlin quadrature points
        init_quadrature_points(quad_nodes, quad_weights);
        numthetahalf(numfour);
        numthetafour(numphys);

        // number of 'physical modes' required for the inner integral.
        // i.e. the plane-wave integration formula involving alpha
        // (see Greengard & Huang 2002 page 649; eqn 26) which is
        // a trapezium rule integration from zero to 2pi, carried out at each
        // step of the Yarvin/Rokhlin quadrature (which is the outer integral),
        // Basically we need to know how many points are needed to get
        // this integral correct, which depends on the value of beta, hence
        // this slightly non-clear test below where the number of physical
        // modes is optimized.
        for (int i=0; i < NLAMBS; ++i)
        {
            double lambda_i = quad_nodes(i);
            double uu = sqrt(lambda_i*lambda_i+2.0*lambda_i*scale/2.0);
            int indd=i;
            int mmax=numphys(i);
            for (int jj=i; jj < NLAMBS; ++jj)
            {
                if (uu <= quad_nodes(jj)) {
                    indd=jj;
                } else if (numphys(jj) > mmax) {
                    mmax=numphys(jj);
                }
            }
            numphys(i) = mmax > numphys(indd) ? mmax : numphys(indd);
        }

        // add up total numbers of fourier and physical modes
        nexptot = 0;
        nexptotp = 0;
        for(unsigned short i=0; i < NLAMBS; ++i)
        {
            nexptot  += numfour(i);
            nexptotp += numphys(i);
        }
        nexptotp = nexptotp/2;  // because of symmetry

        fsize = calc_fsize(numfour, numphys);
        fexpe.realloc(ltl::Shape<1>(Range(0,fsize/2)));
        fexpo.realloc(ltl::Shape<1>(Range(0,fsize/2)));
        fexpback.realloc(ltl::Shape<1>(Range(0,fsize)));

        xs.realloc(ltl::Shape<2>(Range(0,3),Range(0,nexptotp)));
        ys.realloc(ltl::Shape<2>(Range(0,3),Range(0,nexptotp)));
        zs.realloc(ltl::Shape<2>(Range(0,3),Range(0,nexptotp)));

        // init plane wave conversion matrices
        ymkfexp(numfour, numphys, fexpe, fexpo, fexpback);

        yrlscini(scale/2.0, betascal, rlsc, quad_nodes);
        ymkexps(betascal, quad_nodes, numphys, xs, ys, zs);

        // shift coefficients for this level
        {
        	// nb: this r0 (shift distance) is different from that for
        	// multipole shifting since the shift is being carried out on a
        	// different level in each case
			double r0 = sqrt(3.0) * level_edge_length / 4.0;
			ylcshftcoef(scale,
						beta,
						r0,
						local_shift_coefficients);
        }
    }

    static void init_quadrature_points(DblMatrix1D &x, DblMatrix1D &w)
    {
    	// straight outta Yarvin- ftp.cs.yale.edu/pub/yarvin/vwts.f

    	assert(NTERMS==2 || NTERMS==4 || NTERMS == 9 || NTERMS == 18);
        if (NTERMS==2)
        {
			x(1) = 0.31283512406199648347993047536874656e+00;
			x(2) = 0.10894076617764414383060511681833304e+01;
			w(1) = 0.76974593066577023936503110235207714e+00;
			w(2) = 0.84374978799072331003827684980933554e+00;
        }
        else if(NTERMS == 4)
        {
			x(1)=0.21626008006106806069723802465887275e+00;
			x(2)=0.85186788888114162165976495089125820e+00;
			x(3)=0.16705322595388429895990611839806661e+01;
			x(4)=0.25585029089804702806532077374868095e+01;
			w(1)=0.50471599046036208502385989049798809e+00;
			w(2)=0.76071728649476022532383012730861083e+00;
			w(3)=0.87370748029129807754600278713041916e+00;
			w(4)=0.93475133313404745738495194018469192e+00;
    	}
        else if (NTERMS==9)
        {
            x(0) = .099273996739714473469540223504736787;
            x(1) = .47725674637049431137114652301534079;
            x(2) = 1.0553366138218296388373573790886439;
            x(3) = 1.7675934335400844688024335482623428;
            x(4) = 2.5734262935147067530294862081063911;
            x(5) = 3.4482433920158257478760788217186928;
            x(6) = 4.3768098355472631055818055756390095;
            x(7) = 5.3489575720546005399569367000367492;
            x(8) = 6.3576578531337464283978988532908261;
            w(0) = .24776441819008371281185532097879332;
            w(1) = .49188566500464336872511239562300034;
            w(2) = .65378749137677805158830324216978624;
            w(3) = .76433038408784093054038066838984378;
            w(4) = .84376180565628111640563702167128213;
            w(5) = .90445883985098263213586733400006779;
            w(6) = .9537861313683345665381807521043811;
            w(7) = .99670261613218547047665651916759089;
            w(8) = 1.0429422730252668749528766056755558;
        }
        else if (NTERMS==18)
        {
            x(0) = .052788527661177607475107009804560221;
            x(1) = .26949859838931256028615734976483509;
            x(2) = .6322035317468939208396250251098536;
            x(3) = 1.1130756427760852833586113774799742;
            x(4) = 1.6893949614021379623807206371566281;
            x(5) = 2.3437620046953044905535534780938178;
            x(6) = 3.0626998290780611533534738555317745;
            x(7) = 3.8356294126529686394633245072327554;
            x(8) = 4.6542473432156272750148673367220908;
            x(9) = 5.5120938659358147404532246582675725;
            x(10) = 6.4042126837727888499784967279992998;
            x(11) = 7.3268800190617540124549122992902994;
            x(12) = 8.2774009925823861522076185792684555;
            x(13) = 9.2539718060248947750778825138695538;
            x(14) = 10.255602723746401139237605093512684;
            x(15) = 11.282088297877740146191172243561596;
            x(16) = 12.334067909676926788620221486780792;
            x(17) = 13.414920240172401477707353478763252;
            w(0) = .13438265914335215112096477696468355;
            w(1) = .29457752727395436487256574764614925;
            w(2) = .42607819361148618897416895379137713;
            w(3) = .53189220776549905878027857397682965;
            w(4) = .61787306245538586857435348065337166;
            w(5) = .68863156078905074508611505734734237;
            w(6) = .74749099381426187260757387775811367;
            w(7) = .79699192718599998208617307682288811;
            w(8) = .83917454386997591964103548889397644;
            w(9) = .8757009228374531550898041132313665;
            w(10) = .90792943590067498593754180546966381;
            w(11) = .93698393742461816291466902839601971;
            w(12) = .96382546688788062194674921556725167;
            w(13) = .98932985769673820186653756536543369;
            w(14) = 1.0143828459791703888726033255807124;
            w(15) = 1.0400365437416452252250564924906939;
            w(16) = 1.0681548926956736522697610780596733;
            w(17) = 1.1090758097553685690428437737864442;
        }
    }

    static void numthetahalf(UShrtMatrix1D& numtets)
    {

        // this routine returns the number of fourier modes needed in the
        // phi integral for each of the discrete lambda values given
        // by Yarvin/Rokhlin quadrature:
        //
        // input arguments:
        //       NLAMBS - number of nodes in the lambda quadrature.  this must
        //         be either 9 or 18.
        //
        // output arguments:
        //       numtets(i) - number of fourier modes needed for phi
        //                    integral with lambda_i.

        assert(NLAMBS == 2 || NLAMBS == 4 || NLAMBS == 9 || NLAMBS == 18);
        if (NLAMBS == 2)
        {
            numtets(0) = 1;
            numtets(1) = 2;
        }
        else if(NLAMBS == 4)
        {
            numtets(0) = 2;
            numtets(1) = 4;
            numtets(2) = 6;
            numtets(3) = 4;
    	}
        else if (NLAMBS == 9) {
            numtets(0) = 2;
            numtets(1) = 4;
            numtets(2) = 4;
            numtets(3) = 6;
            numtets(4) = 6;
            numtets(5) = 4;
            numtets(6) = 6;
            numtets(7) = 4;
            numtets(8) = 2;
        }
		else if (NLAMBS == 18)
		{
            numtets(0) = 4;
            numtets(1) = 6;
            numtets(2) = 6;
            numtets(3) = 8;
            numtets(4) = 8;
            numtets(5) = 8;
            numtets(6) = 10;
            numtets(7) = 10;
            numtets(8) = 10;
            numtets(9) = 10;
            numtets(10) = 12;
            numtets(11) = 12;
            numtets(12) = 12;
            numtets(13) = 12;
            numtets(14) = 12;
            numtets(15) = 12;
            numtets(16) = 8;
            numtets(17) = 2;

        } else {
            std::cerr << "ERROR! Bad NTERMS in numthetahalf" << std::endl;
        }

        return;
    }

    static void numthetafour(UShrtMatrix1D& numtets)
    {

        // this routine returns the number of fourier modes needed in the
        // phi integral for each of the discrete lambda values given
        // by Yarvin/Rokhlin quadrature:
        //
        // input arguments:
        //       nlambs - number of nodes in the lambda quadrature.  this must
        //         be an integer in the range (2,39).
        //
        // output arguments:
        //       numtets(i) - number of fourier modes needed for phi
        //                    integral with lambda_i.

        assert(NLAMBS == 2 || NLAMBS == 4 || NLAMBS == 9 || NLAMBS == 18);
        if (NLAMBS == 2) {
            numtets(0) = 2;
            numtets(1) = 4;
        }
        else if(NLAMBS == 4)
		{
			numtets(0) = 4;
			numtets(1) = 8;
			numtets(2) = 12;
			numtets(3) = 8;
		}
        else if (NLAMBS == 9) {
            numtets(0) = 4;
            numtets(1) = 8;
            numtets(2) = 12;
            numtets(3) = 16;
            numtets(4) = 20;
            numtets(5) = 20;
            numtets(6) = 24;
            numtets(7) = 8;
            numtets(8) = 2;
        } else if (NLAMBS == 18) {
            numtets(0) = 6;
            numtets(1) = 8;
            numtets(2) = 12;
            numtets(3) = 16;
            numtets(4) = 20;
            numtets(5) = 26;
            numtets(6) = 30;
            numtets(7) = 34;
            numtets(8) = 38;
            numtets(9) = 44;
            numtets(10) = 48;
            numtets(11) = 52;
            numtets(12) = 56;
            numtets(13) = 60;
            numtets(14) = 60;
            numtets(15) = 52;
            numtets(16) = 4;
            numtets(17) = 2;
        } else {
            std::cerr << "ERROR! Bad NTERMS in numthetafour" << std::endl;
        }

        return;
    }


    void ymkexps(const double beta,
                 const DblMatrix1D& quad_nodes,
                 const UShrtMatrix1D& numphys,
                 CmplxMatrix2D& xs,
                 CmplxMatrix2D& ys,
                 DblMatrix2D& zs)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //     this subroutine computes the tables of exponentials needed
        //     for translating exponential representations of harmonic
        //     functions, discretized via norman's quadratures.
        //
        //     u   = \int_0^\infty e^{-(lambda+beta) z}
        //     \int_0^{2\pi} e^{i \dsqrt(lambda*lambda+2*beta*lambda)
        //                         (x cos(u)+y sin(u))}
        //           mexpphys(lambda,u) du dlambda
        //
        //     mexpphys(*):  discrete values of the moment function
        //                   m(\lambda,u), ordered as follows.
        //
        //         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
        //              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
        //         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
        //              m(\lambda_2,0),...,
        //                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
        //         etc.
        //         note : in the current version, only half of the modes are
        //                stored because of the symmetry of u and u+pi.
        //
        //  on input:
        //
        //     beta : the scaled frequency.
        //     quad_nodes(NLAMBS)   discretization points in lambda integral
        //     NLAMBS          number of discret. pts. in lambda integral
        //     numphys(j)     number of nodes in u integral needed
        //                    for corresponding lambda =  lambda_j.
        //     nexptotp        sum_j numphys(j)
        //
        //  on output:
        //
        //        define w1=\lambda_j+beta, and w2= sqrt(lambda_j**2+2*beta*lambda_j)
        //     xs(1,nexptotp)   e^{i w2 (cos(u_k)}  in above ordering
        //     xs(2,nexptotp)   e^{i w2 (2 cos(u_k)}  in above ordering.
        //     xs(3,nexptotp)   e^{i w2 (3 cos(u_k)}  in above ordering.
        //     ys(1,nexptotp)   e^{i w2 (sin(u_k)}  in above ordering.
        //     ys(2,nexptotp)   e^{i w2 (2 sin(u_k)}  in above ordering.
        //     ys(3,nexptotp)   e^{i w2 (3 sin(u_k)}  in above ordering.
        //     zs(1,nexptotp)   e^{-w1}     in above ordering.
        //     zs(2,nexptotp)   e^{-2 w1}   in above ordering.
        //     zs(3,nexptotp)   e^{-3 w1}   in above ordering.
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //

    //
    //-----loop over each lambda value
    //

        int ntot = 0;
        for (int nl=0; nl < NLAMBS; ++nl)
        {
            double w1= quad_nodes(nl)+beta;
            double w2=sqrt( quad_nodes(nl)*(quad_nodes(nl)+beta*2.0) );
            double hu= 2.0*pi / double(numphys(nl));
            for (int mth=0; mth < numphys(nl)/2; ++mth)
            {
                double u = hu*mth;
                int ncurrent = ntot+mth;

                //std::cout << "fmm::ymkexps : ncurrent: " << ncurrent << std::endl;

                zs(0,ncurrent) = exp(-w1);
                zs(1,ncurrent) = zs(0,ncurrent)*zs(0,ncurrent);
                zs(2,ncurrent) = zs(1,ncurrent)*zs(0,ncurrent);
                xs(0,ncurrent) = exp(ima*w2*cos(u));
                xs(1,ncurrent) = xs(0,ncurrent)*xs(0,ncurrent);
                xs(2,ncurrent) = xs(1,ncurrent)*xs(0,ncurrent);
                ys(0,ncurrent) = exp(ima*w2*sin(u));
                ys(1,ncurrent) = ys(0,ncurrent)*ys(0,ncurrent);
                ys(2,ncurrent) = ys(1,ncurrent)*ys(0,ncurrent);
            }
            ntot += numphys(nl)/2;
        }

        return;
    }

    void ymkfexp(const UShrtMatrix1D& numfour,
                 const UShrtMatrix1D& numphys,
                 CmplxMatrix1D& fexpe,
                 CmplxMatrix1D& fexpo,
                 CmplxMatrix1D& fexpback)
    {
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //     precomputes the e^(im*alpha) needed in mp->pw->local.
        //
        //     this subroutine computes the tables of exponentials needed
        //     for mapping from fourier to physical domain.
        //     in order to minimize storage, they are organized in a
        //     one-dimenional array corresponding to the order in which they
        //     are accessed by subroutine ftophys.
        //
        //     size of fexpe, fexpo =          40000   for NLAMBS = 39
        //     size of fexpe, fexpo =          15000   for NLAMBS = 30
        //     size of fexpe, fexpo =           4000   for NLAMBS = 20
        //     size of fexpe, fexpo =            400   for NLAMBS = 10
        //
        //
        //  on input :
        //
        //       NLAMBS : the total number of nodes for the outer integral.
        //       numfour : contains the number of nodes for the fourier
        //                 representation.
        //       numphys : contains the number of nodes for the inner integral
        //                 for the plane wave expansion.
        //
        //  on output :
        //
        //       fexpe : the exponentials for the fourier modes. e^(im*alpha)
        //               where m is all the fourier modes and alpha comes
        //               from the physical modes. odd terms.
        //       fexpo : the even terms.
        //       fexpback : the exponentials used for the translation
        //                  from plane wave to local.
        //
        //     functions called :
        //
        //     called from :
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        unsigned int nexte=0, nexto=0;
        for (int i=0; i < NLAMBS; ++i)
        {
            int nalpha = numphys(i);
            int nalpha2 = nalpha/2;
            double halpha=2.0*pi / nalpha;

            for(int j=0; j < nalpha2; ++j)
            {
                double alpha=halpha*j;
                for (int mm=2; mm <= numfour(i); mm+=2)
                {
                    fexpe(nexte++) = exp(ima*double(mm-1)*alpha);
                }

                for(int mm=3; mm <= numfour(i); mm+=2)
                {
                    fexpo(nexto++) = exp(ima*double(mm-1)*alpha);
                }
            }
        }

        unsigned int next = 0;
        for (int i=0; i < NLAMBS; ++i)
        {
            int nalpha = numphys(i);
            int nalpha2 = nalpha/2;
            double halpha = 2.0*pi/nalpha;
            for (int mm=3; mm <= numfour(i); mm+=2)
            {
                for(int j=0; j < nalpha2; ++j)
                {
                    double alpha=halpha*j;
                    fexpback(next++) = exp(- ima * double(mm-1) * alpha);
                }
            }

            for (int mm=2; mm <= numfour(i); mm+=2)
            {
                for(int j=0; j < nalpha2; ++j)
                {
                    double alpha=halpha*j;
                    fexpback(next++) = exp(- ima * double(mm-1) * alpha);
                }
            }
        }

        return;
    }

    inline unsigned int calc_fsize(const UShrtMatrix1D& numfour,
                                       const UShrtMatrix1D& numphys)
    {
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //     calculates the size of the fexpe/fexpo/fexpback arrays which
        //     are precomputed in ymkfexp.
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        unsigned int fsize=0;
        for (int i=0; i < NLAMBS; ++i)
        {
            fsize += numfour(i)*numphys(i);
        }
        return fsize/2;
    }

    inline void ympshiftcoef(double scal,
                      	  	    double beta,
                      	  	    double r0,
                      	  	    DblMatrix3D& c0)
    {
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //    this subroutine precomputes the coefficients for the
        //      shifting of the multipole expansion.
        //
        //  on input :
        //
        //    scal : the scale factor of the child level. the parent
        //           has the scaling factor of 2*scal
        //    beta : the coefficient for the helmholtz equation.
        //    r0 : the real shifting distance.
        //    NTERMS : number of terms in the multipole expansion.
        //
        //  on output :
        //
        //    c0() : the coefficients of the shifting of the multipole expansion.
        //    info : error messages.
        //
        //  workspace :
        //    fact(0:200): assigned locally.
        //    bj(0:200): assigned locally.
        //
        //  function called : in().
        //
        //  note : the new scaled multipole expansion coefficients
        //         m_l^m= sum_{n=m}^{NTERMS} m_n^m *c(m,l,n)
        //         where m_n^m is also the scaled multipole expansion
        //         of the child (by scal^n).
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

    	// original version
    	// create factorials.
        DblMatrix1D facts(Range(0,2*NTERMS+1));
        facts(0)=1.0;
        for (int n=1; n <= 2*NTERMS; ++n) {
            facts(n) = facts(n-1)*double(n);
        }

        double r0k = r0 * beta;
        DblMatrix1D bj(Range(0,2*NTERMS+1));
        i_n(r0k, r0k, 2*NTERMS, bj);

        for (int m=0; m <= NTERMS; ++m)
        {
            for (int n_dash=m; n_dash <= NTERMS; ++n_dash)
            {
                for (int n=m; n <= NTERMS; ++n)
                {
                    c0(m,n_dash,n)=0.0;
                    int upper = (n <= n_dash) ? n : n_dash;
                    for (int k=m; k <= upper; ++k)
                    {
                        double premult = pow(scal,n-n_dash)*two_pow<double>(-n_dash-k)* ( ( (n_dash+n) % 2 == 0) ? 1.0 : -1.0);
                        double factorials = double(2*n_dash+1)*facts(n_dash-m)/facts(k+m)*facts(n+m)*facts(2*k)/facts(k)/facts(k-m)/facts(n_dash-k)/facts(n-k);
                        c0(m,n_dash,n) += premult*factorials*bj(n_dash+n-k)*pow(r0k,n_dash+n-2*k);
                    }
                }
            }
        }
        return;
    }

    inline void ylcshftcoef(const double scal,
							 const double beta,
							 const double r0,
							 DblMatrix3D& c0)
    {
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //
        //    this subroutine precomputes the coefficients for the
        //      shifting of the local expansion.
        //
        //  on input:
        //    scal: the scale factor.
        //    beta: the coefficient for the helmholtz equation.
        //    NTERMS: number of terms in the local expansion.
        //
        //  on output :
        //    c0(): the coefficients of the shifting of the local expansion.
        //    info: error messages.
        //
        //  workspace:
        //    none
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        DblMatrix1D facts(Range(0,2*NTERMS+1)); // this is better than raw array since we get bounds checking in debug mode
        facts(0)=1.0;
        for (int n=1; n <= 2*NTERMS; ++n) {
            facts(n) = facts(n-1)*double(n);
        }
        const double r0k = r0*beta;
        DblMatrix1D bj(Range(0,2*NTERMS+1));
        i_n(r0k, r0k, 2*NTERMS, bj);

        for(int m=0; m <= NTERMS; ++m)
        {
            for(int n_dash=m; n_dash <= NTERMS; ++n_dash)
            {
                for (int n=m; n <= NTERMS; ++n)
                {
                    c0(m,n_dash,n)=0.0;
                    const int upper = (n <= n_dash) ? n : n_dash;
                    for (int k=m; k <= upper; ++k)
                    {
                        double premult = pow(scal,n_dash-n)*pow(2.0,-n_dash-k);
                        double factorials = double(2*n_dash+1)*facts(n_dash-m)/facts(k+m)*facts(n+m)*facts(2*k)/facts(k)/facts(k-m)/facts(n_dash-k)/facts(n-k);
                        c0(m,n_dash,n) += premult*factorials*bj(n_dash+n-k)*pow(r0k,n_dash+n-2*k);
                    }
                }
            }
        }

        return;
    }

    void yrlscini(double scal,
                  double beta,
                  DblMatrix3D& rlsc,
                  const DblMatrix1D& quad_nodes)
    {

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        //
        //  purpose:
        //    precomputes the coefficients for mpoletoexp for different levels.
        //    note: for different level, this subroutine should be
        //          called because of the change of the scaled frequency.
        //    in return, the rlsc will contain the associated legendre polynomial
        //    p_n^m( (beta+lamda_i)/beta )*scal^n
        //
        //  on input :
        //    scal : the scaling factor for the p_n^m. otherwise p_n^m will
        //           overflow.
        //    beta : the scaled frequency/beta for the current level.
        //    NLAMBS : the total number of lamdas for the first integral.
        //    quad_nodes : the nodes/lamdas.
        //    NTERMS : the total number of terms in the multipole expansion.
        //
        //  on output :
        //    rlsc(n,m,nlamda) : the required information for the whole level.
        //
        //  subroutine called :
        //
        //  called from : brfrc()
        //
        //  note : the rlsc() will be used in the subroutine mpoletoexp().
        //
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd

        for (int nl=0; nl < NLAMBS; ++nl)
        {
            LegendreHolder<NTERMS> rlsc_temp;
            const double ul = quad_nodes(nl)/beta + 1.0;
            lgndrgt1(scal, ul, rlsc_temp);

            // copy rlsc_temp into the rlsc matrix passed in as argument
            // bit cumbersome but we don't do this very often...
            for (int i=0; i <= NTERMS; ++i)
            {
                for (int j=0; j <= i; ++j)
                {
                    //std::cout << i << " " << j << std::endl;
                    //std::cout << rlsc_temp(i,j) << std::endl;
                    rlsc(i,j,nl) = rlsc_temp(i,j);
                }
            }
        }
        return;
    }

};

} // end namespace

#endif /* FMM_GLOBALS_H_ */
