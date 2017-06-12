#include "prerequisites.h"

#include "fmm_globals_nodegroup.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "rhs_handler.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "iteration_handler.decl.h"
#include "gmres.decl.h"
#include "main.decl.h"

#include "rhs_handler.h"
#include "fh_values_nodegroup.h"
#include "iteration_handler.h"
#include "gmres.h"

#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>

// BLAS functions
#include <gsl/gsl_cblas.h>

// GMRES parameters
#define GMRES_PRECISION 1e-4
#define GMRES_RESTART 19

extern /* readonly */ CProxy_Main MainProxy;
extern /* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
extern /* readonly */ CProxy_GMRES GMRESProxy;
extern /* readonly */ CProxy_IterationHandler IterationHandlerProxy;
extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;

GMRES::GMRES()
{
    return;
}

void GMRES::solve(std::string output_filename,
                    double kappa,
                    double Dsolvent,
                    unsigned int total_num_patches,
                    FH_Values rhs)
{

// ============================================================================
//
//  GMRES nach Saad, Schultz
//     GMRES: a generalized minimal residual algorithm for solving nonsymmetric
//     linear systems
//     SIAM J Sci Stat Comput 7, 856-869 (1986)
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//                                           Modified by David Fallaize, 2010
//
//  Obtained from the linsolver package, see:
//  http://aam.mathematik.uni-freiburg.de/IAM/Research/projectskr/lin_solver/
// 
//
// ============================================================================
//
    //std::cout << "Starting GMRES." << std::endl;

    const int m=GMRES_RESTART; // retries
    const int n = static_cast<int>(total_num_patches*2); // dimensions of the problem

    double *b = new double[n]; // rhs
    double *x = new double[n]; // guess
    memcpy(b, rhs.get(), sizeof(double)*n);
    //memcpy(x, b, sizeof(double)*n);
    memset(x,0,sizeof(double)*n);

    const double eps = GMRES_PRECISION; // precision
    const bool detailed = true;
    {

        const double fscale = 1.0;
        const double hscale = 1.0;
        for (unsigned int xxx=0; xxx < total_num_patches; ++xxx)
        {
            b[2*xxx] /= fscale;
            b[2*xxx + 1] /= hscale;
        }

        typedef double *doubleP;
        double *V  = new double[n*(m+1)];
        double *U  = new double[m*(m+1)/2];
        double *r  = new double[n];
        double *y  = new double[m+1];
        double *c  = new double[m];
        double *s  = new double[m];
        double **v = new doubleP[m+1];
        for ( int i=0; i<=m; ++i ) v[i]=V+i*n;
        int its=-1;
        {
            double gmres_beta, h, rd, dd, nrm2b;
            int j, io, uij, u0j;
            nrm2b=cblas_dnrm2(n,b,1);

            io=0;
            do  { // "aussere Iteration
            ++io;

            {
                // PREPARE DATA
                FH_Values fh_vals(n, x);
                RHS_Message* return_msg;

                // BLOCKING CALL -- PERFORM ONE FMM/BEM ITERATION!
                IterationHandlerProxy.do_bem_fmm_iteration(CkCallbackResumeThread((void*&)return_msg), fh_vals);
                //std::cout <<"Resumed thread" << std::endl;

                // EXTRACT RESULTS
                memcpy(r, return_msg->data, sizeof(double)*n);
                delete return_msg;
            }

            cblas_daxpy(n,-1.,b,1,r,1);
            gmres_beta=cblas_dnrm2(n,r,1);
            cblas_dcopy(n,r,1,v[0],1);
            cblas_dscal(n,1./gmres_beta,v[0],1);

            y[0]=gmres_beta;
            j=0;
            uij=0;
            do { // innere Iteration j=0,...,m-1

                u0j=uij;

                double *xxx = v[j];
                double *rrr = v[j+1];

                // PREPARE DATA
                FH_Values fh_vals(n, xxx);
                RHS_Message* return_msg;

                // BLOCKING CALL -- PERFORM ONE FMM/BEM ITERATION!
                IterationHandlerProxy.do_bem_fmm_iteration(CkCallbackResumeThread((void*&)return_msg), fh_vals);
                //std::cout <<"Resumed thread" << std::endl;
                // EXTRACT RESULTS
                memcpy(rrr, return_msg->data, sizeof(double)*n);
                delete return_msg;

                //mult(v[j],v[j+1]);
                cblas_dgemv(CblasColMajor,CblasTrans,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);
                cblas_dgemv(CblasColMajor,CblasNoTrans,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);
                h=cblas_dnrm2(n,v[j+1],1);
                cblas_dscal(n,1./h,v[j+1],1);
                for ( int i=0; i<j; ++i ) { // rotiere neue Spalte
                double tmp = c[i]*U[uij]-s[i]*U[uij+1];
                U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
                U[uij]     = tmp;
                ++uij;
                }
                { // berechne neue Rotation
                rd     = U[uij];
                dd     = sqrt(rd*rd+h*h);
                c[j]   = rd/dd;
                s[j]   = -h/dd;
                U[uij] = dd;
                ++uij;
                }
                { // rotiere rechte Seite y (vorher: y[j+1]=0)
                y[j+1] = s[j]*y[j];
                y[j]   = c[j]*y[j];
                }
                ++j;
                if ( detailed )
                cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<y[j]<<"\t" << nrm2b*eps <<endl;
            } while ( j<m && fabs(y[j])>=eps*nrm2b );

            cblas_dtpsv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,j,U,y,1);
            cblas_dgemv(CblasColMajor,CblasNoTrans,n,j,-1.,V,n,y,1,1.,x,1);

            //  } while ( fabs(y[j])>=eps*nrm2b );
            } while (false);

            // R"uckgabe: Zahl der inneren Iterationen
            its = m*(io-1)+j;
        }

        // DEBUG
        std::cout << "GMRES reached convergence in " << its << " iterations." << std::endl;

        // GET FINAL RESULTS
        FH_Values final_results(n, x);
        FH_Values_NodeGroupProxy.set(final_results, CkCallbackResumeThread());
        const double *fvals = FH_Values_NodeGroupProxy.ckLocalBranch()->fvals();
        const double *hvals = FH_Values_NodeGroupProxy.ckLocalBranch()->hvals();
/*        for (int ii=0; ii < 10; ++ii)
        {
            std::cout << fvals[ii] << " " << hvals[ii] << std::endl;
        }*/
        
        //MeshWorkingCopyProxy.calculate_energy(kappa, Dsolvent);
        MeshWorkingCopyProxy.write_output(output_filename);

        delete[] V;
        delete[] U;
        delete[] r;
        delete[] y;
        delete[] c;
        delete[] s;
        delete[] v;
    }

    delete[] x;
    delete[] b;

    return;
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
GMRES::GMRES(CkMigrateMessage *msg) {
    std::cout << "AARGH! GMRES is migrating!" << std::endl;
}

#include "gmres.def.h"
