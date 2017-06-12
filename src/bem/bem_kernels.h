/*
* bem_kernels.h
*
*  Created on: 22 Jul 2010
*      Author: david
*/

#ifndef BEM_KERNELS_H_
#define BEM_KERNELS_H_

#include "../common/math_vector.h"
#include "../common/matrix_types.h"
#include "node_patch.h"
#include "quad_point.h"

class BasicNodePatch;

template<typename PrecisionType>
void singular_BEM_kernels(double kappa, 
                          const BasicNodePatch& patch,
                          PrecisionType& A_out,
                          PrecisionType& B_out,
                          PrecisionType& C_out,
                          PrecisionType& D_out)
{
    const PrecisionType sigm = 1e-9;

    A_out = 0;
    B_out = 0;
    C_out = 0;
    D_out = 0;
    
    PrecisionType inv_epsilon = static_cast<PrecisionType>(1.0 / patch.get_dielectric_ratio());

    boost::shared_ptr<QuadList> extreme_pts = patch.get_galerkin_points();
    PrecisionType normalise_qual_wt=0;

    Matrix_3x3<PrecisionType> gdd;
    Matrix_3x3<PrecisionType> udd;
    
    //#pragma omp parallel for
    for (int src_ctr=0; src_ctr < extreme_pts->size(); ++src_ctr)
    {
        const QuadPoint& src_qp = (*extreme_pts)[src_ctr];
        const VectorT<PrecisionType> src_pt = src_qp.pt();
        const VectorT<PrecisionType> n0 = src_qp.normal();
        normalise_qual_wt += src_qp.weight();
         
        PrecisionType A = 0;
        PrecisionType B = 0;
        PrecisionType C = 0;
        PrecisionType D = 0;
        
        for (int targ_ctr=0; targ_ctr < extreme_pts->size(); ++targ_ctr)
        {
            const QuadPoint& qp = (*extreme_pts)[targ_ctr];

            PrecisionType wt = qp.weight();
            VectorT<PrecisionType> dx = qp.pt() - src_pt;
            VectorT<PrecisionType> normal = qp.normal();
            
            // use inverse multiplication instead of division
            PrecisionType r2 = dx.length2();
            if (fabs(r2) < sigm) { continue; }
            PrecisionType r = sqrt(r2);
            PrecisionType ir2 = 1.0 / r2;
            PrecisionType ir = 1.0 / r;
            PrecisionType ir3 = ir2*ir;
            PrecisionType ir4 = ir2*ir2;
            PrecisionType ir5 = ir4*ir;

            PrecisionType gr0 = ONE_OVER_4PI;
            PrecisionType gr1 = gr0* ir;
            //PrecisionType gr2 = gr0/r2;
            PrecisionType gr3 = gr0*ir3;
            //PrecisionType gr4 = gr0/r4;
            
            PrecisionType ur0=exp(-kappa*r) * ONE_OVER_4PI;
            PrecisionType ur1=ur0*ir;
            PrecisionType ur2=ur0*ir2;
            PrecisionType ur3=ur0*ir3;
            PrecisionType ur3ur2=ur3+kappa*ur2;

            VectorT<PrecisionType> gd = -dx*gr3;
            VectorT<PrecisionType> ud = -dx*(ur3ur2);

            A += wt*(gr1 - ur1);
            B += wt * normal.dot(inv_epsilon*gd - ud);
            C += -wt * n0.dot(gd - ud*inv_epsilon);
            if (kappa != 0)
            {
                PrecisionType gr5 = gr0*ir5;
                PrecisionType ur4=ur0*ir4;
                PrecisionType ur5=ur0*ir5;
                
                PrecisionType pur4ur3=kappa*(3.0*ur4+kappa*ur3)+3.0*ur5;

                for (int jj=0; jj < 3; ++jj)
                {
                    for (int kk=0; kk < 3; ++kk)
                    {
                        PrecisionType kronecker = (jj == kk) ? 1.0 : 0.0;
                        PrecisionType dxjk = dx(jj)*dx(kk);
                        gdd(jj,kk) = 3.0*dxjk*gr5 - kronecker*gr3;
                        udd(jj,kk) = dxjk*pur4ur3 - kronecker*ur3ur2;                        
                    }
                }
                
                for (int jj=0; jj < 3; ++jj)
                {
                    for (int kk=0; kk < 3; ++kk)
                    {
                        D -= wt * n0(kk) * normal(jj) * (gdd(jj,kk) - udd(jj,kk));
                    }
                }
            }
        }
        
        PrecisionType qual_wt = src_qp.weight();
        A_out += A * qual_wt; 
        B_out += B * qual_wt;
        C_out += C * qual_wt;
        D_out += D * qual_wt;
    }
    
    assert(normalise_qual_wt != 0);
    A_out /= normalise_qual_wt;
    B_out /= normalise_qual_wt;
    C_out /= normalise_qual_wt;
    D_out *= inv_epsilon / normalise_qual_wt;
    
}

template<typename PrecisionType>
void evaluate_local_BEM_kernels(PrecisionType kappa,
                                const BasicNodePatch& src_patch,
                                const BasicNodePatch& targ_patch,
                                PrecisionType& A_out,
                                PrecisionType& B_out,
                                PrecisionType& C_out,
                                PrecisionType& D_out)
{
    const PrecisionType sigm = 1e-9;

    A_out = 0;
    B_out = 0;
    C_out = 0;
    D_out = 0;
    
    // skip self-patch interactions -- they're done elsewhere
    if (src_patch.get_idx() == targ_patch.get_idx()) { return; }
    
    PrecisionType inv_epsilon = static_cast<PrecisionType>(1.0 / targ_patch.get_dielectric_ratio());

    boost::shared_ptr<QuadList> qual_points = src_patch.get_qualocation_points();
    boost::shared_ptr<QuadList> quad_points = targ_patch.get_quadrature_points();

    Matrix_3x3<PrecisionType> gdd;
    Matrix_3x3<PrecisionType> udd;
    
    for (std::vector<QuadPoint>::const_iterator src_it=qual_points->begin(), src_end=qual_points->end();
         src_it != src_end;
         ++src_it)
    {
        const VectorT<PrecisionType> src_pt = src_it->pt();
        const VectorT<PrecisionType> n0 = src_it->normal();
        
        PrecisionType A = 0;
        PrecisionType B = 0;
        PrecisionType C = 0;
        PrecisionType D = 0;

        for (std::vector<QuadPoint>::const_iterator it=quad_points->begin(), end=quad_points->end();
            it != end;
            ++it)
        {
            const QuadPoint& qp = *it;

            PrecisionType wt = qp.weight();
            VectorT<PrecisionType> dx = qp.pt() - src_pt;
            VectorT<PrecisionType> normal = qp.normal();

            // use inverse multiplication instead of division
            PrecisionType r2 = dx.length2();
            if (fabs(r2) < sigm) { continue; } 
            PrecisionType r = sqrt(r2);
            PrecisionType ir2 = 1.0 / r2;
            PrecisionType ir = 1.0 / r;
            PrecisionType ir3 = ir2*ir;
            PrecisionType ir4 = ir2*ir2;
            PrecisionType ir5 = ir4*ir;

            PrecisionType gr0 = ONE_OVER_4PI;
            PrecisionType gr1 = gr0* ir;
            //PrecisionType gr2 = gr0/r2;
            PrecisionType gr3 = gr0*ir3;
            //PrecisionType gr4 = gr0/r4;
            PrecisionType ur0=exp(-kappa*r) * ONE_OVER_4PI;
            PrecisionType ur1=ur0*ir;
            PrecisionType ur2=ur0*ir2;
            PrecisionType ur3=ur0*ir3;
            PrecisionType ur3ur2=ur3+kappa*ur2;

            VectorT<PrecisionType> gd = -dx*gr3;
            VectorT<PrecisionType> ud = -dx*ur3ur2;

            A += wt*(gr1 - ur1);
            B += wt * normal.dot(inv_epsilon*gd - ud);
            C += -wt * n0.dot(gd - ud*inv_epsilon);
            if (kappa != 0)
            {
                
                PrecisionType gr5 = gr0*ir5;
                PrecisionType ur4=ur0*ir4;
                PrecisionType ur5=ur0*ir5;
                
                PrecisionType pur4ur3=kappa*(3.0*ur4+kappa*ur3)+3.0*ur5;

                for (int jj=0; jj < 3; ++jj)
                {
                    for (int kk=0; kk < 3; ++kk)
                    {
                        PrecisionType kronecker = (jj == kk) ? 1.0 : 0.0;
                        PrecisionType dxjk = dx(jj)*dx(kk);
                        gdd(jj,kk) = 3.0*dxjk*gr5 - kronecker*gr3;
                        udd(jj,kk) = dxjk*pur4ur3 - kronecker*ur3ur2;                        
                    }
                }
                
                for (int jj=0; jj < 3; ++jj)
                {
                    for (int kk=0; kk < 3; ++kk)
                    {
                        D -= wt * n0(kk) * normal(jj) * (gdd(jj,kk) - udd(jj,kk));
                    }
                }

            }
        }
        
        PrecisionType qual_wt = src_it->weight();
        A_out += A * qual_wt; 
        B_out += B * qual_wt;
        C_out += C * qual_wt;
        D_out += D * qual_wt;
    }
    
    D_out *= inv_epsilon;
    

}

void FMM_results_to_BEM(double& f_result,
                            double& h_result,
                            const Vector& normal,
                            DblMatrix1D& potentials,
                            DblMatrix2D& fields,
                            DblMatrix3D& grad_fields);


double fGeometricCorrection(double geometry_correction, double epsilon);
double hGeometricCorrection(double geometry_correction, double epsilon);

double Gpt(const Vector& source_point, const Vector& target_point, double kappa=0.0);
double upt(const Vector& source_point, const Vector& target_point, double kappa);
double dGpt_dn(const Vector& source_point,
                    const Vector& target_point,
                    const Vector& normal);
double dupt_dn(const Vector& source_point,
                    const Vector& target_point,
                    const Vector& normal,
                    double kappa);
double d2Gpt_dndn0_minus_d2upt_dndn0_small_r(const Vector& source_point,
                                                const Vector& target_point,
                                                const Vector& n0,
                                                const Vector& n,
                                                double kappa);
double d2Gpt_dndn0_minus_d2upt_dndn0(const Vector& source_point,
                                        const Vector& target_point,
                                        const Vector& n0,
                                        const Vector& n,
                                        double kappa);

#endif /* BEM_KERNELS_H_ */

