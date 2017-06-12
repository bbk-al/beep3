/*
* bem_kernels.cpp
*
*  Created on: 22 Jul 2010
*      Author: david
*/

#include "bem_kernels.h"
#include "node_patch.h"
#include "quad_point.h"

void FMM_results_to_BEM(double& f_result,
                            double& h_result,
                            const Vector& normal,
                            DblMatrix1D& potentials,
                            DblMatrix2D& fields,
                            DblMatrix3D& grad_fields)
{
#if 1
    double f_accum = 0.0;
    double h_accum = 0.0;
    // kappa = actual vals
    {
        // fcomps
        double f_nx_field_x = fields(0,0);
        double f_ny_field_y = fields(1,1);
        double f_nz_field_z = fields(2,2);
        f_accum += (f_nx_field_x + f_ny_field_y + f_nz_field_z);

        // hcomps
        double hcomps = 0.0;
        for (int jj=0; jj < 3; ++jj)
        {
            hcomps += (grad_fields(jj,0,4)+grad_fields(0,jj,4)) * normal(jj);
            hcomps += (grad_fields(jj,1,5)+grad_fields(1,jj,5)) * normal(jj);
            hcomps += (grad_fields(jj,2,6)+grad_fields(2,jj,6)) * normal(jj);
        }
        h_accum += hcomps / 2.0;

        f_accum += potentials(3);
        h_accum += (fields(0,7)*normal.x + fields(1,7)*normal.y + fields(2,7)*normal.z);
    }

    // kappa = kappa0 (1e-10) vals
    {
        // fcomps
        double f_nx_field_x = fields(0,8);
        double f_ny_field_y = fields(1,9);
        double f_nz_field_z = fields(2,10);
        f_accum -= f_nx_field_x + f_ny_field_y + f_nz_field_z;

        // hcomps
        double hcomps = 0.0;
        for (int jj=0; jj < 3; ++jj)
        {
            hcomps += (grad_fields(jj,0,8)+grad_fields(0,jj,8)) * normal(jj);
            hcomps += (grad_fields(jj,1,9)+grad_fields(1,jj,9)) * normal(jj);
            hcomps += (grad_fields(jj,2,10)+grad_fields(2,jj,10)) * normal(jj);
        }
        h_accum -= hcomps / 2.0;

        f_accum -= potentials(11);
        h_accum -= (fields(0,11)*normal.x + fields(1,11)*normal.y + fields(2,11)*normal.z);
    }

    f_result += f_accum * ONE_OVER_4PI;
    h_result += h_accum * ONE_OVER_4PI;
    return;
#else
    //std::cout << f_accum << " " << h_accum << std::endl;

    double f_accum = 0.0;
    double h_accum = 0.0;
    // kappa = actual vals
    {
        // fcomps
        double f_nx_field_x = fields(0,0);
        double f_ny_field_y = fields(1,1);
        double f_nz_field_z = fields(2,2);
        f_accum -= (f_nx_field_x + f_ny_field_y + f_nz_field_z);

        // hcomps
        double hcomps = 0;
        for (int jj=0; jj < 3; ++jj)
        {
            hcomps += normal(jj) * (grad_fields(0,jj,4) + grad_fields(1,jj,5) + grad_fields(2,jj,6));
        }
        h_accum += hcomps;
        
        f_accum += potentials(3);
        h_accum -= (fields(0,7)*normal.x + fields(1,7)*normal.y + fields(2,7)*normal.z);
    }
 
    // kappa = kappa0 (1e-10) vals
    {
        // fcomps
        double f_nx_field_x = fields(0,8);
        double f_ny_field_y = fields(1,9);
        double f_nz_field_z = fields(2,10);
        f_accum += f_nx_field_x + f_ny_field_y + f_nz_field_z;

        // hcomps
        double hcomps = 0;
        for (int jj=0; jj < 3; ++jj)
        {
            hcomps += normal(jj) * (grad_fields(0,jj,8) + grad_fields(1,jj,9) + grad_fields(2,jj,10));
        }
        h_accum -= hcomps;

        f_accum -= potentials(11);
        h_accum += (fields(0,11)*normal.x + fields(1,11)*normal.y + fields(2,11)*normal.z);
    }

    f_result += f_accum * ONE_OVER_4PI;
    h_result += h_accum * ONE_OVER_4PI;

    //std::cout << f_accum << " " << h_accum << std::endl;

    return;
#endif
}

double fGeometricCorrection(double geometry_correction, double epsilon)
{
    //assert(geometry_correction == 0.5);
    double cg = 1.0 - geometry_correction;
    return cg + (geometry_correction / epsilon);

}

double hGeometricCorrection(double geometry_correction, double epsilon)
{
    //assert(geometry_correction == 0.5);
    double cg = 1.0 - geometry_correction;
    return geometry_correction + (cg / epsilon);
}

double Gpt(const Vector& source_point, const Vector& target_point, double kappa)
{
    double r = distance(source_point, target_point);
    return ONE_OVER_4PI / r;
}

double upt(const Vector& source_point, const Vector& target_point, double kappa)
{
    double r = distance(source_point, target_point);
    return EXPONENTIAL(-kappa * r) * ONE_OVER_4PI / r;
}

double dGpt_dn(const Vector& source_point,
                    const Vector& target_point,
                    const Vector& normal)
{
    // analytical method (Juffer et al.)
    Vector diff = subtract(target_point, source_point);
    double ct = cosTheta(normal, diff);
    double r2 = distance2(source_point, target_point);
    double analytical = -ct * ONE_OVER_4PI / r2;
    //std::cout << " src: " << source_point << " targ: " << target_point << " diff: " << diff << " ct: "<< ct << " r2: " << r2 << " analytical: " << analytical << std::endl;

    return analytical;
}

double dupt_dn(const Vector& source_point,
                    const Vector& target_point,
                    const Vector& normal,
                    double kappa)
{
    // if kappa is zero then dGpt_dn is the same thing, and avoids a few flops.
    if (kappa==0.0f) {
        return dGpt_dn(source_point, target_point, normal);
    }

    // analytical method (Juffer et al.)
    double r2 = distance2(source_point, target_point);
    double kappa_r = kappa * SQUARE_ROOT(r2);
    double lhs = (1.0f + kappa_r) * exp(-kappa_r);

    Vector diff = subtract(target_point, source_point);
    double ct = cosTheta(normal, diff);
    double rhs = ct * ONE_OVER_4PI / r2;
    double analytical = -lhs*rhs;

    return analytical;
}

double d2Gpt_dndn0_minus_d2upt_dndn0_small_r(const Vector& source_point,
                                                const Vector& target_point,
                                                const Vector& n0,
                                                const Vector& n,
                                                        double kappa)
{
    // see Juffer et al. pg 151, eqn 2.15 (g(r;s))
    Vector r_s = subtract(target_point, source_point);
    double r2 = distance2(source_point, target_point);

    double ct = cosTheta(n ,r_s);
    double ct0 = cosTheta(n0, r_s);

    double kap2_pi = kappa*kappa / pi;
    double ndotn = n0.dot(n);

    double r_inv = R_SQUARE_ROOT(r2);
    double left = (ct*ct0 - ndotn) * r_inv / 8.0f;
    double right = kappa * ndotn / 12.0f;

    double result = kap2_pi*(left + right);
    return result;

}

double d2Gpt_dndn0_minus_d2upt_dndn0(const Vector& source_point,
                                        const Vector& target_point,
                                        const Vector& n0,
                                        const Vector& n,
                                        double kappa)
{
    // see Juffer et al. pg 150, eqn 2.14+
    double r2 = distance2(source_point, target_point);
    if (r2 < 1e-2f) {
        return d2Gpt_dndn0_minus_d2upt_dndn0_small_r(source_point, target_point, n0, n, kappa);
    }
    Vector r_s = subtract(target_point, source_point);
    double ct0ct = cosTheta(n ,r_s) * cosTheta(n0, r_s);
    double r = SQUARE_ROOT(r2);

    double r3 = r2 * r;
    double kappa_r = kappa * r;
    double kappa2 = kappa*kappa;
    double ndotn = n0.dot(n);

    double d2fdn0dn = (ndotn - 3.0f*ct0ct) * ONE_OVER_4PI / r3;
    double ekr = expf(-kappa_r);
    double left = (1.0f + kappa_r)* d2fdn0dn;
    double right = (kappa2 * ct0ct) * ONE_OVER_4PI / r;

    double result = d2fdn0dn - ekr*(left - right);

    return result;

}

// weights and abcissae for 10pt Gaussian quadrature
//(+/- these quad points symmetrical about the middle of the integration range)
double gauss_quads[5] = {0.1488743389816312108848260f,
                        0.4333953941292471907992659f,
                        0.6794095682990244062343274f,
                        0.8650633666889845107320967f,
                        0.9739065285171717200779640f};
double gauss_weights[5] = {0.2955242247147528701738930f,
                        0.2692667193099963550912269f,
                        0.2190863625159820439955349f,
                        0.1494513491505805931457763f,
                        0.0666713443086881375935688f};

double
A_integrand_kernel_rho(double rho, double kappa, double param)
{
    double t = param;
    double r = 2.0*t*(t-1.0) + 1.0;
    return (1.0 - expf(-kappa * rho * r)) / (4.0 * pi * r);
}

double
gaussian_quadrature_10_rho(double kappa, double param, double a, double b)
{
    // See Numerical Recipes 3rd edition 4.6 pp180.
    double range = (b-a)*0.5f;
    double midpt = (b+a)*0.5f;
    double dx;
    double accum=0.0f;
    for (unsigned int i=0; i < 5; ++i)
    {
        dx = range*gauss_quads[i];
        accum += gauss_weights[i]*(A_integrand_kernel_rho(midpt+dx, kappa, param) +
                                A_integrand_kernel_rho(midpt-dx, kappa, param));
    }
    return accum*range;
}

double
A_integrand_kernel_t(double t, double kappa)
{

    // integrate rho between 0 and 1
    double result = gaussian_quadrature_10_rho(kappa,
                                        t,
                                        0.0f,
                                        1.0f);

    return result;
}

double
gaussian_quadrature_10_t(double kappa, double a, double b)
{
    // See Numerical Recipes 3rd edition 4.6 pp180.
    double range = (b-a)*0.5f;
    double midpt = (b+a)*0.5f;
    double dx;
    double accum=0.0f;
    for (unsigned int i=0; i < 5; ++i)
    {
        dx = range*gauss_quads[i];
        accum += gauss_weights[i]*(A_integrand_kernel_t(midpt+dx, kappa) +
                                A_integrand_kernel_t(midpt-dx, kappa));
    }
    return accum*range;
}
