/* Test program:  get baseline values for gsl_sf_legendre_sphPlm_deriv_array */
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

const double root2 = 1.4142135623730951;

int main(int argc, const char *argv[]) {
	// Ynm_dtheta
	// Test values
	int nmax = 5;
	double cos_theta = 0.5;
	double sin_theta = sqrt(3)/2;
	double phi = M_PI/5;

    assert(nmax >= 0);
    double *result_array = new double[(nmax+1)*(nmax+1)];
    double *result_deriv_array = new double[(nmax+1)*(nmax+1)];

    for (int m=-nmax; m <= nmax; ++m)
    {
        if (m >= 0)
        {
            gsl_sf_legendre_sphPlm_deriv_array(nmax, m, cos_theta,
											result_array, result_deriv_array);
            int res_ctr=0;
            for (int nctr=fabs(m); nctr <= nmax; ++nctr)
            {
                std::cout << nctr << " " << m << " "
						<< result_deriv_array[res_ctr++] * root2
							* cos(phi*fabs(m)) * -sin_theta
						<< std::endl;
            }
        }
        else
        {
            gsl_sf_legendre_sphPlm_deriv_array(nmax, -m, cos_theta,
											result_array, result_deriv_array);
            int res_ctr=0;
            for (int nctr=fabs(m); nctr <= nmax; ++nctr)
            {
                std::cout << nctr << " " << m << " "
						<< result_deriv_array[res_ctr++] * root2
							* sin(phi*fabs(m)) * -sin_theta
						<< std::endl;
            }
        }
    }

    delete[] result_array;
    delete[] result_deriv_array;

	return 0;
}
