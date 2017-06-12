/* Test program:  get target values for gsl_sf_legendre_deriv_array_e */
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
    double *result_array = new double[gsl_sf_legendre_array_n(nmax)];
    double *result_deriv_array = new double[gsl_sf_legendre_array_n(nmax)];

	gsl_sf_legendre_deriv_array_e(GSL_SF_LEGENDRE_SPHARM, nmax, cos_theta,
									-1, result_array, result_deriv_array);
	for (int nctr=0; nctr <= nmax; ++nctr) {
        std::cout << nctr << " " << 0 << " "
				<< result_deriv_array[gsl_sf_legendre_array_index(nctr,0)]
								* root2 * -sin_theta
				<< std::endl;
	}
	for (int m=1; m <= nmax; ++m) {
		for (int nctr=m; nctr <= nmax; ++nctr) {
			std::cout << nctr << " " << m << " "
					<< result_deriv_array[gsl_sf_legendre_array_index(nctr,m)]
								* root2 * cos(phi*m) * -sin_theta
					<< std::endl;
			std::cout << nctr << " " << -m << " "
					<< result_deriv_array[gsl_sf_legendre_array_index(nctr,m)]
								* root2 * sin(phi*m) * -sin_theta
					<< std::endl;
		}
	}

	delete[] result_array;
	delete[] result_deriv_array;

	return 0;
}
