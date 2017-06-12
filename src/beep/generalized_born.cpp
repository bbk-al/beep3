/*
* generalized_born.cpp
*
*  Created on: 15 Mar 2011
*      Author: david
*/

#include "generalized_born.h"

int main(int argc, char *argv[])
{

    // Usage: generalized_born gtsfile chargesfile internal_dielectric external_dielectric

    // get input vars
    std::string gts_filename(argv[1]);
    std::string charges_filename(argv[2]);
    std::string eps_int_str(argv[3]);
    std::string eps_ext_str(argv[4]);

    double eps_int = ::atof(eps_int_str.c_str());
    double eps_ext = ::atof(eps_ext_str.c_str());

    std::cout << "Reading input files..." << std::endl;
    GeneralizedBorn gb(gts_filename, charges_filename);

    std::cout << "Integrating over surface to find Born radii..." << std::endl;
    gb.calculate_born_radii();

    std::cout << "Evaluating energy..." << std::endl;
    double e = gb.evaluate_energy(eps_int, eps_ext);
    std::cout << "GB Energy = " << e*beep_constants::convert_energy_to_kj_per_mol << " kJ/mol" << std::endl;

    return 0;

}
