/*
 * constants.h
 *
 *  Created on: 30 Jul 2010
 *      Author: david
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace beep_constants
{

static const double Avogadro = 6.02214179e23;
static const double Boltzmann_SI = 1.3806503e-23; // k in SI units (Joules per Kelvin)
static const double Temperature_in_Kelvin = 300.0; // T

// kT -- boltzmann constant times temperature
static const double kT = Boltzmann_SI * Temperature_in_Kelvin; // Joules

// permittivity of free space is in Farads per metre
static const double epsilon0 = 8.85418782e-12;
static const double elementary_charge = 1.60217646e-19;  // e in coulombs

static const double Angstroms = 1.0e-10;  // metres per Angstrom

// Electric field (in internal units) is in elementary charge units per
// Angstrom squared, and has only been scaled by the relative permittivity not
// the real permittivity of free space. E = - q / (4pi_eps_r*r)
static const double convert_E_field_to_SI = elementary_charge / (epsilon0 * Angstroms * Angstroms);

// the Maxwell stress tensor elements are in units of electric field elements,
// squared, multiplied by the permittivity.
static const double convert_MST_to_SI = (convert_E_field_to_SI*convert_E_field_to_SI) * epsilon0;

// stress tensor (force per unit area) is multiplied by area to yield force -->
// convert from square Angstroms to square metres through two multiplications
// by metres per Angstrom.  This gives force in Newtons (= Joules per metre)
static const double force_conversion_factor_SI = convert_MST_to_SI * Angstroms * Angstroms;
static const double force_conversion_Newtons = force_conversion_factor_SI;
static const double convert_force_to_picoNewtons = 1.0e12 * elementary_charge*elementary_charge/(epsilon0 * Angstroms * Angstroms);
static const double convert_force_to_kJ_per_mol_Angstrom = elementary_charge*elementary_charge * Avogadro /(1000. * epsilon0 * Angstroms);

static const double kJ_in_kcals = 0.239005736;

// Newtons are Joules per metre; we want kT per Angstrom, so divide
// by kT in Joules, and multiply by metres per Angstrom
static const double force_conversion_kT_per_Angstrom=force_conversion_factor_SI * Angstroms /kT;

// Convert energy from angstroms, e, without epsilon0, to kj/mol
static const double convert_energy_to_kj_per_mol = Avogadro * elementary_charge * elementary_charge / (epsilon0 * Angstroms * 1000.0);
static const double convert_potential_to_kT_per_e = elementary_charge * elementary_charge / (epsilon0 * Angstroms * kT);

static double convert_monovalent_conc_to_kappa(double concentration, double solvent_dielectric)
{
	double screening_length = sqrt(solvent_dielectric*epsilon0*kT/(2.0*Avogadro*elementary_charge*elementary_charge*concentration)) / Angstroms;
	return 1.0 / screening_length;
}

}

#endif /* CONSTANTS_H_ */
