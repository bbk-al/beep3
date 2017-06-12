# vim: set fileencoding=UTF8 :
#
# constants.py
#
# Useful physical constants and conversions. Internal units are generally
# Angstroms, elementary charge, and permittivities are relative.
#
# Author: David Fallaize, University College London, 2009
# E-mail: drf33@cantab.net
#

# pi is fairly useful...
from math import pi

Avogadro = 6.02214179e23
Boltzmann_SI = 1.3806503e-23 # k in SI units (Joules per Kelvin)
Temperature_in_Kelvin = 298.15 # T

# kT -- boltzmann constant times temperature
kT = Boltzmann_SI * Temperature_in_Kelvin # Joules

# permittivity of free space is in Farads per metre
epsilon0 = 8.854e-12
elementary_charge = 1.602e-19 # e in coulombs

Angstroms = 1.0e-10 # metres per Angstrom

# Electric field (in internal units) is in elementary charge units per
# Angstrom squared, and has only been scaled by the relative permittivity not
# the real permittivity of free space. E = - q / (4pi_eps_r*r)
convert_E_to_SI = elementary_charge / (epsilon0 * Angstroms * Angstroms)

# the Maxwell stress tensor elements are in units of electric field elements,
# squared, multiplied by the permittivity.
convert_MST_to_SI = (convert_E_to_SI ** 2) * epsilon0 

# stress tensor (force per unit area) is multiplied by area to yield force -->
# convert from square Angstroms to square metres through two multiplications
# by metres per Angstrom.  This gives force in Newtons (= Joules per metre)
force_conversion_factor_SI = convert_MST_to_SI * Angstroms * Angstroms
force_conversion_Newtons = force_conversion_factor_SI

# Newtons are Joules per metre; we want kT per Angstrom, so divide
# by kT in Joules, and multiply by metres per Angstrom
#force_conversion_kT_per_Angstrom=force_conversion_factor_SI * Angstroms /kT
force_conversion_kT_per_Angstrom = 7045.0738534336688  # 1 Newton in units kT/A

#
# TESTS
#
if __name__=="__main__":
    
    # check that the cached version of conversion factor is correct
    assert(force_conversion_kT_per_Angstrom == force_conversion_factor_SI * Angstroms /kT)

    from math import pi
    
    # force between two like charges, 1 Angstrom apart, in SI units (Newtons)
    actual = elementary_charge * elementary_charge / (4.0*pi*epsilon0*1.0e-10*1.0e-10)
    
    # same thing expressed in internal units (unscaled permittivity,
    # Angstroms, elementary charges)
    internal = 1.0*1.0 / (4.0*pi*1.0*1.0*1.0)

    # internal force units converted to SI units
    converted = internal*force_conversion_Newtons
    
    # check that our conversion factor gets it right to good precision
    assert(abs(actual - converted)/actual < 1e-15)
