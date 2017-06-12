#!/usr/bin/env python3

import constants
from math import pi

def born_solvation_energy(charge, radius, Dint, Dext):
    """The analytic born solvation energy in kJ/mol (charge in e, radius in Angstroms)"""

    prefac = (charge*charge) / (8.0*pi*radius)
    energy = -(1.0/Dint - 1.0/Dext)*prefac
    return energy * constants.convert_energy_to_kj_per_mol

if __name__=="__main__":
    import sys
    argc = len(sys.argv)
    assert(argc >= 3 and argc <= 5)

    charge = float(sys.argv[1])
    radius = float(sys.argv[2])
    Dint = 2.0
    Dext = 80.0
    
    if (argc >= 4): Dint=float(sys.argv[3])
    if (argc == 5): Dext=float(sys.argv[4])
    
    e = born_solvation_energy(charge, radius, Dint, Dext)
    print(e, "kJ/mol")

