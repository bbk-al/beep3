#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from libBEEP import *

if __name__=="__main__":

    import sys
    from mesh import Mesh
    from math import sqrt, pi
    #from random import uniform
    import constants

    mesh_filename = "coulomb.mtz"
    fh_filename = "coulomb-%s.fh" %(sys.argv[1])
    m1 = Mesh(mesh_filename)
    m2 = Mesh(mesh_filename)
    
    eps_int = 80.0
    eps_ext = 80.0
    kappa = 0.0
    sep = float(sys.argv[1])
    analytic = constants.force_conversion_kJ_mol_A / (eps_int*4.0*pi*sep*sep)
       
#    mzero.load_fh_vals("coulomb-iso.fh.0")
    m1.load_fh_vals(fh_filename+".0")
    m2.load_fh_vals(fh_filename+".1")
   
    #energy = m.calculate_energy(kappa, eps_int, eps_ext)
    #print "Total energy: ", energy
#    fzero = mzero.calculate_qe_force(eps_int, eps_ext)
    f1 = m1.calculate_qe_force(eps_int, eps_ext)
    f2 = m2.calculate_qe_force(eps_int, eps_ext)
    total_qe_force = f1+f2
#    print "Total qE by patch-charge: %s (%s / %s)" %(total_qe_force, f1, f2)
    #print f1, f2, analytic
    qe_err= max([100*abs((abs(f1.z) - analytic)/analytic), 100*abs((abs(f2.z) - analytic)/analytic)])
#    bfzero = mzero.calculate_boundary_force(kappa, eps_int, eps_ext) 
    bf1 = m1.calculate_boundary_force(kappa, eps_int, eps_ext)
    bf2 = m2.calculate_boundary_force(kappa, eps_int, eps_ext)
    total_boundary_force = bf1+bf2
    mst_err= max([100*abs((abs(bf1.z) - analytic)/analytic), 100*abs((abs(bf2.z) - analytic)/analytic)])
    print(qe_err, mst_err)
    # print bf1, bf2, bf1+bf2

#print "Total qE by MST: %s (%s / %s)" %(total_boundary_force, bf1, bf2)

