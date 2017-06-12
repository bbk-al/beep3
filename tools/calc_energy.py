#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pybeep import *

if __name__=="__main__":

    import sys
    from math import sqrt
    from random import uniform
    import constants

    mesh_filename = sys.argv[1]
    fh_filename = sys.argv[2]
    m = Mesh(mesh_filename)
    #m.init_energy_precalcs()
    fh_vals = [xx.split()[:2] for xx in open(fh_filename, 'r').readlines()]
    fh_vals = [(float(ff),float(hh)) for (ff,hh) in fh_vals]

    eps_int = 2.0
    eps_ext = 80.0
    epsilon_ratio = eps_ext/eps_int
    qe_force = Vector(0,0,0)
    energy = 0.0
    kahan = 0.0
    kahan_force = Vector(0,0,0)
    for (f,h),np in zip(fh_vals, m.node_patches):

        # add energy using kahan compensated addition
        energy_frag = (np.energy_coefficient_h*epsilon_ratio*h - np.energy_coefficient_f*f);
        y = energy_frag - kahan;
        t = energy + y;
        kahan = (t - energy) - y;
        energy = t;

        # add force using kahan compensated addition
        qe_force_frag = (np.force_coefficient_h*epsilon_ratio*h - np.force_coefficient_f*f);
        y_vec = qe_force_frag - kahan_force
        t_vec = qe_force + y_vec
        kahan_force = (t_vec - qe_force) - y_vec
        qe_force = t_vec

    print(energy * 0.5 * constants.convert_energy_to_kj_per_mol)

