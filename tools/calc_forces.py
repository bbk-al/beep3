#!/usr/bin/python
# -*- coding: utf-8 -*-
from libBEEP import *

def test_force_on_planar_triangle():
    v1 = Vector(0,0,0)
    v2 = Vector(1,0,0)
    v3 = Vector(1,1,0)
    n1 = Vector(0,0,1)
    n2 = Vector(0,0,1)
    n3 = Vector(0,0,1)
    h=HybridBezierTriangle(v1,n1,v2,n2,v3,n3)
    print h.calculate_force(Vector(0,0,0),Vector(1,1,1),0)

if __name__=="__main__":

    import sys
    from mesh import Mesh
    from math import sqrt
    from random import uniform
    import constants

    mesh_filename = sys.argv[1]
    fh_filename = sys.argv[2]
    m = Mesh(mesh_filename)
    fh_vals = [xx.split()[:2] for xx in open(fh_filename, 'r').readlines()]
    fh_vals = [(float(ff),float(hh)) for (ff,hh) in fh_vals]

    eps_int = 4.0
    eps_ext = 80.0
    kappa = 3.0
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
    
    qe_force *= constants.force_conversion_picoNewtons

    print "Total energy: ", energy * 0.5 * constants.convert_energy_to_kj_per_mol
    print "Total qE force: ", qe_force

    def trial():
        #normalisation = sqrt(k_para**2 + k_para**2 + k_perp**2)
        #k_para /= normalisation
        #k_perp /= normalisation
        force = Vector(0,0,0)
        force_kahan = Vector(0,0,0)
        for ctr,tri in enumerate(m.triangles):
            v1 = m.get_vertex(tri.v1_idx())
            n1 = v1.normal()
            v1 = v1.vertex()
            f1,h1 = fh_vals[tri.v1_idx()]

            v2 = m.get_vertex(tri.v2_idx())
            n2 = v2.normal()
            v2 = v2.vertex()
            f2,h2 = fh_vals[tri.v2_idx()]

            v3 = m.get_vertex(tri.v3_idx())
            n3 = v3.normal()
            v3 = v3.vertex()
            f3,h3 = fh_vals[tri.v3_idx()]

            n0 = tri.normal
            #h=HybridBezierTriangle(tri)
            h=PNG1_Triangle(tri)
            #force += h.calculate_force(Vector(f1,f2,f3),Vector(h1,h2,h3),eps_int, eps_ext,0,2)
            force_frag = tri.calculate_force(Vector(f1,f2,f3),Vector(h1,h2,h3),eps_int, eps_ext,kappa,0)
            y_vec = force_frag - force_kahan
            t_vec = force + y_vec
            force_kahan = (t_vec - force) - y_vec
            force = t_vec

        return force 
    #while True:
    #deltas = [0,0,prev_winner]
    #for qn, choice in enumerate(choices):
    #    force = trial(scale_para*choice[0], scale_perp*choice[1])
    #    print force, force.length()
    #    deltas[qn] = force.length()

    #winner = 2
    #for n,r in enumerate(deltas):
    #    if r < deltas[winner]:
        #        winner = n
    #	prev_winner = r
    #if winner == 2:
    #    delta /= 2.0
    #else:
    #    scale_para *= choices[winner][0]
    #    scale_perp *= choices[winner][1]
    #    normalisation = sqrt(scale_para**2 + scale_para**2 + scale_perp**2)
    #    scale_para /= normalisation
    #    scale_perp /= normalisation

    force = trial() * constants.force_conversion_picoNewtons
 
    print "boundary force: ", force
    force += qe_force
    print force, force.length()
    #break


