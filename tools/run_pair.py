#!/usr/bin/python
# -*- coding: utf-8 -*-
from pybeep import BEEP, Mesh, Vector, Quaternion
from geometry import _rand_rot
from constants import calculate_kappa
import sys

#running_output = open("ache-fas-pair.txt",'w')
#running_output.close()

Dsolvent = 80.0
Dprotein = 2.0
kappa = float(sys.argv[1])
GMRES_tolerance = 1e-6
GMRES_max_iterations = 100

no_rotation = Quaternion(1,0,0,0)
centre_crystal_ache = Vector(35.3081,17.9041,169.832)
centre_crystal_fas = Vector(6.91911,25.184,168.836)
axis = (centre_crystal_fas - centre_crystal_ache)
crystal_separation = axis.length() 
print(crystal_separation)
axis.normalise()
minimum_separation = 10.165
ache_location = Vector(0,0,0)

ache = "ache-fas-png1/1MAH-ache.7.mtz"
fas = "ache-fas-png1/1MAH-fas.10.mtz"

qual_pts = 0
quad_pts = 0
nbsize = 2200
for dummy in [0]:
    beep = BEEP(Dsolvent, kappa, quad_pts, qual_pts, nbsize)
    beep.load_library_mesh(ache)
    beep.load_library_mesh(fas)

    results_file = "results-salt-%f.dat" %(kappa)
    results = open(results_file,'a')
    print("Centre-centre separation, Total Energy" , file=results)
    results.close()

    for iteration in range(0,100,2):
        
        fas_location = ache_location + axis*(crystal_separation+minimum_separation + iteration)
        beep.clear_mesh_instances()
        beep.insert_mesh_instance(0, ache_location, no_rotation, Dprotein)
        beep.insert_mesh_instance(1, fas_location, no_rotation, Dprotein)
        beep.solve(GMRES_tolerance, GMRES_max_iterations)
        beep.reset_library_fh_vals();
        total_energy = beep.calculate_energies()

        results = open(results_file,'a')
        print("%f %f" %(minimum_separation + iteration, total_energy), file=results)
        results.close()
