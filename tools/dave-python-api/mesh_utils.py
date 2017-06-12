#!/usr/bin/python
# -*- coding: utf-8 -*-

from mesh import CPP_Mesh
from string import atof
import gts_utils
import kintools

def dump_xyzn(gts_filename, xyzn_filename):

    mesh = CPP_Mesh()

    vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    # add vertices to the mesh (as Vertex objects)
    for v in vertices:
        mesh.add_vertex(v)

    # define the triangles in the mesh
    num_triangles = len(triangles)
    for i,t in enumerate(triangles):
        mesh.define_triangle(t[0], t[1], t[2]);

    # init NodePatch objects
    mesh.create_node_patches()

    fout = open(xyzn_filename,'w')
    for np in mesh.node_patches:
        x,y,z = np.node()
        xn,yn,zn = np.normal()
        print >>fout, "%f %f %f %f %f %f" %(x,y,z, xn,yn,zn)
    fout.close()

    return

def colour_mesh_results(gts_filename, results_filename, kin_filename):

    mesh = CPP_Mesh()

    vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    # add vertices to the mesh (as Vertex objects)
    for v in vertices:
        mesh.add_vertex(v)

    # define the triangles in the mesh
    num_triangles = len(triangles)
    for i,t in enumerate(triangles):
        mesh.define_triangle(t[0], t[1], t[2]);

    results = open(results_filename, 'r').readlines()
    assert(len(results) == mesh.num_vertices)

    for line, v in zip(results, mesh.vertices):
        f,h = [atof(xx) for xx in line.split()]
        v.f = f
        v.h = h

    max_f = max([abs(v.f) for v in mesh.vertices])
    max_h = max([abs(v.h) for v in mesh.vertices])

    from math import log10, ceil
    scale_f = 10 ** ceil(log10(max_f))
    scale_h = 10 ** ceil(log10(max_h))

    print "f/h scale factors: %f / %f (max abs f/h vals: %f / %f)" %(scale_f, scale_h, max_f, max_h)

    kin_out = open(kin_filename, 'w')
    print >>kin_out, "@kinemage\n"
    print >>kin_out, kintools.red_blue_colours()
    print >>kin_out, mesh.kinemage_fh_vals(1e-3, 1e-3, 100)
    kin_out.close()

def precalc_energy_vector(gts_filename, pqr_filename, energy_filename):

    mesh = CPP_Mesh()

    vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    # add vertices to the mesh (as Vertex objects)
    for v in vertices:
        mesh.add_vertex(v)

    # define the triangles in the mesh
    num_triangles = len(triangles)
    for i,t in enumerate(triangles):
        mesh.define_triangle(t[0], t[1], t[2]);

    # add charges
    import pqrtools
    import _BEM
    cpp_charges = [_BEM.Charge(c.charge,c.position)
                    for c in pqrtools.getCharges(pqr_filename)]
    self.add_charges(cpp_charges)


    fvals = [0.0 for xx in range(len(vertices))]
    hvals = [0.0 for xx in range(len(vertices))]

    mesh.precalc_energy_vector(fvals, hvals)

    fout = open(energy_filename, 'w')
    for ff,hh in zip(fvals, hvals):
        print >>fout, "%f %f" %(ff, hh)
    fout.close()

def calc_energy(energy_precalc_filename, fh_filename, epsilon_ratio):
    """Calculate the energy using a results file and a pre-calc file."""
    from string import atof

    energy_precalc_lines = open(energy_precalc_filename,'r').readlines()
    fh_val_lines = open(fh_filename, 'r').readlines()

    assert(len(energy_precalc_lines) == len(fh_val_lines))

    energy = 0.0
    for e_line, fh_line in zip(energy_precalc_lines, fh_val_lines):
        ef,eh = [atof(xx) for xx in e_line.split()]
        f,h = [atof(xx) for xx in fh_line.split()]
        energy += f*ef + h*eh*epsilon_ratio;

    import constants
    con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
      / (constants.epsilon0 * constants.Angstroms * 1000.0)
    print "Energy is %f" %(energy * con)

    return

if __name__ == "__main__":
    import sys
    from string import atof, atoi
    executable = sys.argv[0].split("/")[-1]
    if (executable == "gts2xyzn"):
        dump_xyzn(sys.argv[1], sys.argv[2])
    elif (executable == "colour_results"):
        gts_filename, results_filename, kin_filename = sys.argv[1:]
        colour_mesh_results(gts_filename, results_filename, kin_filename)
    elif (executable == "precalc_energy"):
        gts_filename, pqr_filename, energy_filename = sys.argv[1:]
        precalc_energy_vector(gts_filename, pqr_filename, energy_filename)
    elif (executable == "evaluate_energy"):
        energy_filename, fh_filename, eps_str = sys.argv[1:]
        calc_energy(energy_filename, fh_filename, atof(eps_str))

