#!/usr/bin/python
# -*- coding: utf-8 -*-

from pybeep import Mesh
import gts_utils
import kintools

def dump_xyzn(gts_filename, xyzn_filename):

    mesh = Mesh()

    vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    # add vertices to the mesh (as Vertex objects)
    for v in vertices:
        mesh.add_vertex(v)

    # define the triangles in the mesh
    num_triangles = len(triangles)
    for i,t in enumerate(triangles):
        mesh.define_triangle(t[0], t[1], t[2]);

    # init NodePatch objects
    #mesh.create_node_patches()

    fout = open(xyzn_filename,'w')
    for np in mesh.node_patches:
        x,y,z = np.node()
        xn,yn,zn = np.normal()
        print("%f %f %f %f %f %f" %(x,y,z, xn,yn,zn), file=fout)
    fout.close()

    return

def colour_delta_results(mtz_filename, base_results_filename, results_filename, kin_filename, scale_f=0, scale_h=0):
    
    try:
        scale_f = float(scale_f)
    except:
        pass
    try:
        scale_h = float(scale_h)
    except:
        pass
    
    mesh = Mesh(mtz_filename)
    #vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    base_results = open(base_results_filename, 'r').readlines()
    results = open(results_filename, 'r').readlines()
    assert(len(results) == mesh.num_vertices)

    for line_base, line, v in zip(base_results, results, mesh.node_patches):
        f0,h0 = [float(xx) for xx in line_base.split()]
        f,h = [float(xx) for xx in line.split()]
        v.f = f - f0
        v.h = h - h0

    max_f = max([abs(v.f) for v in mesh.node_patches])
    max_h = max([abs(v.h) for v in mesh.node_patches])

    #from math import log10, ceil
    #scale_f = 10 ** ceil(log10(max_f))
    #scale_h = 10 ** ceil(log10(max_h))
    if (scale_f == 0): scale_f = max_f
    if (scale_h == 0): scale_h = max_h

    print("f/h scale factors: %f / %f (max abs f/h vals: %f / %f)" %(scale_f, scale_h, max_f, max_h))

    kin_out = open(kin_filename, 'w')
    print(kin_out, "@kinemage\n")
    print(kintools.hundred_red_blue_colours(), file=kin_out)
    print(mesh.kinemage_fh_vals(scale_f, scale_h, 100), file=kin_out)
    kin_out.close()
    
def colour_mesh_results(mtz_filename, results_filename, kin_filename, scale_f=0, scale_h=0):
    
    try:
        scale_f = float(scale_f)
    except:
        pass
    try:
        scale_h = float(scale_h)
    except:
        pass
    
    mesh = Mesh(mtz_filename)
    #vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    results = open(results_filename, 'r').readlines()
    assert(len(results) == mesh.num_vertices)

    for line, v in zip(results, mesh.node_patches):
        f,h = [float(xx) for xx in line.split()]
        v.f = f
        v.h = h

    max_f = max([abs(v.f) for v in mesh.node_patches])
    max_h = max([abs(v.h) for v in mesh.node_patches])

    #from math import log10, ceil
    #scale_f = 10 ** ceil(log10(max_f))
    #scale_h = 10 ** ceil(log10(max_h))
    if (scale_f == 0): scale_f = max_f
    if (scale_h == 0): scale_h = max_h

    print("f/h scale factors: %f / %f (max abs f/h vals: %f / %f)" %(scale_f, scale_h, max_f, max_h))

    kin_out = open(kin_filename, 'w')
    print("@kinemage\n", file=kin_out)
    print(kintools.hundred_red_blue_colours(), file=kin_out)
    print(mesh.kinemage_fh_vals(scale_f, scale_h, 100), file=kin_out)
    kin_out.close()
    return

def precalc_energy_vector(gts_filename, xyzq_filename, energy_filename):

    mesh = Mesh()

    vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

    # add vertices to the mesh (as Vertex objects)
    for v in vertices:
        mesh.add_vertex(v)

    # define the triangles in the mesh
    num_triangles = len(triangles)
    for i,t in enumerate(triangles):
        mesh.define_triangle(t[0], t[1], t[2]);

    # add charges
    from _BEM import Charge, Vector
    cpp_charges = []
    for line in open(xyzq_filename,'r').readlines():
        x,y,z,q,r = [float(xx) for xx in line.split()]
        cpp_charges.append(Charge(q, Vector(x,y,z)))
    mesh.add_charges(cpp_charges)

    fvals = [0.0 for xx in range(len(vertices))]
    hvals = [0.0 for xx in range(len(vertices))]

    mesh.precalc_energy_vector(fvals, hvals)

    fout = open(energy_filename, 'w')
    for ff,hh in zip(fvals, hvals):
        print("%e %e" %(ff, hh), file=fout)
    fout.close()
    return

def precalc_peak_splitting_surface_integrals(gts_filename, 
                                             xyzq_filename, 
                                             peak_filename,
                                             kappa,
                                             Dext,
                                             Dint):
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
    from _BEM import Charge, Vector
    cpp_charges = []
    for line in open(xyzq_filename,'r').readlines():
        x,y,z,q,r = [float(xx) for xx in line.split()]
        cpp_charges.append(Charge(q, Vector(x,y,z)))
    mesh.add_charges(cpp_charges)

    fvals = [0.0 for xx in range(len(vertices))]
    hvals = [0.0 for xx in range(len(vertices))]

    #print "Calculating peaks:"
    mesh.precalc_peak_splitting(fvals, hvals, kappa, Dext, Dint)

    #print "Writing peaks to ", peaks_filename
    fout = open(peak_filename, 'w')
    for ff,hh in zip(fvals, hvals):
        print("%e %e" %(ff, hh), file=fout)
    fout.close()
    return

def precalc_rhs(gts_filename, 
                xyzq_filename, 
                rhs_filename,
                Dext):
    
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
    from _BEM import Charge, Vector
    cpp_charges = []
    for line in open(xyzq_filename,'r').readlines():
        x,y,z,q,r = [float(xx) for xx in line.split()]
        cpp_charges.append(Charge(q, Vector(x,y,z)))
    mesh.add_charges(cpp_charges)

    fvals = [0.0 for xx in range(len(vertices))]
    hvals = [0.0 for xx in range(len(vertices))]

    mesh.precalc_rhs(fvals, hvals, Dext)

    fout = open(rhs_filename, 'w')
    for ff,hh in zip(fvals, hvals):
        print("%e %e" %(ff, hh), file=fout)
    fout.close()
    return

def precalc_peak_energy(gts_filename, xyzq_filename, kappa, Dext, Dint):
    
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
    from _BEM import Charge, Vector
    cpp_charges = []
    for line in open(xyzq_filename,'r').readlines():
        x,y,z,q,r = [float(xx) for xx in line.split()]
        cpp_charges.append(Charge(q, Vector(x,y,z)))
    mesh.add_charges(cpp_charges)

    fvals = [0.0 for xx in range(len(vertices))]
    hvals = [0.0 for xx in range(len(vertices))]

    import constants
    con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
      / (constants.epsilon0 * constants.Angstroms * 1000.0)
    energy = mesh.precalc_peak_energy(kappa, Dext, Dint) * con
    return energy

def calc_energy(energy_precalc_filename, fh_filename, epsilon_ratio):
    """Calculate the energy using a results file and a pre-calc file."""

    energy_precalc_lines = open(energy_precalc_filename,'r').readlines()
    fh_val_lines = open(fh_filename, 'r').readlines()

    assert(len(energy_precalc_lines) == len(fh_val_lines))

    energy = 0.0
    for e_line, fh_line in zip(energy_precalc_lines, fh_val_lines):
        ef,eh = [float(xx) for xx in e_line.split()]
        f,h = [float(xx) for xx in fh_line.split()]
        energy += f*ef + h*eh*epsilon_ratio;

    import constants
    con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
      / (constants.epsilon0 * constants.Angstroms * 1000.0)
    print("Energy is %f" %(energy * con))

    return

if __name__ == "__main__":
    import sys
    
    mode = sys.argv[1]
    if (locals().has_key(mode)):
        locals()[mode](*sys.argv[2:])
    else:
        print("Can't run that: try one of these... ", locals())
