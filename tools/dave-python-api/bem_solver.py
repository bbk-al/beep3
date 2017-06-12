#!/usr/bin/python
# -*- coding: utf-8 -*-
from mesh import *
from scipy.linalg import iterative

GMRES_PRINT_DEBUG = True
GMRES_TOLERANCE = 1.0e-9
GMRES_RESTART = 50

NUM_SPHARM_TERMS = 9

def solveFullElectrostatics(mesh, full_matrix):

    kappa = full_matrix.kappa
    Dext = full_matrix.Dext
    Dint = full_matrix.Dint
    epsilon = Dext / Dint
    num_unknowns = mesh.num_node_patches

    inverse = full_matrix.matrix_inverse

    rhs, fh_peaks = get_fh_peaks(mesh, epsilon, kappa, Dext, Dint)

    solution = inverse * matrix(rhs).T

def E_to_f(E, normal, kappa, fval):

    E2 = E.dot(E)
    Ex = E.x
    Ey = E.y
    Ez = E.z

    #print np.f, fval, np.h, hval, node, normal, E, E.length(), area, E2

    # x stress
    stress_xx = Ex * Ex - 0.5*E2 - 0.5*kappa*kappa*fval*fval
    stress_xy = Ex * Ey
    stress_xz = Ex * Ez

    fx = normal.dot(Vector(stress_xx, stress_xy, stress_xz))

    # y stress
    stress_yx = Ey * Ex
    stress_yy = Ey * Ey - 0.5*E2 - 0.5*kappa*kappa*fval*fval
    stress_yz = Ey * Ez

    fy = normal.dot(Vector(stress_yx, stress_yy, stress_yz))

    # z stress
    stress_zx = Ez * Ex
    stress_zy = Ez * Ey
    stress_zz = Ez * Ez - 0.5*E2 - 0.5*kappa*kappa*fval*fval

    fz = normal.dot(Vector(stress_zx, stress_zy, stress_zz))

    return Vector(fx, fy, fz)

def solveElectrostatics(mesh, kappa, Dext, Dint):
    """Use BEM to solve electrostatic parameters f and h on surface of all
    diffusing entities."""

    # epsilon is ratio of internal to external dielectric
    epsilon = float(Dext) / float(Dint)
    num_unknowns = mesh.num_node_patches

    rhs, fh_peaks = get_fh_peaks(mesh, epsilon, kappa, Dext, Dint)
    GMRES_Matrix = BEM_GMRES_matrix(mesh, epsilon, kappa)

    #rhs_out = open("2OZO-rhs.txt",'w')
    #for i,(rhsf,rhsh) in enumerate(zip(rhs[:num_unknowns],rhs[num_unknowns:])):
        #print >>rhs_out, "Patch: %d rhs f/h: %f %f" %(i, rhsf, rhsh)
    #rhs_out.close()
    x = array([0.0 for xx in range(len(rhs))])
    residual = max([r*r for r in rhs])
    print "Residual: ", residual
    if residual > 1e-16:

        b = array(rhs)
        print b
        restart_param = min([num_unknowns, GMRES_RESTART])

        if GMRES_PRINT_DEBUG:
            print "starting gmres - %d unknowns" %(num_unknowns)

        x, info = iterative.gmres(GMRES_Matrix,
                                  b,
                                  x0=None,
                                  tol=GMRES_TOLERANCE,
                                  restrt=restart_param,
                                  maxiter=None,
				  xtype='d')

        if info != 0: raise Exception
        if GMRES_PRINT_DEBUG:
            print "done gmres - %d iterations" %(GMRES_Matrix.num_iterations)

    print x[:5], x[num_unknowns:num_unknowns+5]
    locs = []
    fvals = []
    hvals = []
    for (np,f,h,fpeak,hpeak) in zip(mesh.node_patches,x[:num_unknowns],
                                    x[num_unknowns:],
                                    fh_peaks[:num_unknowns],
                                    fh_peaks[num_unknowns:]):
        np.f = float(f) + float(fpeak)
        np.h = float(h) + float(hpeak)

    #mesh.interpolate()

    #for np in mesh.node_patches:
        #qp_list = []
        #for tri_idx in range(np.num_triangles):
            #tri = np.get_triangle(tri_idx)
            #tri.quad_points(qp_list)
        #locs.extend([qp.pt for qp in qp_list])
        #fvals.extend([qp.f for qp in qp_list])
        #hvals.extend([qp.h for qp in qp_list])

    #import _SPHARM as sp
    #origin = Vector(0,0,0)
    #locs = [np.node() for np in mesh.node_patches]
    #fvals = [np.f for np in mesh.node_patches]
    #hvals = [np.h for np in mesh.node_patches]
    #sph_f = sp.least_squares_spherical_harmonic(locs, fvals, origin, 11)
    #sph_h = sp.least_squares_spherical_harmonic(locs, hvals, origin, 11)

    #print "Solved spherical harmonic approximation of f on surface:"
    #force = Vector(0,0,0)
    #flux = 0.0
    #total_area = 0.0

    #rf_force = Vector(0,0,0)

    #from math import pi, cos, sin
    #num_theta_pts = 90
    #num_phi_pts = 2 * num_theta_pts
    #inc_theta = pi / float(num_theta_pts)
    #inc_phi = 2.0*pi / float(num_phi_pts)
    #for theta_ctr in range(num_theta_pts):
        #theta = theta_ctr*inc_theta
        #for phi_ctr in range(num_phi_pts):
            #phi = phi_ctr*inc_phi
            #node = Vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))

            #normal = node.normal()
            #area = sph_f.dArea(node, origin).length()*inc_theta*inc_phi
            #fval = sph_f.evaluate(node, origin)
            #hval = sph_h.evaluate(node, origin)

            #dtheta = sph_f.dtheta(node, origin)
            #dphi = sph_f.dphi(node, origin)

            #flux += -area * hval * epsilon
            #total_area += area
            #Eout = (dtheta + dphi + normal*hval) * -1.0
            #Ein = (dtheta + dphi + normal*hval*epsilon) * -1.0
            #force -= normal * area * (normal.dot(E_to_f(Ein, normal, 0, fval)) - normal.dot(E_to_f(Eout, normal, kappa, fval)))

    #print "Total area:" , total_area
    #print "Net force: ", force
    #print "Net flux: ", flux


    #for np in mesh.node_patches:
        #for tri_idx in range(np.num_triangles):

            #qp_list = []
            #tri = np.get_triangle(tri_idx)

            #fval, hval, dtheta, dphi = 0.0, 0.0, 0.0, 0.0
            #node = tri.centre
            #normal = tri.normal
            #area = tri.area
            #fval = sph_f.evaluate(node, origin)
            #hval = sph_h.evaluate(node, origin)
            #flux += hval*area*epsilon
            #dtheta = sph_f.dtheta(node, origin)
            #dphi = sph_f.dphi(node, origin)
            #print dtheta, dphi
            #E =  (dtheta + dphi + np.normal()*hval) * -1.0
            #E2 = E.dot(E)
            #Ex = E.x
            #Ey = E.y
            #Ez = E.z

            ##print np.f, fval, np.h, hval, node, normal, E, E.length(), area, E2

            ## this is the MST premultiplier, in internal units (Angstroms,
            ## relative permittivities etc.)
            #mult = normal * area

            ## x stress
            #stress_xx = Ex * Ex - 0.5*E2 - 0.5*kappa*kappa*fval*fval
            #stress_xy = Ex * Ey
            #stress_xz = Ex * Ez

            #fx = mult.dot(Vector(stress_xx, stress_xy, stress_xz))

            ## y stress
            #stress_yx = Ey * Ex
            #stress_yy = Ey * Ey - 0.5*E2 - 0.5*kappa*kappa*fval*fval
            #stress_yz = Ey * Ez

            #fy = mult.dot(Vector(stress_yx, stress_yy, stress_yz))

            ## z stress
            #stress_zx = Ez * Ex
            #stress_zy = Ez * Ey
            #stress_zz = Ez * Ez - 0.5*E2 - 0.5*kappa*kappa*fval*fval

            #fz = mult.dot(Vector(stress_zx, stress_zy, stress_zz))

            #force += Vector(fx, fy, fz)

    #print "Net force: ", force
    #print "Net flux: ", flux
    return

def hack_solveElectrostatics(mesh, fh_results):
    """Assign BEM results from fh file."""

    flines = open(fh_results,'r').readlines()
    for np, line in zip(mesh.node_patches, flines):
        bits = line.split()
        np.f = float(bits[0])
        np.h = float(bits[1])
    return

def calc_potential(mesh, pt, kappa, Dext, Dint):

    # unpack the charges and charge positions
    charges = [c.charge for c in mesh.charges]
    charge_positions = [list(c.position) for c in mesh.charges]

    pot = mesh.calc_potential_at_point(pt,
                                       charges,
                                       charge_positions,
                                       kappa,
                                       Dint,
                                       Dext)
    return pot

def solve_energy(mesh, kappa, Dext, Dint):

    # unpack the charges and charge positions
    charges = [c.charge for c in mesh.charges]
    charge_positions = [list(c.position) for c in mesh.charges]

    # this gets the reaction field energy contribution
    E = mesh.calculate_energy(charges,
                              charge_positions,
                              kappa,
                              Dint,
                              Dext)

    # convert to kJ
    import constants
    con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
      / (constants.epsilon0 * constants.Angstroms * 1000.0)
    #print E,con,E*con

    return E*con

if __name__=="__main__":

    #import constants
    #print constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
      #/ (constants.epsilon0 * constants.Angstroms * 1000.0)

    import sys

    Dext = 80.0
    Dint = 2.0
    kappa = 0.0
    radius = 3.0
    charge = +1.0
    run_model = False
    fh_file = "results.fh"
    node_patch_kinemage = "node_patch_mesh.kin"
    print "creating mesh ...",
    if (len(sys.argv) == 2):
        from config_file import *
        cfg = ConfigFile(sys.argv[1])
        mesh = StaticEnsemble(sys.argv[1])
        Dext = cfg.Dext
        Dint = cfg.Dint
        kappa = cfg.kappa
        
        base_filename = sys.argv[1]
        idx = base_filename.find(".")
        base_filename = base_filename[:idx]
        fh_file = cfg.fh_file = "%s.fh" %(base_filename)
        node_patch_kinemage = "%s.kin" %(base_filename)
        print "Base filename: ", base_filename
    elif (len(sys.argv) >= 3):

        exe, gts_file, pqr_file = sys.argv[:3]
        mesh = GTS_Mesh(gts_file,pqr_file)

        if (len(sys.argv) == 5):
            run_model = True
            mesh_fh = sys.argv[3]
            ellipse_model = sys.argv[4]

    else:
        mesh = Born_Ion_Mesh(Dext, Dint, num_subdivides=4, radius=radius, charge=charge)

    print "done"

    #BEM_matrix_pickle_file = os.path.splitext(sims_file)[0]+"_BEM_full_matrix.pickle"

    #if os.path.exists(BEM_matrix_pickle_file):
        #print "%s exists! loading full BEM matrix ... " %(BEM_matrix_pickle_file),
        #fin = open(BEM_matrix_pickle_file,'rb')
        #full_matrix = cPickle.load(fin)
        #print "done"

    #else:
        #print "%s doesn't exist - calculating everything from scratch" %(BEM_matrix_pickle_file)

        #print "calculating full BEM matrix ... ",
        #full_matrix = BEM_Full_Matrix(mesh, kappa, Dext, Dint)
        #print "done"

        #print "pickling ... ",
        #fout = open(BEM_matrix_pickle_file, 'wb')
        #cPickle.dump(full_matrix, fout, protocol=2)
        #print "done"

    #print "validating quadrature rules"
    #mesh.validate_quadrature()
    #print "done"
    #del full_matrix.matrix_inverse
    #print linalg.cond(full_matrix.matrix)
    #origin = _BEM.Vector(0.0,0.0,0.0)
    #sa = sum( [v.solid_angle(origin, v.normal) for v in mesh] )
    #print sa, pi*4.0

    # analytic solution

    mesh.kinemage(node_patch_kinemage)
    #solveElectrostatics(mesh, fh_file)
    if run_model:
        fh_vals = [[float(xx) for xx in line.split()]
                   for line in open(mesh_fh,'r').readlines()]
        fvals = [x[0] for x in fh_vals]
        hvals = [x[1] for x in fh_vals]
        mesh.set_fh(fvals, hvals)

        E = solve_energy(mesh, kappa, Dext, Dint)
        print "Solvation energy: ", E

        vertices = gts_utils.get_vertices_from_gts(ellipse_model)
        for v in vertices:
            pot = calc_potential(mesh, v, kappa, Dext, Dint)
            print "Potential at %s is %f" %(v, pot)

    else:
        from constants import epsilon0
        solveElectrostatics(mesh, kappa, Dext, Dint)
        mesh.write_fh(fh_file)

        #to_kcal = 96485.3383 * 1e-3 * 0.2388458966275 * 1e10 * 1.6e-19 / constants.epsilon0
        #for np in mesh.node_patches:
        #    print np.f*to_kcal, np.h*to_kcal

        mesh.calc_energies(Dext, Dint, kappa)
        try:
            
            mesh.calc_fluxes(Dext/Dint)
            rf_force = mesh.calc_reaction_field_forces(Dext, Dint, kappa)
            boundary_force = mesh.boundary_force(Dext, Dint, kappa)
            print "Sum of forces: ", rf_force + boundary_force
        except AttributeError:
            pass

        if (type(mesh) == Born_Ion_Mesh):
            print "born ion analytic solvation energy: ", mesh.analytic_solvation_energy
        #print mesh.area, 4.0*pi*radius*radius

        #_BEM.get_max_recursions()
        #solveFullElectrostatics(mesh, full_matrix)
