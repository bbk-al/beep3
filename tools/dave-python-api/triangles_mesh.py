import _BEM
from _BEM import Vector
import constants
import pqrtools
from numpy import zeros, array, ones, matrix, linalg
from scipy.linalg import iterative
import time
import os
import cPickle
import geometry
import gts_utils
from math import pi

GMRES_PRINT_DEBUG = True
GMRES_TOLERANCE = 1.0e-6
GMRES_RESTART = 50

class CPP_Mesh(_BEM.Mesh, object):
    
    def __init__(self):
        super(CPP_Mesh, self).__init__()

    # generator function to provide triangle iterator
    @property
    def triangles(self):
        for x in xrange(self.num_triangles):
            yield self.get_triangle(x)
    
    # generator function to provide node patch iterator
    @property    
    def node_patches(self):
        for x in xrange(self.num_node_patches):
            yield self.get_node_patch(x)

    # generator function to provide charge iterator
    @property
    def charges(self):
        for x in xrange(self.num_charges):
            yield self.get_charge(x)
            
    def kinemage(self, filename="mesh.kin"):
        
        output = ["@kinemage"]
        from kintools import randomColour
        
        output.append("@trianglelist {triangles}")    
        for t in self.triangles:
            a = (ax, ay, az) = t.v1
            b = (bx, by, bz) = t.v2
            c = (cx, cy, cz) = t.v3
            colour = randomColour()
            output.append("X %f %f %f %f %f %f %s %f %f %f" %(ax, ay, az,
                                                   bx, by, bz, colour,
                                                   cx, cy, cz ))
            
        output.append("@vectorlist {centres} off")    
        for t in self.triangles:
            (x1, y1, z1) = t.centre
            (x2, y2, z2) = t.centre + t.normal
            output.append("{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2))

        #output.append("@vectorlist {quad_points} off")    
        #for t in self.triangles:
            #qp_list = []
            #t.quad_points(qp_list)
            #for qp in qp_list:
                #(x1, y1, z1) = qp.pt
                #(x2, y2, z2) = qp.pt + t.normal
                #output.append("{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2))

        output.append("@vectorlist {vertex_normals}")
        for np in self.node_patches:
            v = np.node()
            vn = np.normal()
            (x1, y1, z1) = v
            (x2, y2, z2) = v + vn
            output.append("{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2))

        output.append("@trianglelist {node_patches} off")
        for np in self.node_patches:
            colour = randomColour()
            for ii in xrange(np.num_triangles):
                
                t = np.get_triangle(ii)
                
                a = (ax, ay, az) = t.v1
                b = (bx, by, bz) = t.v2
                c = (cx, cy, cz) = t.v3
                output.append("X %f %f %f %f %f %f %s %f %f %f" %(ax, ay, az,
                                                                  bx, by, bz, colour,
                                                                  cx, cy, cz ))

        # don't bother doing this for very large meshes
        if self.num_node_patches < 1000:
            output.append("@vectorlist {quad_points} off")    
            for np in self.node_patches:
                for ii in xrange(np.num_triangles):
                    
                    t = np.get_triangle(ii)
                    qp_list = []
                    t.quad_points(qp_list)
                    for qp in qp_list:
                        (x1, y1, z1) = qp.pt
                        (x2, y2, z2) = qp.pt + t.normal
                        output.append("{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2))

        output.append("@spherelist {charges} color=purple off")
        for ch in self.charges:
            x,y,z = ch.position
            output.append("r=%f {} %f %f %f" %(ch.radius,x,y,z))
                        
        f = open(filename, 'w')
        print >>f, "\n".join(output)
        f.close()
        
        return
    
    def calc_fluxes(self, epsilon):
        """Calculate the flux E.dot(n) over the surface of a meshed entity."""
        try:
    
            for k in self.entity_mappings.keys():
                [(first, num), (first_tri, num_tris), (first_ch, num_ch) ] = self.entity_mappings[k]
                flux = self.calculate_net_flux(first, num, epsilon)
                print "Flux on mesh %d = %f" %(k, flux)
        
        except AttributeError:
            # for meshes with no entity mapping of object to node patch range
            flux = self.calculate_net_flux(0, self.num_node_patches, epsilon)
            print "Flux = %f" %(flux)
        return
    
    def calc_energies(self, Dext, Dint, kappa):
        
        import constants
        con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
          / (constants.epsilon0 * constants.Angstroms * 1000.0)
        
        try:
    
            for k in self.entity_mappings.keys():
                [(first, num), (first_tri, num_tris), (first_ch, num_ch) ] = self.entity_mappings[k]
                energy = self.calculate_energy(first, 
                                               num,
                                               first_ch,
                                               num_ch,
                                               kappa,
                                               Dint, 
                                               Dext) * con
                print "Solvation energy for mesh %d = %f" %(k, energy)
        
        except AttributeError:
            # for meshes with no entity mapping of object to node patch range
            energy = self.calculate_energy(0, self.num_node_patches, epsilon) * con
            print "Solvation energy = %f" %(energy)
        
        return
    
    def net_force(self, kappa):
        """Calculate net force over all triangular elements in the mesh."""

        try:
    
            for k in self.entity_mappings.keys():
                [(first, num), (first_tri, num_tris), (first_ch, num_ch) ] = self.entity_mappings[k]
                net_force = self.force_on_triangles(first_tri, num_tris, kappa)
                print "Net force on mesh %d = %s" %(k, net_force)
        
        except AttributeError:
            pass
        
        net_force = self.force_on_triangles(0, self.num_triangles, kappa)
        print "Total net force = %s" %(net_force)
            
        return    

class StaticEnsemble(CPP_Mesh):
    
    def __init__(self, config_filename):

        # call base class constructor
        super(StaticEnsemble, self).__init__()
                
        from config_file import *
        cfg = ConfigFile(config_filename)
        
        self.entity_mappings = {}
        
        from gts_utils import get_vertices_triangles_from_gts
        
        for ctr, sim_def in enumerate(cfg.sim_defs):

            gts_filename = sim_def.mesh_def.gts_file
            centre_filename = sim_def.mesh_def.centre_file
            vertices, triangles = get_vertices_triangles_from_gts(gts_filename)
            vertex_numbering_offset = self.num_vertices

            self.entity_mappings[ctr] = [ (vertex_numbering_offset, len(vertices)) ]
            
            x,y,z = [atof(xx) 
                     for xx in open(centre_filename,'r').readline().split()]
            offset = sim_def.xyz_offset - Vector(x,y,z)
            
            for v in vertices:
                self.add_vertex(v + offset)
            
            self.entity_mappings[ctr].append( (self.num_triangles, len(triangles)) )
            for t in triangles:
                self.define_triangle(t[0] + vertex_numbering_offset,
                                     t[1] + vertex_numbering_offset,
                                     t[2] + vertex_numbering_offset)
            
            charge_file = open(sim_def.mesh_def.xyzq_file,'r').readlines()
            charge_list = []
            for line in charge_file:
                x,y,z,q,r = [atof(xx) for xx in line.split()]
                new_xyz = Vector(x,y,z) + offset
                charge_list.append(_BEM.Charge(q, new_xyz, r))
            self.entity_mappings[ctr].append( (self.num_charges, len(charge_list)) )
            self.add_charges(charge_list)
            
        # done adding all vertices, triangles, charges for all meshed bits
        self.create_node_patches()
        
        return
     
class GTS_Mesh(CPP_Mesh):

    def __init__(self, gts_filename, pqr_filename):
        """Create a mesh from GTS format file."""
                
        # call base class constructor
        super(GTS_Mesh, self).__init__()
        
        import gts_utils
        vertices, triangles = gts_utils.get_vertices_triangles_from_gts(gts_filename)

        # add vertices to the mesh (as Vertex objects)
        for v in vertices:
            self.add_vertex(v)

        # define the triangles in the mesh
        num_triangles = len(triangles)
        for i,t in enumerate(triangles):
            self.define_triangle(t[0], t[1], t[2]);
            
        # init NodePatch objects
        self.create_node_patches()
            
        # add charges
        cpp_charges = [_BEM.Charge(c.charge,c.position) 
                       for c in pqrtools.getCharges(pqr_filename)]
        self.add_charges(cpp_charges)
        
        return

class Spherical_Mesh(CPP_Mesh):

    def __init__(self, radius=1.0, num_subdivides=3, xyz=Vector(0.0,0.0,0.0), rot=None):
        """Create a spherical mesh based on a subdivided icosahedron."""
                
        # call base class constructor
        super(Spherical_Mesh, self).__init__()

        # use gts_utils to get a spherical mesh
        vertices, triangles = gts_utils.get_spherical_mesh(num_subdivides)

        # rotate vertices
        if rot is not None:
            vertices = [geometry.apply_quaternion_to_vector(rot, v).normal()
                        for v in vertices]

        # store vertex normals
        vnormals = vertices[:]
            
        # apply radial extension and xyz offset
        vertices = [v*radius + xyz for v in vertices]

        # add vertices to the mesh (as Vertex objects)
        for v in vertices:
            self.add_vertex(v)

        # define the triangles in the mesh
        num_triangles = len(triangles)
        for i,t in enumerate(triangles):

            self.define_triangle(t[0], t[1], t[2]);        

        # init NodePatch objects
        self.create_node_patches()
            
        self.vnormals = [(v,vn) for v,vn in zip(vertices, vnormals)]
        
        return

class Born_Ion_Mesh(Spherical_Mesh):

    def __init__(self, Dext, Dint, radius=1.0, num_subdivides=3, xyz=Vector(0.0,0.0,0.0), charge=1.0):
        """Create a Born Ion mesh. Basically a SphericalMesh with a charge."""

        # randomize the orientation
        rot = geometry._rand_rot()
        
        # call Spherical_Mesh constructor
        super(Born_Ion_Mesh, self).__init__(radius, num_subdivides, xyz, rot)
        
        # Add a charge
        self.add_charges([_BEM.Charge(charge, xyz, radius)])

        import constants
        premult = -constants.Avogadro * (constants.elementary_charge ** 2) \
            / (1000.0*8.0*constants.pi*constants.epsilon0*constants.Angstroms)

        self.__analytic = premult * (1.0/Dint - 1.0/Dext) * (charge*charge) /radius
        
        return
    
    @property
    def analytic_solvation_energy(self):
        """The anayltical Born solvation energy for this cavity."""
        return self.__analytic
    
def get_fh_peaks(mesh, epsilon, kappa, Dext, Dint):
    """Get the fh peaks vector (fast varying component approximation)."""
    
    mesh.set_rhs(kappa, Dext, Dint)

    fh_peaks = [np.f for np in mesh.node_patches]
    h_peaks = [np.h for np in mesh.node_patches]
    fh_peaks.extend(h_peaks)
    
    rhsf = [np.rhs_f for np in mesh.node_patches]
    rhsh = [np.rhs_h for np in mesh.node_patches]
    rhsf.extend(rhsh)
    #print rhsf
    return array(rhsf), array(fh_peaks)

class BEM_GMRES_matrix(object):
    
    def __init__(self, 
                 mesh, 
                 epsilon,
                 kappa):

        self.mesh = mesh
        self.kappa = kappa
        self.epsilon = epsilon
        self.N = mesh.num_node_patches
        #print "DEBUG BEM_GMRES_MATRIX: N=", self.N
        self.num_iterations = 0
        
        return

    def matvec(self, x):
        """Matrix-vector multiplication -- required function for GMRES."""

        self.num_iterations += 1
        # TODO: debug mode only -- use logging module
        if GMRES_PRINT_DEBUG:
            print "iteration: %d (%d)" %(self.num_iterations, self.N)
        
        # sanity check the incoming vector
        assert(len(x) == 2*self.N)

        result = self._matrix_vector_multiply(x)
        
        # all done; return the result
        return result
    
    def _matrix_vector_multiply(self, xx):
        """Call C++ implementation of matrix-vector multiplication."""

        start = time.clock()
        res = [0.0 for i in range(self.N*2)]
        
        for valf, valh ,np in zip(xx[:self.N], xx[self.N:], self.mesh.node_patches):
            np.f = float(valf)
            np.h = float(valh)
            
        mesh.integrate_node_patches(self.epsilon, self.kappa)
        
        if GMRES_PRINT_DEBUG:
            print "MultiplyVector: %f" %(time.clock() - start)

        for i,np in enumerate(self.mesh.node_patches):
            res[i] = np.f_accum
            res[i+self.N] = np.h_accum
            
        return array(res)

def solveFullElectrostatics(mesh, full_matrix):
    
    kappa = full_matrix.kappa
    Dext = full_matrix.Dext
    Dint = full_matrix.Dint
    epsilon = Dext / Dint
    num_unknowns = mesh.num_node_patches
    
    inverse = full_matrix.matrix_inverse
    
    rhs, fh_peaks = get_fh_peaks(mesh, epsilon, kappa, Dext, Dint)
    
    solution = inverse * matrix(rhs).T
    
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
        
    if sum([r*r for r in rhs]) != 0.0:

        b = array(rhs)
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
    for (np,f,h,fpeak,hpeak) in zip(mesh.node_patches,x[:num_unknowns],x[num_unknowns:], fh_peaks[:num_unknowns], fh_peaks[num_unknowns:]):
        np.f = float(f) + float(fpeak)
        np.h = float(h) + float(hpeak)

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
    print "creating mesh ...",
    if (len(sys.argv) == 2):
        from config_file import *
        cfg = ConfigFile(sys.argv[1])
        mesh = StaticEnsemble(sys.argv[1])
        Dext = cfg.Dext
        Dint = cfg.Dint
        kappa = cfg.kappa
        
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
    
    mesh.kinemage("node_patch_mesh.kin")
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
        solveElectrostatics(mesh, kappa, Dext, Dint)
        mesh.write_fh("results.fh")
    
        mesh.calc_energies(Dext, Dint, kappa)
        mesh.calc_fluxes(Dext/Dint)
        mesh.net_force(kappa)
        
        if (type(mesh) == Born_Ion_Mesh):
            print "born ion analytic solvation energy: ", mesh.analytic_solvation_energy
        #print mesh.area, 4.0*pi*radius*radius
        
        #_BEM.get_max_recursions()
        #solveFullElectrostatics(mesh, full_matrix)
