import _SIMS_BEM
from vector import *
from constants import *
import pqrtools
from numpy import zeros, array, ones, matrix, linalg
from scipy.linalg import iterative
import time
import os
import cPickle
import geometry

SIMS_PRINT_DEBUG = True
GMRES_TOLERANCE = 1.0e-3
GMRES_RESTART = 5

#class BEM_Full_Matrix(object):
    
    #def __init__(self, mesh, kappa, Dext, Dint):

        ##self.mesh = mesh
        #self.kappa = kappa
        #self.Dext = Dext
        #self.Dint = Dint

        ## create BEM matrix
        #self.matrix = self.init_BEM_matrix(mesh, self.kappa, self.Dext, self.Dint)
        #self.matrix_inverse = self.matrix.I
    
    #@staticmethod
    #def init_BEM_matrix(mesh, kappa, Dext, Dint):

        #N = len(mesh.vertices)
        #bem_matrix = matrix(zeros((N*2,N*2),'f'))
        
        #epsilon = float(Dext)/float(Dint)
        #res = [0.0, 0.0, 0.0, 0.0]
        #for i, row_vert in enumerate(mesh.vertices):
            #for j, col_vert in enumerate(mesh.vertices):
                #_SIMS_BEM.calculateMatrixElements(res,epsilon,kappa,
                                                  #row_vert,col_vert)
                #A,B,C,D = res
                #bem_matrix[i,j] = B
                #bem_matrix[i,j+N] = A
                #bem_matrix[i+N,j] = D
                #bem_matrix[i+N,j+N] = C
        #return bem_matrix

class CPP_Mesh(object):
    
    def __init__(self):
        self.__mesh = None
        self.__charges = None

    # Charges property
    def __get_charges(self):
        return self.__charges
    def __set_charges(self, new_charges):
        self.__charges = new_charges
        return
    charges = property(__get_charges, __set_charges, doc="Get/set Charge objects within Mesh.")

    # Mesh property
    def __get_mesh(self):
        return self.__mesh
    def __set_mesh(self, new_mesh):
        self.__mesh = new_mesh
        return
    mesh = property(__get_mesh, __set_mesh, doc="Get/set the C++ BEM Mesh object.")
    
    @property
    def num_vertices(self):
        return self.mesh.num_points

    def __getitem__(self, idx):
        return self.mesh.get(idx)

    def __len__(self):
        return self.num_vertices

class SIMS_Mesh(CPP_Mesh):

    def __init__(self, sims_filename, pqr_filename):
        
        # call base class constructor
        super(CPP_Mesh, self).__init__()
        
        # set mesh
        self.mesh = self.load_sims_mesh(sims_filename)
        
        # set charges
        self.charges = pqrtools.getCharges(pqr_filename)
        
        return
    
    @staticmethod
    def load_sims_mesh(filename):

        f = open(filename, "r")
        all_lines = [line for line in f.readlines() if line.strip()[0] != "#"]

        num_vertices = len(all_lines)
        sims_mesh = _SIMS_BEM.Mesh(num_vertices)

        for i,line in enumerate(all_lines):

            # decode the line of the SIMS format file
            area, x,y,z, nx,ny,nz = [float(val) for val in line.split()[-7:]]

            # use C++ function to set the contents of the struct
            sims_mesh.set_quad_point(i, area, x,y,z, nx,ny,nz, 0.5);
            
        return sims_mesh

class Spherical_Mesh(CPP_Mesh):

    def __init__(self, radius=1.0, num_subdivides=3, xyz=Vector(0.0,0.0,0.0), rot=None):
        """Create a spherical mesh based on a subdivided icosahedron."""
                
        # call base class constructor
        super(Spherical_Mesh, self).__init__()
        
        import _GTSMESH as gts
        import tempfile
        spherical_file = tempfile.mktemp(suffix=".gts", 
                                         prefix="gts_sphere_", 
                                         dir=".")
        gts.write_spherical_mesh(num_subdivides, spherical_file)

        import gts_utils
        vertices, triangles = gts_utils.get_vertices_triangles_from_gts(spherical_file)
        
        num_vertices = len(triangles)        
        self.mesh = _SIMS_BEM.Mesh(num_vertices)
        # assign quadrature points at centres of triangles
        for i,t in enumerate(triangles):

            v1 = vertices[t[0]]
            v2 = vertices[t[1]]
            v3 = vertices[t[2]]
            
            # apply rotation if passed in
            if rot is not None:
                #rot = geometry._rand_rot()
                v1 = geometry.apply_quaternion_to_vector(rot, v1).normal()
                v2 = geometry.apply_quaternion_to_vector(rot, v2).normal()
                v3 = geometry.apply_quaternion_to_vector(rot, v3).normal()

            # apply xyz offset
            v1 = v1*radius + xyz
            v2 = v2*radius + xyz
            v3 = v3*radius + xyz
            
            tt = geometry.Triangle(v1, v2, v3)
            
            x,y,z = tt.centre
            nx,ny,nz = tt.normal
            
            # sanity check
            assert((tt.centre - xyz).normal() * tt.normal.normal() - 1.0 < 1e-9)
                
            self.mesh.set_quad_point(i, tt.area, x,y,z, nx,ny,nz, 0.5);
        
        return

class Born_Ion_Mesh(Spherical_Mesh):

    def __init__(self, radius=1.0, num_subdivides=3, xyz=Vector(0.0,0.0,0.0), charge=1.0):
        """Create a Born Ion mesh. Basically a SphericalMesh with a charge."""

        # randomize the orientation
        rot = geometry._rand_rot()
        
        # call Spherical_Mesh constructor
        super(Born_Ion_Mesh, self).__init__(radius, num_subdivides, xyz, rot)
        
        # Add a charge
        from entities import Charge
        self.charges = [Charge(xyz, charge, radius=radius, name="Born Ion")]
        
        return
    
def get_fh_peaks(mesh, epsilon, kappa, Dext):
    """Get the fh peaks vector (fast varying component approximation)."""
    
    # unpack the charges and charge positions
    charges = [c.charge for c in mesh.charges]
    charge_positions = [list(c.position) for c in mesh.charges]
    
    num_unknowns = mesh.num_vertices
    
    fh_peaks = [0.0 for xx in range(2*num_unknowns)]
    rhs = [0.0 for xx in range(2*num_unknowns)]
    _SIMS_BEM.calculate_fh_peaks(fh_peaks, 
                                 rhs, 
                                 charges, 
                                 charge_positions, 
                                 mesh.mesh, 
                                 epsilon, 
                                 kappa, 
                                 Dext)

    return array(rhs), array(fh_peaks)

class BEM_GMRES_matrix(object):
    
    def __init__(self, 
                 SIMS_Mesh, 
                 epsilon,
                 kappa):

        self.SIMS_Mesh = SIMS_Mesh
        self.kappa = kappa
        self.epsilon = epsilon
        self.N = SIMS_Mesh.num_points
        #print "DEBUG BEM_GMRES_MATRIX: N=", self.N
        self.num_iterations = 0
        
        return

    def matvec(self, x):
        """Matrix-vector multiplication -- required function for GMRES."""

        self.num_iterations += 1
        # TODO: debug mode only -- use logging module
        if SIMS_PRINT_DEBUG:
            print "iteration: %d (%d)" %(self.num_iterations, self.N)
        
        # sanity check the incoming vector
        assert(len(x) == 2*self.N)

        # convert numpy array into a list
        # for C++ CUDA call (TODO: clean up interface)
        xx = [val for val in x]
        
        result = self._matrix_vector_multiply(xx)
        
        # all done; return the result
        return result
    
    def _matrix_vector_multiply(self, xx):
        """Call C++ implementation of matrix-vector multiplication."""

        start = time.clock()
        res = [0.0 for i in range(2*self.N)]
        xxx = [float(x) for x in xx]
        _SIMS_BEM.multiplyVector(xxx, 
                                 res, 
                                 self.SIMS_Mesh, 
                                 self.epsilon,
                                 self.kappa)
    
        if SIMS_PRINT_DEBUG:
            print "MultiplyVector: %f" %(time.clock() - start)
        
        return array(res)

def solveFullElectrostatics(mesh, full_matrix):
    
    kappa = full_matrix.kappa
    Dext = full_matrix.Dext
    Dint = full_matrix.Dint
    epsilon = Dext / Dint
    num_unknowns = mesh.num_vertices
    
    inverse = full_matrix.matrix_inverse
    
    rhs, fh_peaks = get_fh_peaks(mesh, epsilon, kappa, Dext)
    
    solution = inverse * matrix(rhs).T
    
def solveElectrostatics(mesh, kappa, Dext, Dint):
    """Use BEM to solve electrostatic parameters f and h on surface of all
    diffusing entities."""
    
    # epsilon is ratio of internal to external dielectric
    epsilon = float(Dext) / float(Dint)
    num_unknowns = mesh.num_vertices

    rhs, fh_peaks = get_fh_peaks(mesh, epsilon, kappa, Dext)
    GMRES_Matrix = BEM_GMRES_matrix(mesh.mesh, epsilon, kappa)
    
    rhs_offset = GMRES_Matrix.matvec(fh_peaks)
    rhs -= rhs_offset

    if sum([r*r for r in rhs]) != 0.0:
    
        b = array(rhs)
        restart_param = min([num_unknowns, GMRES_RESTART])

        if SIMS_PRINT_DEBUG:
            print "starting gmres - %d unknowns" %(num_unknowns)

        x, info = iterative.gmres(GMRES_Matrix,
                                  b,
                                  x0=None,
                                  tol=GMRES_TOLERANCE,
                                  restrt=restart_param,
                                  maxiter=None,
                                  xtype='f')
        
        if info != 0: raise Exception
        if SIMS_PRINT_DEBUG:
            print "done gmres - %d iterations" %(GMRES_Matrix.num_iterations)

        # verify that the matrix inversion does appear to have worked!
        #y = GMRES_Matrix.matvec(x)
        #print y, rhs
        #print max([abs(xx) for xx in (y - rhs)])
                
        x += fh_peaks # add in slowly varying part of the peak-splitting solution
        
    else:
        x = fh_peaks
        
        if SIMS_PRINT_DEBUG:
            print "fh peak splitting gave exact solution"
            if epsilon != 1.0: print "HIGHLY SUSPICIOUS!!"
    
    for i,(f,h) in enumerate(zip(x[:num_unknowns],x[num_unknowns:])):
        vert = mesh[i]
        vert.f = float(f)
        vert.h = float(h)

    return

def solve_energy(mesh, kappa, Dext, Dint):
    
    # unpack the charges and charge positions
    charges = [c.charge for c in mesh.charges]
    charge_positions = [list(c.position) for c in mesh.charges]
   
    # this gets the reaction field energy contribution
    E = _SIMS_BEM.calculate_energy(mesh.mesh, 
                                   charges, 
                                   charge_positions, 
                                   kappa, 
                                   Dint, 
                                   Dext)
    
    # add 'self' energy term
    #E += de.mesh.self_energy(kappa, Dint)

    # convert to kJ
    import constants
    E *= constants.elementary_charge*constants.elementary_charge \
      / (constants.epsilon0 * constants.__Angstroms * 1000.0)
    
    return E

if __name__=="__main__":
    
    import sys
    exe, sims_file, pqr_file = sys.argv

    Dext = 80.0
    Dint = 2.0
    kappa = 0.0

    print "creating mesh ...",
    mesh = SIMS_Mesh(sims_file,pqr_file)
    #mesh = Born_Ion_Mesh(num_subdivides=4, radius=3.0)
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

    #del full_matrix.matrix_inverse
    #print linalg.cond(full_matrix.matrix)
    solveElectrostatics(mesh, kappa, Dext, Dint)
    E = solve_energy(mesh, kappa, Dext, Dint)
    print E * Avogadro
    #solveFullElectrostatics(mesh, full_matrix)
    
 
