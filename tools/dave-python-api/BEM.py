from vector import Vector, Rotation, Quaternion
import _BEM
import math
from numpy import zeros, array, ones
from scipy.linalg import iterative
import time

# Global parameters for GMRES
GMRES_RESTART = 500
GMRES_TOLERANCE = 1e-6

# global parameter to control use of geometric correction factors
BEM_USE_GEOMETRIC_CORRECTIONS = False # no improvement in accuracy at the moment (4/5/09)

# debug messages
BEM_PRINT_DEBUG = True

def VectorPoint(scientific_vector_obj):
    """Convert Scientific.Geometry.Vector object into _BEM.VectorPoint struct"""

    return _BEM.VectorPoint(scientific_vector_obj.x(),
                            scientific_vector_obj.y(),
                            scientific_vector_obj.z())

def get_charges_and_positions(diffusing_entities):
    """Return lists of charges and charge positions."""
    
    charges = []
    charge_positions = []

    for de in diffusing_entities:
        c,p = de.universe_charges_and_positions()
        charges.extend(c)
        charge_positions.extend(p)
    
    return charges, charge_positions

def get_RHS(diffusing_entities, c_patches, Dext, kappa):
    """Return the RHS for the matrix-vector equation (as Numpy array).

    Underlying call is to _BEM.so C++ library object."""
    
    # get the charges and charge positions
    charges, charge_positions = get_charges_and_positions(diffusing_entities)
    
    num_unknowns = len(c_patches)
    rhs = [0.0 for xx in range(2*num_unknowns)]
    _BEM.calculateRHS(rhs, charges, charge_positions, c_patches, Dext, kappa)
    return array(rhs)

def get_fh_peaks(diffusing_entities, c_patches, epsilon, kappa, Dext):
    """Get the fh peaks vector (fast varying component approximation)."""
    
    # get the charges and charge positions
    charges, charge_positions = get_charges_and_positions(diffusing_entities)
    
    num_unknowns = len(c_patches)
    fh_peaks = [0.0 for xx in range(2*num_unknowns)]
    rhs = [0.0 for xx in range(2*num_unknowns)]
    _BEM.calculate_fh_peaks(fh_peaks, rhs, charges, charge_positions, c_patches, epsilon, kappa, Dext)

    return array(rhs), array(fh_peaks)

class BEM_matrix(object):

    """Base class for BEM matrix -- do not use this; use a subclass with
    appropriate implementation of the _matrix_vector_multiply function."""
    
    def __init__(self, 
                 patches, 
                 kappa, 
                 epsilon, 
                 rhs):

        self.patches = patches
        self.kappa = kappa
        self.epsilon = epsilon
        self.N = len(patches)
        self.num_iterations = 0
        self.rhs = rhs[:] # useful for debug
        
        return

    def matvec(self, x):
        """Matrix-vector multiplication -- required function for GMRES."""

        self.num_iterations += 1
        # TODO: debug mode only -- use logging module
        if BEM_PRINT_DEBUG:
            print "iteration: %d (%d)" %(self.num_iterations, self.N)
        
        # sanity check the incoming vector
        assert(len(x) == 2*self.N)

        # convert numpy array into a list
        # for C++ CUDA call (TODO: clean up interface)
        xx = [val for val in x]
        
        # underlying call to possibly specialist (CUDA) implementation not
        # this is not defined in the base class -- the subclasses implement
        # this as they see fit (e.g. by calls to CUDA library, or FMM, or
        # maybe even CUDA-FMM)
        result = self._matrix_vector_multiply(xx)
        
        # the rhs vector is normalised to ones
        #for i, scaling in enumerate(self.rhs):
        #    result[i] /= scaling
        
        # all done; return the result
        return result

class FMM_BEM_matrix(BEM_matrix):
    
    def __init__(self, fmm, *args):
        
        super(FMM_BEM_matrix, self).__init__(*args)

        # create BEM/FMM tree to do the real work
        self.fmm_tree = fmm
        
        return

    def _matrix_vector_multiply(self, xx):
        """Call FMM implementation of matrix-vector multiplication."""

        # Non-CUDA version
        fmm_start = time.clock()
        results, fmm_time, local_time =  self.fmm_tree.matvec(xx,
                                                              self.patches,
                                                              self.kappa, 
                                                              self.epsilon)
        if BEM_PRINT_DEBUG:
            print """FMM: %f (fmm=%f, local=%f)""" \
                  %(time.clock()-fmm_start,fmm_time,local_time)

        return results
    
class CUDA_BEM_matrix(BEM_matrix):
    
    def __init__(self, thread_block_dims, *args):
        
        self._thread_block_dims = thread_block_dims
        super(CUDA_BEM_matrix, self).__init__(*args)
        
        return
    
    def _matrix_vector_multiply(self, xx):
        """Call CUDA implementation of matrix-vector multiplication."""

        block_width, block_height = self._thread_block_dims
        
        # CUDA version
        start_cuda = time.clock()
        res = [0.0 for i in range(2*self.N)]
        _BEM.CUDA_multiplyVector(xx, 
                                 res, 
                                 self.patches, 
                                 self.kappa, 
                                 self.epsilon,
                                 block_width,
                                 block_height)
        #_BEM.multiplyVector(xx, 
        #                    res, 
        #                    self.patches, 
        #                    self.kappa, 
        #                    self.epsilon)
    
        
        if BEM_PRINT_DEBUG:
            print "CUDA (%dx%d): %f" %(block_width, 
                                       block_height, 
                                       time.clock() - start_cuda)
        return array(res)
        
    def set_thread_block_dimensions(self, width, height):
        """Set the CUDA thread block dimensions."""
        
        self._thread_block_dims = (width, height)
        return
    
def solveElectrostatics(diffusing_entities, kappa, Dext, Dint):
    """Use BEM to solve electrostatic parameters f and h on surface of all
    diffusing entities."""
    
    # epsilon is ratio of internal to external dielectric
    epsilon = Dext / Dint

    patches = get_patches_from_entities(diffusing_entities)
    num_unknowns = len(patches)

    # get a list of C++'ed NodePatch objects -- NB: This is where we choose
    # whether to employ solid-angle based geometric correction of the dBEM
    # equations
    geo = BEM_USE_GEOMETRIC_CORRECTIONS # (for brevity...)
    from bem_fmm import convert_NodePatch_to_C_Obj as cnvt
    
    c_patches = []
    for de in diffusing_entities:
        de._cpatches = [cnvt(p, geo) for p in de._patches]
        c_patches.extend(de._cpatches)

    # get the RHS vector of the matrix-vector equation 
    # (using peak splitting method -- see Juffer 1990)
    #rhs = get_RHS(diffusing_entities, c_patches, Dext, kappa)
    rhs, fh_peaks = get_fh_peaks(diffusing_entities, c_patches, epsilon, kappa, Dext)
    #fh_peaks = array([0.0 for i in range(2*num_unknowns)])
    BEM_GMRES_Matrix = CUDA_BEM_matrix((2,2),c_patches, kappa, epsilon, array([1.0 for i in range(2*num_unknowns)]))
    
    #import bem_fmm
    #fmm = bem_fmm.BEM_FMM_Tree(patches)
    #BEM_GMRES_Matrix = FMM_BEM_matrix(fmm, c_patches, kappa, epsilon, rhs)

    rhs_offset = BEM_GMRES_Matrix.matvec(fh_peaks)
    rhs_saved = array(rhs)
    #print rhs
    #print rhs_offset
    rhs -= rhs_offset
    #print rhs
    print "fpeak: ", max([abs(r) for r in fh_peaks[:num_unknowns]])
    print "hpeak: ", max([abs(r) for r in fh_peaks[num_unknowns:]])
    BEM_GMRES_Matrix.rhs = array(rhs)
    
    if sum([r*r for r in rhs]) != 0.0:
    
        fvals = rhs[:num_unknowns]
        hvals = rhs[num_unknowns:]
        meanf = sum([abs(f) for f in fvals]) / num_unknowns
        meanh = sum([abs(h) for h in hvals]) / num_unknowns
        
        #b = array([1.0 for i in range(2*num_unknowns)])
        b = array(rhs)
        restart_param = min([num_unknowns, GMRES_RESTART])

        if BEM_PRINT_DEBUG:
            print "starting gmres - %d unknowns" %(num_unknowns)

        x, info = iterative.gmres(BEM_GMRES_Matrix,
                                  b,
                                  x0=array(b),
                                  tol=GMRES_TOLERANCE,
                                  restrt=restart_param,
                                  maxiter=None,
                                  xtype='d')
        
        if info != 0: raise Exception
        if BEM_PRINT_DEBUG:
            print "done gmres - %d iterations" %(BEM_GMRES_Matrix.num_iterations)
            
        #y = BEM_GMRES_Matrix.matvec(x)
        #print y, rhs
        #print max([abs(xx) for xx in (y - rhs)])
                
        x += fh_peaks # add in slowly varying part of the peak-splitting solution
        
    else:
        x = fh_peaks
        
        if BEM_PRINT_DEBUG:
            print "fh peak splitting gave exact solution"
            if epsilon != 1.0: print "HIGHLY SUSPICIOUS!!"
    
    #BEM_GMRES_Matrix.rhs = array(rhs_saved)
    #b = array(rhs_saved)
    #x, info = iterative.gmres(BEM_GMRES_Matrix,
                              #b,
                              #x0=array(x),
                              #tol=GMRES_TOLERANCE,
                              #restrt=restart_param,
                              #maxiter=None,
                              #xtype='d')
            
    #y = BEM_GMRES_Matrix.matvec(x)
    ##print y, rhs_saved
    #print max([abs(xx) for xx in (y - rhs_saved)])
            
    #if BEM_PRINT_DEBUG:
        #print x[:5],x[num_unknowns:num_unknowns+5]
    
    # TODO: clean up this bit...
    
    # set the results at the mesh vertices
    set_results(patches, x)

    # calculate the Electric field (by fitting)
    # this also modifies the normal component of electric field (h)
    # by a fitting process.
    for de in diffusing_entities:
        de.mesh.calculate_electric_fields()
    
    #for v in de.mesh.vertices[:5]:
        #print v.f_kcal_per_mole, v.h_kcal_per_mole_angstrom
        
    return

def solve_energy(diffusing_entities, kappa, Dext, Dint):
    
    # for each entity calculate the electrostatic energy
    for de in diffusing_entities:
    
        # get the patches which make up just this entity
        # HACK: these have been cached as de._patches
        de_c_patches = de._cpatches
        charges, charge_positions = de.universe_charges_and_positions()
        potentials = [v.f for v in de.mesh.vertices]
        normal_derivs = [v.h for v in de.mesh.vertices]
        
        # this gets the reaction field energy contribution
        E = _BEM.calculate_energy(de_c_patches, 
                                  charges, 
                                  charge_positions, 
                                  potentials,
                                  normal_derivs,
                                  kappa, 
                                  Dint, 
                                  Dext)
        
        # add 'self' energy term
        #E += de.mesh.self_energy(kappa, Dint)

        # convert to kJ
        import constants
        E *= constants.elementary_charge*constants.elementary_charge \
          / (constants.epsilon0 * constants.Angstroms * 1000.0)
        
        # store result
        de.energy = E 
        
        # destroy cached patches
        del(de._patches)
        del(de._cpatches)
        
    return
        
def write_fh_state(diffusing_entities, output_filename):
    """Write the electrostatic parameters to file."""

    outfile = open(output_filename,'w')
    for de in diffusing_entities:
        for v in de.mesh.vertices:
            outfile.write("""%e %e\n""" %(v.f, v.h))
        
    outfile.close()
    return

def create_fmm_tree(diffusing_entities):
    """Given a list of diffusing entities, create FMM Octree."""
    
    # iterate over all diffusing entities to get the maximum x,y,z
    # coordinates.
    mins = [None, None, None]
    maxes = [None, None, None]

    for de in diffusing_entities:
        lower_left, upper_right = de.universe_bounding_box
        # cast the vectors into lists
        lower_left = list(lower_left)
        upper_right = list(upper_right)

        # get lowest x,y,z coordinates
        for i, (val, current_lowest) in enumerate(zip(lower_left, mins)):
            if current_lowest is None or val < current_lowest:
                mins[i] = val
        # get highest x,y,z coordinates
        for i, (val, current_highest) in enumerate(zip(upper_right, maxes)):
            if current_highest is None or val > current_highest:
                maxes[i] = val

    # maximum edge length for box
    lower_left_total = Vector(mins)
    upper_right_total = Vector(maxes)
    real_centre = (lower_left_total + upper_right_total) / 2.0
    isotropic_edge_length = max(list(upper_right_total - lower_left_total))

    # create an octree -- unit cube to match FMM code. Use the
    # isotropic_edge_length as a scaling factor for the position of each
    # element.
    import bem_fmm
    root = bem_fmm.BEM_FMM_Tree(isotropic_edge_length)

    # create the node patch list from the list of all diffusing entities
    patches = get_patches_from_entities(diffusing_entities)
    #cull_patches(patches) # cull poorly defined patches

    # iterate over all elements and insert them into the octree
    for node_patch in patches:
        pos = (node_patch.xyz - real_centre) / isotropic_edge_length
        root.insertObject(node_patch, pos)
    root.optimize()
    
    return root

def get_patches_from_entities(diffusing_entities):
    """Create NodePatch objects from the diffusing entities.
    
    (Or more accurately- create NodePatch objects from the mesh objects which
    underlie the diffusing entities."""
    
    # get definitive list of all boundary elements (node patches)
    patch_indices = []
    patch_lengths = []
    patches = []
    for de in diffusing_entities:
        patch_indices.append(len(patches))        
        patch_lengths.append(len(de.mesh.vertices))
        de._patches = de.mesh.construct_node_patches(ref_pt_local=de.centre_of_diffusion,
                                                     ref_pt_universe=de.xyz,
                                                     rotation=de.rotation)

        patches.extend(de._patches)
    return patches

def cull_patches(patches):
    """Remove bad patches from the list of patches.
    
    Exact definition of bad is open to debate and subject to
    modification..."""
    
    total_surface_area = sum([p.area for p in patches])
    average_patch_area = total_surface_area / len(patches)
    
    print "before cull: total area=%e, average area=%e, num patches=%d" \
           %(total_surface_area, average_patch_area, len(patches))
    
    #areas = [(p.area,i) for i,p in enumerate(patches)]
    #areas.sort()
    #cull_list = [area_idx[1] for area_idx in areas[:150]]
    #patches = [p for i,p in enumerate(patches) if i not in cull_list]
    
    #patches = [p for p in patches if p.area > 0.1 * average_patch_area]
    #patches = [p for p in patches if abs(0.5 - p._real_vertex.solid_angle_real) > 0.05]
    
    total_surface_area = sum([p.area for p in patches])
    print "after cull: total area=%e, average area=%e, num patches=%d" \
          %(total_surface_area, total_surface_area / len(patches), len(patches))    
        
    return

def set_results(patches, results):
    """Set the f and h values in results array to vertices of the original mesh."""
    
    num_unknowns = len(patches)
    
    ## set the potential and derivative values in the triangles
    for patch, f, h in zip(patches, results[:num_unknowns], results[num_unknowns:]):
        patch.set_BEM_results(f, h)
        
    return

