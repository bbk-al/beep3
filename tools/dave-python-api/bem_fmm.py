from numpy import zeros,matrix
import octree
import _BEM
from BEM import VectorPoint
import math
import fmm
from vector import Vector

class BEM_FMM_Tree(octree.Octree):

    """Octree structure for the purposes of BEM/FMM."""

    def __init__(self, patches):

        from Scientific.Geometry import Vector
        
        # process patches to figure out scaling factor
        mins, maxs = self.get_bounds(patches)
        self.box_scale_factor = max([maxs[coord] - mins[coord] 
                                     for coord in range(3)])
        centre = (Vector(maxs) + Vector(mins)) / 2.0

        # delegate to parent constructor; force unit cube (to match FMM code).
        super(BEM_FMM_Tree, self).__init__(worldCentre=centre,
                                           worldSize=1.0,
                                           max_objects_per_node=25)

        for p in patches:

            #not much point in this assertion because the octree will throw
            #and exception if the patch is outside this boundary anyway
            #assert( (p.xyz / self.box_scale_factor).length() < 0.5*math.sqrt(3.0) )
            
            # insert patch into the FMM tree
            obj = self.insertObject(p, centre + ((p.xyz - centre)/self.box_scale_factor))

        self.optimize()
            
        # instantiate the Fast Multipole Method class
        self.fmm = fmm.FMM(self)
            
    @staticmethod
    def get_bounds_low_mem(patches):
        """Get the maximum and minimum bounds of set of patches, without
        requiring possibly large lists."""

        mins = [None,None,None]
        maxs = [None,None,None]
        
        for p in patches:
            for coord in (0,1,2): # loop over x,y,z
                val = p.xyz[coord]
                if mins[coord] is None or val < mins[coord]:
                    mins[coord] = val
                if maxs[coord] is None or val > maxs[coord]:
                    maxs[coord] = val
        return mins, maxs

    @staticmethod
    def get_bounds(patches):
        """Get the maximum and minimum bounds of set of patches.
        
        Moderately more elegant, without regard for memory consumption. Could
        we use a generator in the min/max function?"""

        mins, maxs = [0,0,0],[0,0,0]
        for coord in (0,1,2):
            maxs[coord] = max([p.xyz[coord] 
                               for p in patches])
            mins[coord] = min([p.xyz[coord] 
                               for p in patches])

        return mins, maxs
    
    def hack_set_charges(self, charges):
        """A blatant hack into the underlying octree to reset the magnitudes of
        the charges within the FMM octree.
        
        This function is required because in BEM/FMM method we run the FMM 8
        times -- for equivalent charges f*nx, f*ny, f*nz and h; for kappa=0
        and kappa=<actual_val>."""
        
        assert(len(charges) == len(self.master_list))
        for new_charge,wrapped_obj in zip(charges, self.master_list):
            wrapped_obj.charge = new_charge
            
        return

    @property
    def num_unknowns(self):
        """Return the number of unknowns in the BEM/FMM system (num
        elements)."""

        # return the total number of elements stored in the octree.
        return len(self.master_list)

    def matvec(self, f_h_vector, c_patches, kappa, epsilon):

        import time
        fmm_start = time.clock()
        fmm_res = self.do_fast_multipoles(f_h_vector, kappa, epsilon)
        fmm_time = time.clock() - fmm_start
        local_start = time.clock()
        local = self.do_locals(f_h_vector, c_patches, kappa, epsilon)
        local_time = time.clock() - local_start
        
        # for testing -- compare FMM to explicit calcs
        #explicit = self.do_farfield_explicitly(f_h_vector, c_patches, 
        #                                       kappa, epsilon)
        
        result = local+fmm_res
        #print local[:10],local[-10:]        
        #print fmm_res[:10],fmm_res[-10:]        
        #print explicit[:10],explicit[-10:]        
        #print result[:10],result[-10:]

        return result, fmm_time, local_time

    def do_fast_multipoles(self, f_h_vector, kappa_actual, epsilon):
        """Run the fast multipole method.

        Runs the fast multipole method for each of the 8 charge distributions
        required by the BEM. See Lu2007."""
        
        # this takes care of case when kappa=0
        # since kappa=0 leads to singular results from the FMM
        # code, so approximate zero with a very small number
        # such as 1e-8
        if kappa_actual < 1e-9:
            kappa = 1.0e-9
        else:
            kappa = kappa_actual

        scaled_kappa = kappa * self.box_scale_factor
        scaled_kappa0 = 1.0e-9 * self.box_scale_factor

        num_unknowns = self.num_unknowns
        results = zeros(2*num_unknowns,'d')
        nlev = self.max_depth

        patches = [wrapped_data.obj for wrapped_data in self.master_list]

        # lists to hold the effective charges
        charges_nx_f = []
        charges_ny_f = []
        charges_nz_f = []
        charges_h = []
        #locations = zeros((3,num_unknowns),'d')

        for idx, wrapped_data in enumerate(self.master_list):

            np = wrapped_data.obj
            #locations[:,idx] = wrapped_data.position

            a = np.area
            nx, ny, nz = np.normal

            f = f_h_vector[idx]
            h = f_h_vector[idx + num_unknowns]

            charges_nx_f.append(a * f * nx)
            charges_ny_f.append(a * f * ny)
            charges_nz_f.append(a * f * nz)
            charges_h.append(a * h)

        # utility functions
        def run_yukawa_fmm(charges, beta):
            self.hack_set_charges(charges)
            self.fmm.solve(beta, calc_local_interactions=False)
            return
        
        def update_fcomps(field_component, coefficient):
            for idx, wrapped_data in enumerate(self.master_list):
                field = wrapped_data.field
                results[idx] += coefficient*field[0,field_component] / \
                                 self.box_scale_factor

        def update_hcomps(field_component, coefficient):
            for idx, wrapped_data in enumerate(self.master_list):
                field2 = wrapped_data.field2
                nx, ny, nz = wrapped_data.normal
                hcomp = field2[0,field_component] * nx
                hcomp += field2[1,field_component] * ny
                hcomp += field2[2,field_component] * nz
                results[idx+num_unknowns] += coefficient*hcomp / \
                    (self.box_scale_factor*self.box_scale_factor)
                
        # scaled kappa because this is all in a unit cube, rather than the
        # real units of the universe.
        #print "starting FMM calculations (%d unknowns)" %(num_unknowns)
        run_yukawa_fmm(charges_nx_f, scaled_kappa)
        update_fcomps(0,1.0)
        update_hcomps(0,1.0/epsilon)
        run_yukawa_fmm(charges_ny_f, scaled_kappa)
        update_fcomps(1,1.0)
        update_hcomps(1,1.0/epsilon)
        run_yukawa_fmm(charges_nz_f, scaled_kappa)
        update_fcomps(2,1.0)
        update_hcomps(2,1.0/epsilon)
        run_yukawa_fmm(charges_h, scaled_kappa)

        for idx, wrapped_data in enumerate(self.master_list):

            field = wrapped_data.field
            np = wrapped_data.obj
            nx, ny, nz = np.normal
            
            hcomp = (field[0,0]*nx + field[0,1]*ny + field[0,2]*nz) / epsilon
            results[idx] += wrapped_data.potential 
            results[idx+num_unknowns] += hcomp / self.box_scale_factor

        run_yukawa_fmm(charges_nx_f, scaled_kappa0)
        update_fcomps(0,-1.0/epsilon)
        update_hcomps(0,-1.0/epsilon)
        run_yukawa_fmm(charges_ny_f, scaled_kappa0)
        update_fcomps(1,-1.0/epsilon)
        update_hcomps(1,-1.0/epsilon)
        run_yukawa_fmm(charges_nz_f, scaled_kappa0)
        update_fcomps(2,-1.0/epsilon)
        update_hcomps(2,-1.0/epsilon)
        run_yukawa_fmm(charges_h, scaled_kappa0)
        
        for idx, wrapped_data in enumerate(self.master_list):

            field = wrapped_data.field
            np = wrapped_data.obj
            nx, ny, nz = np.normal
            
            hcomp = field[0,0]*nx + field[0,1]*ny + field[0,2]*nz
            results[idx] -= wrapped_data.potential            
            results[idx+num_unknowns] -= (hcomp / self.box_scale_factor)
            
        sf = 4.0*math.pi*self.box_scale_factor
        results /= sf
        
        #print "done local expansions calculations"
        
        return results

    def do_farfield_explicitly(self, f_h_vector, c_patches, kappa, epsilon):
        """Calculate far-field contributions to the matrix-vector product.
        
        This is a test function intended to verify the results of the faster
        FMM method."""
        
        num_unknowns = self.num_unknowns
        results = zeros(2*num_unknowns,'d')

        patches = [wrapped_data.obj for wrapped_data in self.master_list]
        reverse_lookup = {} # lookup table from element object to index into f/h vector
        for i, p in enumerate(patches):
            reverse_lookup[p] = i

        print "Starting explicit far-field interactions."

        # loop over the leaf nodes, and calculate the local neighbour
        # interactions for each element within each node
        for leaf in self.leaves:

            # this gets us a list of NodePatch objects
            patches_in_leaf = [wrapped_data.obj for wrapped_data in leaf._data]

            # this is a list of neighbouring leaves (everything else is
            # far-field by definition)
            neighbour_patches = []
            for nn in leaf.colleagues:
                some_patches = [wrapped_data.obj for wrapped_data in nn._data]
                neighbour_patches.extend(some_patches)

            for source_patch in patches_in_leaf:
                assert(source_patch in neighbour_patches)
 
                idx = reverse_lookup[source_patch]
                source_patch_cpp = c_patches[idx]

                accum_f = 0.0
                accum_h = 0.0
                for target_patch in patches:
                    
                    # reject near-field
                    if target_patch in neighbour_patches:
                        continue

                    # get the NodePatch as a C++ object for the _BEM.so lib
                    idx_target = reverse_lookup[target_patch]
                    target_patch_cpp = c_patches[idx_target]
                    
                    f = f_h_vector[idx_target]
                    h = f_h_vector[idx_target + num_unknowns]

                    # C++ library to do the heavy lifting
                    result = [0.0,0.0,0.0,0.0]
                    _BEM.calculateMatrixElements_ff(result, 
                                                    kappa, 
                                                    epsilon,
                                                    source_patch_cpp,
                                                    target_patch_cpp)

                    A,B,C,D = result
                    
                    # this function should never hit the diagonals
                    assert(idx_target != idx)

                    accum_f += f*B + h*A
                    accum_h += f*D + h*C

                # store the results
                results[idx] = accum_f
                results[idx + num_unknowns] = accum_h

        print "Completed explicit far-field interactions."

        return results        
    
    def do_locals(self, f_h_vector, c_patches, kappa, epsilon, all=False):
        """Calculate the local interactions between all patches.

        Returns a vector of length 2*num_unknowns containing the
        near-neighbour interactions."""

        num_unknowns = self.num_unknowns
        results = zeros(2*num_unknowns,'d')

        magic_epsilon_const = (0.5 + 1.0/(2.0*epsilon))

        patches = [wrapped_data.obj for wrapped_data in self.master_list]
        reverse_lookup = {} # lookup table from element object to index into f/h vector
        for i, p in enumerate(patches):
            reverse_lookup[p] = i

        #print "Starting local interactions."

        # loop over the leaf nodes, and calculate the local neighbour
        # interactions for each element within each node
        if all:
            targs = [self]
        else:
            targs = self.leaves
        for leaf in targs:

            # this gets us a list of NodePatch objects
            patches_in_leaf = [wrapped_data.obj for wrapped_data in leaf._data]

            # this is a list of neighbouring leaves
            neighbour_patches = []
            for nn in leaf.colleagues:
                some_patches = [wrapped_data.obj for wrapped_data in nn._data]
                neighbour_patches.extend(some_patches)

            for source_patch in patches_in_leaf:
                assert(source_patch in neighbour_patches)
 
                idx = reverse_lookup[source_patch]
                source_patch_cpp = c_patches[idx]

                accum_f = 0.0
                accum_h = 0.0
                assert(not all or len(neighbour_patches) == num_unknowns)
                for target_patch in neighbour_patches:

                    # get the NodePatch as a C++ object for the _BEM.so lib
                    idx_target = reverse_lookup[target_patch]
                    target_patch_cpp = c_patches[idx_target]
                    
                    f = f_h_vector[idx_target]
                    h = f_h_vector[idx_target + num_unknowns]

                    # C++ library to do the heavy lifting
                    result = [0.0,0.0,0.0,0.0]
                    _BEM.calculateMatrixElements(result, 
                                                 kappa, 
                                                 epsilon,
                                                 source_patch_cpp,
                                                 target_patch_cpp)

                    A,B,C,D = result
                    
                    # take care of the diagonal magic epsilon constant
                    if idx_target == idx:

                        gc = source_patch.solid_angle / (4.0*math.pi)
                        cg = 1.0 - gc
                        fcorrect = cg + (gc / epsilon)
                        hcorrect = gc + (cg / epsilon)
                        B += fcorrect
                        C += hcorrect
                        #B += magic_epsilon_const
                        #C += magic_epsilon_const

                    accum_f += f*B + h*A
                    accum_h += f*D + h*C
#                    if not all: print idx, idx_target, C, h, h*C, accum_h

                # store the results
                results[idx] = accum_f
                results[idx + num_unknowns] = accum_h

        #print "Completed local interactions."

        return results

def convert_NodePatch_to_C_Obj(python_node_patch, geometry_correction=True):
    """Convert a NodePatch object into the C++ version."""

    from geometry import NodePatch
    assert(isinstance(python_node_patch, NodePatch))

    # if geometry correction is turned on we use a real solid angle in the BEM
    # calculations, rather than the usual assumption of 0.5 (corresponding to
    # a planar dielectric boundary).
    if geometry_correction:
        geometry_term = python_node_patch.solid_angle / (4.0*math.pi)
        central_vertex = python_node_patch.xyz
    else:
        geometry_term = 0.5
        central_vertex = python_node_patch.xyz
    
    cpp_np = _BEM.NodePatch(VectorPoint(central_vertex),
                            VectorPoint(python_node_patch.normal),
                            geometry_term)

    for t in python_node_patch.triangles:
        c_triangle = _BEM.Triangle(VectorPoint(t.a),
                                   VectorPoint(t.b),
                                   VectorPoint(t.c),
                                   VectorPoint(t.normal),
                                   VectorPoint(t.centre),
                                   t.area)
        cpp_np.addTriangle(c_triangle)

    # done converting NodePatch to _BEM.NodePatch
    return cpp_np
