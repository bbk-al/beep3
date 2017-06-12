#
# mesh.py
#
# This file contains the code for handling electrostatic meshes.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

# geometry.py defines the Triangle and Vertex objects (amongst other useful
# things).
import geometry
from _BEM import Vector
from vector import Rotation
import math

class MeshedObject(object):
    
    """An electrostatic mesh. Subclasses of this give specific implementations
    of the mesh, from e.g. MSMS or (for test purposes) subdivided
    icosahedra."""
       
    def __init__(self, triangles=[], vertices=[], charges=[]):

        assert( len(charges) == 0)

        self._triangles = triangles[:]
        self._vertices = vertices[:]
        self._charges = charges[:]
        self._patches = None

    @property
    def triangles(self):
        """Return list of Triangle objects in this Mesh."""
        return self._triangles

    @property
    def vertices(self):
        """Return list of Vertex objects in this Mesh."""
        return self._vertices

    @property
    def bounding_radius(self):
        """Return the maximum radius of this meshed object.
        
        Iterate over all vertices and return the greatest length from the
        origin."""
        
        return max([v.length() for v in self.vertices])
        
    @property
    def charges(self):
        """Return charges (in local coordinates)."""
        return self._charges

    def calculate_force(self, kappa, Dext, si_units=False):
        """Returns the electrostatic force on this mesh in local coords, and
        units kT per Angstrom.
        
        If si_units parm set to True, then force is given in Newtons."""

        F = Vector(0.0, 0.0, 0.0)
        for t in self.triangles:
            F += t.force(kappa, Dext)

        # for testing purposes only
        if si_units:
            from constants import force_conversion_Newtons
            return F * force_conversion_Newtons
            
        # convert internal force units to kT per Angstrom
        from constants import force_conversion_kT_per_Angstrom
        Force_kT_per_Angstrom = F * force_conversion_kT_per_Angstrom
        
        return Force_kT_per_Angstrom

    def calculate_moment(self, kappa, Dext):
        """Returns the electrostatic moment on this mesh in local coords.
        
        Forces are scaled to kT per Angstrom, so the moment is in units of
        kT."""

        M = Vector(0.0, 0.0, 0.0)
        for t in self.triangles:
            M += (t.centre - self._centre_of_rotation).cross(t.force(kappa, Dext))

        # convert internal force units to kT per Angstrom (which in
        # calculation for M above, get multiplied by Angstroms anyway (since
        # they're the internal units) to leave moment in units of kT
        from constants import force_conversion_kT_per_Angstrom
        M_kT = M * force_conversion_kT_per_Angstrom
            
        return M_kT

    def calculate_electric_fields(self):
        """Calculate the Electric Field (E) at each Vertex of the Mesh.

        This function assumes that the f and h attributes of each Vertex has
        already been set (probably by calling set_BEM_results on a NodePatch
        object)."""

        for v in self.vertices:
            v.calculate_electric_field()

    def calculate_potential(self, position):
        raise Exception
            
    #def self_energy(self, kappa, Dint):
        #"""Contribution to electrostatic energy from the charges of this entity.
        
        #WARNING: Calculated once and cached."""

        #from numpy import array, zeros
        
        #try:
            #return self._cached_self_energy
        #except AttributeError:
        
            #if kappa==0.0:
                #use_kappa = 1e-8
            #else:
                #use_kappa = kappa
    
            #nlev = 4
            #natoms = len(self.charges) # natoms for consistency with fortran names
            #zat = array([c.position for c in self.charges]).T
            #source_charges = array([c.charge for c in self.charges]).T
            #pot = zeros(natoms,'d').T
            #field = zeros((natoms,3), 'd').T
    
            #import yukawa
            #yukawa.fmmyuk_uni(use_kappa,
                              #zat,
                              #source_charges,
                              #pot,
                              #field,
                              #nlev)
            
            #E = 0.0
            #for ch, p in zip(self.charges, pot):
                #E += ch.charge*p
            #E /= (8.0*Dint*math.pi)

            ## stash the result
            #self._cached_self_energy = E
            
        #return self._cached_self_energy
        
    def init_fh(self, ref_pt_local, ref_pt_universe, rotation):
        """HACK function: just to get f/h peaks for mesh analysis."""

        from geometry import local_to_universal
        cvt = lambda x: VectorPoint(local_to_universal(x, 
                                                       ref_pt_local, 
                                                       ref_pt_universe, 
                                                       rotation))
        
        patches = self.construct_node_patches(ref_pt_local, ref_pt_universe, rotation)
        from bem_fmm import convert_NodePatch_to_C_Obj as cnvt
        from BEM import VectorPoint
            
        c_patches = [cnvt(p, False) for p in patches]
    
        charges = []
        charge_positions = []
        for c in self.charges:
            charges.append(c.charge)
            charge_positions.append(cvt(c.position))
            
        num_unknowns = len(c_patches)
        fh_peaks = [0.0 for xx in range(2*num_unknowns)]
        rhs = [0.0 for xx in range(2*num_unknowns)]
        import _BEM
        _BEM.calculate_fh_peaks(fh_peaks, rhs, charges, charge_positions, c_patches, 40.0, 0.0, 80.0)

        for v,f,h in zip(self.vertices, fh_peaks[:len(patches)],fh_peaks[len(patches):]):
            v.f = f
            v.h = h
            
        return
    
    # this function is rubbish and does not work
    #def calculate_born_solvation_energy(self, kappa, Dext, Dint):
        #"""Returns the electrostatic energy associated with this mesh."""

        #epsilon = Dext / Dint
        
        #import entities
        #from numpy import array, zeros

        ## source charges
        #source_charges = self.charges

        #nlev = 4
        
        ## reaction field charges
        #rf_charges = [entities.Charge(t.centre, t.area*t.f, 0.0)
                      #for t in self.triangles]
        #all_charges = source_charges[:]
        #all_charges.extend(rf_charges)
        
        #natoms = len(all_charges)
        #charge = array([c.charge for c in all_charges]).T
        #zat = array([c.position for c in all_charges]).T
        #pot = zeros(natoms,'d').T
        #field = zeros((natoms,3), 'd').T
        
        #import yukawa
        #pot,field,ier = yukawa.fmmyuk_uni(1e-10,
                                          #zat,
                                          #charge,
                                          #pot,
                                          #field,
                                          #nlev)
        
        ## electrostatic energy associated with the protein charges
        ## interacting with the solvent reaction field (which includes the
        ## influence of the external charge distribution)
        #E = 0.0
        #for ch, p in zip(self.charges, pot):
            #print ch.charge
            #print pot
            #E += ch.charge*p

        ## scale reaction field charges
        ##E *= (epsilon-1)
                
        #E /= 2.0

        ## convert to kT units (i.e. express Energy as a multiple of kT)
        #import constants
        #E *= constants.elementary_charge*constants.elementary_charge \
          #/ (4.0 * math.pi * constants.epsilon0 * 1e-20)
        #E *= 1.0 - (Dint/Dext)
        #E *= constants.Avogadro / 1000.0
        ##E /= constants.kT
        
        #return E
            
    def construct_node_patches(self, ref_pt_local, ref_pt_universe, rotation):
        """Return a list of NodePatch objects for this Mesh.
        
        Required: rotation/translation params for local->universe coordinate
        transformation.

        Iterate over the Vertex objects in this mesh and create NodePatch for
        each; note that the Universe rotation and translation conversion will
        be done by the NodePatch constructor, so the NodePatch objects we
        return are in Universe coordinates."""

        if self._patches is not None:
            return self._patches
        
        # yes, this could probably be done in one line by a list comprehension,
        # but then it wouldn't be very readable.
        patches = []
        for v in self.vertices:
            patch = geometry.NodePatch(v, ref_pt_local, ref_pt_universe, rotation)
            patches.append(patch)
         
        return patches
    
        # store for later use
        #self._patches = patches
        
        #return self._patches
        
    @staticmethod
    def createVertexList(triangle_list):
        """Given a set of triangles, create a Vertex list to represent all
        unique vertices; modify Triangle objects to hold references to
        vertices within the newly created vertex list."""

        vertex_dict = {}
        for t in triangle_list:
            new_verts = []
            for i, vertex in enumerate(t.vertices):
                key = (vertex.x(), vertex.y(), vertex.z())
                if key not in vertex_dict:
                    new_vertex = geometry.Vertex(vertex)
                    vertex_dict[key] = new_vertex
                vertex_dict[key].triangles.append(t)
                new_verts.append(vertex_dict[key])
            t.a, t.b, t.c = new_verts

        return vertex_dict.values()

    def kinemage(self, ref_pt_local, ref_pt_universe, rotation=None, f=None, normals=True, charges=True):
        """Produce Kinemage for this Mesh.

        Pass in a file object to write directly to a file, otherwise this will
        return a string."""

        output = []
        
        # colour list for kinemage
        for hue,prefix in ((0,"red"),(240,"blue")):
            for i, saturation in enumerate(range(0,101,1)):
                colour_name = "%s_%d" %(prefix, i)
                output.append("@hsvcolor {%s} %d %d 100" %(colour_name, hue, saturation))
        
        output.append("@trianglelist {triangles}")
        for t in self.triangles:
            output.append(t.kinemage(centre_of_rotation=ref_pt_local, 
                                     rotation=rotation, 
                                     translation=ref_pt_universe))

        if normals:
            output.append("@vectorlist {normals} off")
            for t in self.triangles:
                output.append(t.kinemageNormal(centre_of_rotation=ref_pt_local, 
                                               rotation=rotation, 
                                               translation=ref_pt_universe))

            output.append("@vectorlist {vertex_normals} color=cyan off")
            for v in self.vertices:
                output.append(v.kinemageNormal(centre_of_rotation=ref_pt_local, 
                                               rotation=rotation, 
                                               translation=ref_pt_universe))

        
        from geometry import local_to_universal
        centre = Vector(0.0,0.0,0.0)
        for atom in self.charges:
            centre += local_to_universal(atom.position, ref_pt_local, ref_pt_universe, rotation)
        centre /= len(self.charges)
        
        ## get maximum radius
        #r = max([(atom.position - centre).length() + atom.radius
                 #for atom in self.charges])                
                
        #output.append("@spherelist {boundary} color=cyan off")
        #output.append("r=%f {} %f %f %f" %(r,centre[0],centre[1],centre[2]))
        
        if charges:
            output.append("@spherelist {charges} color=purple off")
            for c in self.charges:
                output.append(c.kinemage(centre_of_rotation=ref_pt_local,
                                         rotation=rotation, 
                                         translation=ref_pt_universe))

        patches = self.construct_node_patches(ref_pt_local, ref_pt_universe, rotation)
        max_f = max([abs(p._real_vertex.f) for p in patches])
        max_h = max([abs(p._real_vertex.h) for p in patches])
        output.append("@trianglelist {f_node_patches} off")
        for p in patches:
            output.extend(p.kinemage(max_f, 100,f=True))

        output.append("@trianglelist {h_node_patches} off")
        for p in patches:
            output.extend(p.kinemage(max_h, 100,f=False))
            
            
        #output.append("@vectorlist {patch_normals} color=red off")
        #for p in patches:
            #output.extend(p.kinemage_triangle_norms())
            
        #try:
            #output.append("@vectorlist {field directions}")
            #for p in patches:
                #for t in self.triangles:
                    #output.append(t.kinemageFieldDirection(centre_of_rotation=centre_of_rotation,
                                         #rotation=rotation, 
                                         #translation=translation))
        #except ZeroDivisionError:
            #pass

        # turn the whole thing into a string
        result = "\n".join(output) + "\n"

        # try to write it to a file, otherwise return it as a string
        try:
            f.write(result)
        except AttributeError:
            return result

        return

    def expand(self, expansion):
        """Expand the mesh out along normal vectors by specified distance."""
        
        # expand_to_first_hydration_layer means expand the mesh by 3
        # Angstroms in the outward normal direction
        for v in self._vertices:
            v.__setstate__(v + v.normal * expansion)
        
        for t in self.triangles:
            t._force_recalc()

    def write_gts(self, filename):
        """Write the mesh as a GTS (Gnu Triangulated Surface) format file."""
        
        import _GTSMESH as gts
        
        gts_vertices = [gts.SimpleVertex(v[0],v[1],v[2]) 
                        for v in self.vertices]
        gts_triangles = [gts.SimpleTriangle(self.vertices.index(t.a),
                                            self.vertices.index(t.b),
                                            self.vertices.index(t.c)) 
                         for t in self.triangles]
        
        # convert vertex normals to triangle normals
        gts_tnormals = []
        for t in self.triangles:
            n = t.normal
            gts_tnormals.append(gts.SimpleVertex(n[0], n[1], n[2]))
        
        gts.write_gts_surface(gts_vertices, 
                              gts_triangles, 
                              gts_tnormals, 
                              filename)
        return
                        
    def gts_coarsen(self, coarsen_factor, keep_gts_output=None):
        """Coarsen the mesh using the Gnu Triangulated Surface Library."""
        
        import tempfile
        filename_in = tempfile.mktemp(prefix="gtsmesh_", dir=".")
        filename_out = tempfile.mktemp(prefix="gtsmesh_", dir=".")

        # write the mesh to a temporary file
        self.write_gts(filename_in)
        
        import _GTSMESH as gts
        gts.coarsen(filename_in, filename_out, coarsen_factor)
        
        # reload the triangle/vertex objects from the gts file
        self._vertices, self._triangles = GTS_Mesh.load_gtsfile(filename_out)
        
        if (keep_gts_output is not None):
            import shutil
            shutil.copy(filename_out, keep_gts_output)

#        except:
#            print "Failed to coarsen mesh."
            
#        finally:
            
        # cleanup -- delete the two temporary files created above
        from os import remove, path
        from os.path import exists
        if exists(filename_in): remove(filename_in)
        if exists(filename_out): remove(filename_out)

##
## TODO: This class (ProAct_Mesh) needs implementing.
##
#class ProAct_Mesh(MeshedObject):

    #"""A mesh created from ProAct2 Solvent-Accessible-Surface."""

    #def __init__(self):

        ##
        ## TODO: This class
        ##
        #raise(Exception) # this class isn't ready yet!

        ## call the MeshedObject __init__ to set various internal lists.
        #super(ObjMesh, self).__init__([],[],[])

        ## load the proact solvent accessible surface
        ##proact = open('output-3.proact','r')
        ##for l in proact.readlines():
        ##    [id,x,y,z] = l.split()
         ##   c.append([float(x),float(y),float(z)])
        ##    r.append(2)
        ##proact.close()

class GTS_Mesh(MeshedObject):

    """GTS format mesh."""

    def __init__(self, gts_filename):
        
        # call the MeshedObject __init__ to set various internal lists.
        super(GTS_Mesh, self).__init__([],[],[])
        
        self._vertices, self._triangles = self.load_gtsfile(gts_filename)

        return

    def enforce_no_charge_clashes(self, atomlist, keep_gts_output=None):
        something_moved = False
        
        print "fixing vertex/charge clashes"
        for v in self._vertices:
            
            while True:
                unclean = False
                for at in atomlist:
                    
                    while True:
                        
                        dist = (v - at.position).length()
                        
                        if dist >= at.radius: 
                            break
                        
                        unclean = True
                        something_moved = True
                        #print dist, at.radius
                        #print "vertex within atom radius!"
                        direction = ((v-at.position)/dist + v.normal)/2.0
                        offset = v + direction*1e-3
                        #print offset - v

                        v.x = offset.x
                        v.y = offset.y
                        v.z = offset.z
                        for t in v.triangles:
                            t._force_recalc()

                # loop until all vertices missing all atoms
                if not unclean:
                    break
        if keep_gts_output is not None:
            self.write_gts(keep_gts_output)
            
        return something_moved
    
    @staticmethod
    def load_gtsfile(filename):
        
        from gts_utils import get_next_line_not_comment    
        f = open(filename,"r")    
        
        vertices, triangles = [], []
        try:
        
            # get the first line, ignoring comment lines (starting with # or !)
            first_line = get_next_line_not_comment(f)
                
            # number of vertices, edges and faces
            nv, ne, nf = [int(xx) for xx in first_line.split()[:3]]
            #print nv, ne, nf
            while (len(vertices) < nv):
                
                # convert to Vertex object
                vx, vy, vz = [float(xx) for xx in get_next_line_not_comment(f).split()]
                vertices.append(geometry.Vertex(vx, vy, vz))
    
            # edges are defined as an ordered pair of vertices
            edges = []
            while (len(edges) < ne):
                edge_verts = get_next_line_not_comment(f).split()
                edges.append( set([int(edge_verts[0])-1, int(edge_verts[1])-1]) )
            
            # triangle defs come after the edge defs
            while (len(triangles) < nf):
                e1, e2, e3 = [edges[int(xx)-1] for xx in get_next_line_not_comment(f).split()]
                
                # get the 3 vertices, in order
                v1_idx = e1.intersection(e2).pop()
                v2_idx = e2.intersection(e3).pop()
                v3_idx = e3.intersection(e1).pop()
                
                new_t = geometry.Triangle(vertices[v1_idx],vertices[v2_idx],vertices[v3_idx])
                for v in new_t.vertices:
                    v.triangles.append(new_t)
                triangles.append(new_t)
    
        except IOError:
            print "Bad GTS file."
            vertices, triangles = [], []
    
        finally:
            f.close()
        
        return vertices, triangles    
        
class MeshedSphere(MeshedObject):

    """A spherical surface defined by a set of triangles. Subclasses of this
    give implementations in terms of subdivided octahedron or icosahedron.
    
    This is mainly intended for testing the Boundary Element Method
    electrostatics code."""

    def __init__(self,
                 radius=1.0,
                 xyz=Vector(0.0,0.0,0.0),
                 rot=None,
                 num_subdivides=3):

        # call the MeshedObject __init__ to set various internal lists.
        super(MeshedSphere, self).__init__([],[],[])

        self.radius = radius
        self.centre = xyz

        num_vertices = len(self.vertices)
        num_triangles = len(self.triangles)

        import _GTSMESH as gts
        import tempfile
        spherical_file = tempfile.mktemp(suffix=".gts", 
                                         prefix="gts_sphere_", 
                                         dir=".")
        gts.write_spherical_mesh(num_subdivides, spherical_file)
        self._vertices, self._triangles = GTS_Mesh.load_gtsfile(spherical_file)
        
        # cleanup
        import os
        os.remove(spherical_file) 
        
        print num_vertices, num_triangles

        analytic_vol = 4.0/3.0 * math.pi
        analytic_sur = 4.0 * math.pi
        
        # invent a random rotation if one wasn't passed in
        if rot is None:
            rot = geometry._rand_rot()
        
        # apply rotation, set radius, and xyz offset
        for v in self._vertices:
            pt = geometry.apply_quaternion_to_vector(rot, v).normal()
            v.__setstate__(pt*radius + xyz)

        # update triangle centres
        for t in self._triangles:
            t._centre = t._calc_centre()

        # HACK: set triangle normals as same as vector from centre to midpoint
        # of triangles
        for t in self.triangles:
            t._normal_vector.__setstate__((t._centre - xyz).normal())
            t._area = t._calc_area()

        return

#
# TESTS
#
if __name__ == "__main__":
        
    def born_ion():
        
        import entities
        from vector import Vector, Rotation, Quaternion
    
        rad = 3.0
        detail = 5 
        Dext = 80.0
        Dint = 1.0
        charge = 1.0
        kappa = 0.0
        dist = 3.0
        #kappa = 0.0728 # 50mM salt; 13.7A screening length
        #kappa = 1.0 / 6.142
        
        origin = Vector(0.0, 0.0, 0.0) # the origin
    
        from mesh import MeshedSphere
        from math import pi
    
        # create diffusing entity
        de = entities.DiffusingEntity()
        import universe
        #de.set_xyz(universe._rand_xyz())
        de.set_xyz(Vector(-dist/2.0,0.0,0.0))
        de.atoms = [entities.Atom(origin, rad)] 
        de.mesh = MeshedSphere(radius=rad, num_subdivides=detail)
        de.mesh._charges = [entities.Charge(origin, charge, radius=rad)]

        de2 = entities.DiffusingEntity()
        de2.set_xyz(Vector(dist/2.0,0.0,0.0))
        de2.atoms = [entities.Atom(origin, rad)] 
        de2.mesh = MeshedSphere(radius=rad, num_subdivides=detail)
        de2.mesh._charges = [entities.Charge(origin, 2*charge, radius=rad)]
        
        import constants
        premult = -constants.Avogadro * (constants.elementary_charge ** 2) \
                / (1000.0*8.0*math.pi*constants.epsilon0*constants.__Angstroms)
        #analytic = -700.7725 * (1.0 - 1.0/Dext) * (charge*charge/rad) #Nathan's
        analytic = premult * (1.0/Dint - 1.0/Dext) * (charge*charge/rad) 

        
        import BEM
        BEM.solveElectrostatics([de],kappa,Dext,Dint)
        BEM.solve_energy([de],kappa,Dext,Dint)

        # write a kinemage
        f = open("born.kin", 'w')
        print >>f, "@kinemage"
        de.mesh.kinemage(ref_pt_local=de.centre_of_diffusion, ref_pt_universe=de.xyz, f=f)
        #de2.mesh.kinemage(f, translation=de2.xyz)
        f.close()        
        
        print "BORN SOLVATION ENERGY (kJ/mol): ", \
              de.energy*constants.Avogadro , analytic, 100.0 * (de.energy*constants.Avogadro - analytic) / analytic
        print "BORN SOLVATION ENERGY (kCal/mol): ", \
              de.energy*constants.Avogadro*0.2388, analytic * 0.2388

        print de.mesh.calculate_force(kappa,Dext,si_units=True) #* constants.kT * 0.2388 * constants.Avogadro
        #print de2.mesh.calculate_force(kappa,Dext,si_units=True)
        
        #de.mesh.write_gts("born_ion.gts")
        #de.mesh.gts_coarsen(0.9)
        
        return

    def some_molecule(filename):
        
        import entities
        from vector import Vector, Rotation, Quaternion
    
        Dext = 80.0
        Dint = 2.0
        kappa = 0.0
        
        origin = Vector(0.0, 0.0, 0.0) # the origin
        
        #from Scientific.IO.PDB import Structure
        #struct = Structure(filename + ".pdb")
        
        from geometry import _rand_rot
        
        # create diffusing entity
        de = entities.DiffusingEntity()
        de.pqr_filename = filename + ".pqr"
        #de.set_rotation(_rand_rot())
        de.set_charges_from_pqr(filename + ".pqr", set_atoms_from_pqr=True)

        # gts mesh must already exist
        de.create_mesh(filename)
        
        centre = Vector(0.0,0.0,0.0)
        for atom in de.atoms:
            centre += atom.position
        de.set_centre_of_diffusion( centre / len(de.atoms) )
        de.set_xyz(de.centre_of_diffusion)

        import BEM
        import tempfile
        from shutil import move
        #while(len(de.mesh.vertices) > 100):
        for dummy in range(1):
        
            v = len(de.mesh.vertices)
                    
            BEM.solveElectrostatics([de],kappa,Dext,Dint)
            BEM.solve_energy([de],kappa,Dext,Dint)
    
            f = open("%s-%d.kin" %(filename,v),'w')
            print >>f, "@kinemage"
            de.mesh.kinemage(ref_pt_local=de.centre_of_diffusion, ref_pt_universe=de.xyz, rotation=de.rotation, f=f)
            f.close()
            
            import constants
            print "SOLVATION ENERGY (kJ/mol): ", de.energy*constants.Avogadro
            print "SOLVATION ENERGY (kCal/mol): " , de.energy*constants.Avogadro*0.2388

            fresults = open("%s.txt" %(filename), "a")
            print >>fresults, "%d %f" %(v, de.energy*constants.Avogadro)
            fresults.close()

            # write out the converged electrostatic values -- might be useful...
            fh_vals = open("%s-%d.fh"%(filename, len(de.mesh.vertices)), "w")
            for v in de.mesh.vertices:
                
                # local and universe coords should coincide in these cases
                uni_v = x,y,z = de.convert_local_coordinate(v)
                assert((uni_v-v).length() < 1e-12)
                print >>fh_vals, "%f %f %f %9.9e %9.9e" %(x, y, z, v.f, v.h)
            fh_vals.close()
            
            # coarsen mesh - keep the gts file for later analysis
            #gts_filename = tempfile.mktemp(prefix="%s-gts_mesh" %(filename), dir=".")
            #de.mesh.gts_coarsen(200)
            #while de.mesh.enforce_no_charge_clashes(de.atoms, keep_gts_output=gts_filename):
            #    de.mesh.gts_coarsen(0)
            #move(gts_filename, "%s-%d.gts" %(filename, len(de.mesh.vertices)))
            
        return
    
    #import cProfile
    #cProfile.run('born_ion()','born_ion_profile')
    
    born_ion()
    import sys
    some_molecule(sys.argv[1])
    #some_molecule("ubiquitin/ubi_test1")
    #some_molecule("1MAH-ache")
    #some_molecule("2CGA-mono")
    #some_molecule("ecm_model")
    #some_molecule("2OZO")
