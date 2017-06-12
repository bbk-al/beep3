#!/usr/bin/python
# -*- coding: utf-8 -*-
import _BEM
from _BEM import Vector
import constants
import pqrtools
from numpy import zeros, array, ones, matrix, linalg
from scipy.linalg import iterative
import time
import os
import geometry
import gts_utils
from math import pi
import string
GMRES_PRINT_DEBUG = True

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

    # generator function to provide node patch iterator
    @property
    def vertices(self):
        for x in xrange(self.num_vertices):
            yield self.get_vertex(x)
            
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
                sph = self.get_spherical_harmonic_approx_h(self.entity_mappings[k])
                flux = self.calculate_net_flux(self.entity_mappings[k], epsilon, sph)
                print "Flux on mesh %d = %f" %(k, flux)
        except AttributeError:
            pass

        return

    def calc_energies(self, Dext, Dint, kappa):

        epsilon = Dext / Dint
        
        import constants
        con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
          / (constants.epsilon0 * constants.Angstroms * 1000.0)

        try:

            for k in self.entity_mappings.keys():
                energy = self.calculate_energy(self.entity_mappings[k],
                                               kappa,
                                               Dint,
                                               Dext) * con
                print "Solvation energy for mesh %d = %f" %(k, energy)

        except AttributeError:
            # for meshes with no entity mapping of object to node patch range
            energy = self.calculate_energy_whole_mesh(kappa, Dint, Dext) * con
            print "Solvation energy = %f" %(energy)

        return

    def get_spherical_harmonic_approx_for_molecule(self, k):
        
        spf = self.get_spherical_harmonic_approx_f(self.entity_mappings[k])
        sph = self.get_spherical_harmonic_approx_h(self.entity_mappings[k])
        
        return spf, sph

    def get_spherical_harmonic_approx_f(self, offs):
        
        locs = []
        fvals = []

        start = offs.vertex_numbering_offset
        end = start + offs.num_vertices
        patch_list = [self.mesh.get_node_patch(x) for x in range(start,end)]
        
        import _SPHARM as sp
        locs = [np.node() for np in patch_list]
        fvals = [np.f for np in patch_list]
        sph_f = sp.least_squares_spherical_harmonic(locs, fvals, offs.origin, NUM_SPHARM_TERMS)
        
        return sph_f

    def get_spherical_harmonic_approx_h(self, offs):
        
        locs = []
        hvals = []

        start = offs.vertex_numbering_offset
        end = start + offs.num_vertices
        patch_list = [self.mesh.get_node_patch(x) for x in range(start,end)]
        
        import _SPHARM as sp
        locs = [np.node() for np in patch_list]
        hvals = [np.h for np in patch_list]
        sph_h = sp.least_squares_spherical_harmonic(locs, hvals, offs.origin, NUM_SPHARM_TERMS)
        
        return sph_h
    
    def calc_reaction_field_forces(self, Dext, Dint, kappa):
        
        import constants
        con = constants.Avogadro*constants.elementary_charge*constants.elementary_charge \
          / (constants.epsilon0 * 1000.0)

        for k in self.entity_mappings.keys():
            
            # get a spherical harmonic approximation for the f and h vals
            sph_f, sph_h = self.get_spherical_harmonic_approx_for_molecule(k)
            
            rf_force = self.calc_reaction_field_force(self.entity_mappings[k],
                                                      kappa,
                                                      Dint,
                                                      Dext,
                                                      sph_h)
            
            print "Net reaction field force on charges for mesh %d = %s" %(k, rf_force)

        return rf_force

    def boundary_force(self, Dext, Dint, kappa):
        """Calculate net force over all triangular elements in the mesh."""

        for k in self.entity_mappings.keys():
            
            # get a spherical harmonic approximation for the f and h vals
            sph_f, sph_h = self.get_spherical_harmonic_approx_for_molecule(k)
            
            boundary_force = self.calc_boundary_force(self.entity_mappings[k],
                                                      Dint,
                                                      Dext,
                                                      kappa, 
                                                      sph_f, 
                                                      sph_h)
            print "Boundary force on mesh %d = %s" %(k, boundary_force)

        return boundary_force
 
class StaticEnsemble(CPP_Mesh):

    def __init__(self, config_filename):

        # call base class constructor
        super(StaticEnsemble, self).__init__()

        from config_file import ConfigFile
        cfg = ConfigFile(config_filename)

        self.entity_mappings = {}

        from gts_utils import get_vertices_triangles_from_gts

        for ctr, sim_def in enumerate(cfg.sim_defs):

            self.entity_mappings[ctr] = _BEM.Offsets()
            offs = self.entity_mappings[ctr]
            offs.origin = sim_def.xyz_offset # this is the centre of the molecule in universe coords
            
            gts_filename = sim_def.mesh_def.gts_file
            centre_filename = sim_def.mesh_def.centre_file
            vertices, triangles = get_vertices_triangles_from_gts(gts_filename)
            offs.vertex_numbering_offset = self.num_vertices
            offs.num_vertices = len(vertices)

            x,y,z = [string.atof(xx)
                     for xx in open(centre_filename,'r').readline().split()]
            offset = sim_def.xyz_offset - Vector(x,y,z)

            for v in vertices:
                self.add_vertex(v + offset)

            offs.triangle_numbering_offset = self.num_triangles
            offs.num_triangles = len(triangles)

            for t in triangles:
                self.define_triangle(t[0] + offs.vertex_numbering_offset,
                                     t[1] + offs.vertex_numbering_offset,
                                     t[2] + offs.vertex_numbering_offset)

            charge_file = open(sim_def.mesh_def.xyzq_file,'r').readlines()
            charge_list = []
            for line in charge_file:
                x,y,z,q,r = [string.atof(xx) for xx in line.split()]
                new_xyz = Vector(x,y,z) + offset
                charge_list.append(_BEM.Charge(q, new_xyz, r))
            offs.charge_numbering_offset = self.num_charges
            offs.num_charges = len(charge_list)
            
            self.add_charges(charge_list)

            # load fh values
            try:
                fh_file = open(sim_def.mesh_def.fh_vals_file, 'r').readlines()
                #vertex_list = [v for v in self.vertices]
                for line, v in zip(fh_file, self.vertices):
                    f,h = [string.atof(xx) for xx in line.split()]
                    v.f = f
                    v.h = h
            except:
                pass
            
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
    #print "fh peaks:", fh_peaks
    h_peaks = [np.h for np in mesh.node_patches]
    #print "h peaks:" , h_peaks
    fh_peaks.extend(h_peaks)

    rhs = [np.rhs_f for np in mesh.node_patches]
    #print "rhsf:", rhs
    rhsh = [np.rhs_h for np in mesh.node_patches]
    #print "rhsh:", rhsh
    rhs.extend(rhsh)
    #print rhsf
    return array(rhs), array(fh_peaks)

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

    @property
    def shape(self):
        return (2*self.N,2*self.N)

    @property
    def dtype(self):
        return 'd'

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

        self.mesh.integrate_node_patches(self.epsilon, self.kappa)

        if GMRES_PRINT_DEBUG:
            print "MultiplyVector: %f" %(time.clock() - start)

        for i,np in enumerate(self.mesh.node_patches):
            res[i] = np.f_accum
            res[i+self.N] = np.h_accum

        return array(res)

