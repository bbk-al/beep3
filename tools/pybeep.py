# -*- coding: utf-8 -*-
from ctypes import *
cpp = CDLL("libstdc++.so.6", mode=RTLD_GLOBAL | RTLD_LOCAL)
from libBEEP import _Charge, _Mesh, _MeshInstance, _BEEP, BasicTriangle, _Vector, Quaternion, Octree, NodePatchTree

from string import atof

class Charge(_Charge, object):
    
    def __init__(self, *args):
        super(Charge, self).__init__(*args)    

    def replicate_in_universe_coords(self, xyz_local, rotation, xyz_universe):
        """A factory method to produce new charges in a new coordinate frame."""
        
        new_charge = Charge(self.position(), self.charge, self.radius)
        new_charge._change_coordinate_frame(xyz_local, rotation, xyz_universe)
        return new_charge

    @staticmethod
    def getChargesFromFile(filename):
        """Returns a list of Charge objects from an xyzqr file."""
        
        charge_list = []
        f = open(filename, 'r')
        for line in f.readlines():
            try:
                bits = [atof(xx) for xx in line.split()]
                charge_list.append(Charge(Vector(bits[0],bits[1],bits[2]), bits[3], bits[4]))
            except:
                pass
        f.close()
        
        return charge_list

    def __str__(self):
        print self.position()
        
class Vector(_Vector, object):
    
    def __init__(self, *args):
        super(Vector, self).__init__(*args)    

    def normalised(self):
        cp = Vector(self)
        cp.normalise()
        return cp

    def __add__(self, *args):
        return Vector(super(Vector, self).__add__(*args))
    def __sub__(self, *args):
        return Vector(super(Vector, self).__sub__(*args))

    def __iter__(self):
        ctr=0
        while (ctr < 3):
            if (ctr == 0): yield self.x
            if (ctr == 1): yield self.y
            if (ctr == 2): yield self.z
            ctr += 1


class Mesh(_Mesh, object):

    def __init__(self, *args):
        super(Mesh, self).__init__(*args)

        pts = NodePatchTree(10, self.get_centre(), self.get_radius()*3.0)
        for xx in self.node_patches:
            pts.insert(xx)

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

    # generator function to provide vertex iterator
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
            a = t.v1()
            b = t.v2()
            c = t.v3()
            colour = randomColour()
            output.append("X %f %f %f %f %f %f %s %f %f %f" %(a.x, a.y, a.z,
                                                   b.x, b.y, b.z, colour,
                                                   c.x, c.y, c.z ))

        output.append("@vectorlist {centres} off")
        for t in self.triangles:
            a = t.centre()
            b = t.centre() + t.normal()
            output.append("{}P %f %f %f %f %f %f" %(a.x,a.y,a.z,b.x,b.y,b.z))

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
            vn = v + np.normal()
            output.append("{}P %f %f %f %f %f %f" %(v.x, v.y, v.z, vn.x, vn.y, vn.z))

        output.append("@vectorlist {alternative_normals} off")
        for np in self.node_patches:
            v = np.node()
            v_ftc = v + np.alt_normal()
            output.append("{}P %f %f %f %f %f %f" %(v.x, v.y, v.z, v_ftc.x,v_ftc.y,v_ftc.z))

        output.append(self.kinemage_node_patches())

        output.append("@spherelist {charges} color=purple off")
        for ch in self.charges:
            x,y,z = Vector(ch.position())
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

            print "Net reaction field force on charges for mesh %d = %s" %(k, rf_force )

            return rf_force
            #yield self.get_triangle(x)

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

    def write_nodes(self, nodes_filename):
        """Write node patches to file"""
        
        nodes = open(nodes_filename, 'w')
        for np in self.node_patches:
            print >>nodes, np.node(), \
                            np.centroid(), \
                            np.normal(), \
                            np.alt_normal(), \
                            np.gc, \
                            np.planar_area(), \
                            np.bezier_area(), \
                            np.weighted_area(),\
                            np.num_quad_points(), \
                            np.energy_coefficient_f, \
                            np.energy_coefficient_h, \
                            np.force_coefficient_f, \
                            np.force_coefficient_h, \
                            np.f, \
                            np.h
        nodes.close()
        return
    
    def write_quads(self, quads_filename):
        """Write quadrature points to file"""
        
        quads = open(quads_filename, 'w')
        for np in self.node_patches:
            for ii in range(np.num_quad_points()):
                print >>quads, np.get_quad_point(ii)
        quads.close()
        return

    def write_energies(self, energies_filename):
        """Write the precalculated energy and force coefficients to file."""

        energies = open(energies_filename, 'w')
        for np in self.node_patches:
            fx = np.force_coefficient_f.x
            fy = np.force_coefficient_f.y
            fz = np.force_coefficient_f.z
            hx = np.force_coefficient_h.x
            hy = np.force_coefficient_h.y
            hz = np.force_coefficient_h.z
            print >>energies, "%2.16e %2.16e %2.16e %2.16e %2.16e %2.16e %2.16e %2.16e" %(np.energy_coefficient_f, np.energy_coefficient_h, fx, fy, fz, hx, hy, hz)

        energies.close()
        return

class MeshInstance(_MeshInstance, object):
    
    def __init__(self, *args):
        super(MeshInstance, self).__init__(*args)
    
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

    @property
    def mesh(self):
        return Mesh(self.get_mesh())
            
class BEEP(_BEEP, object):
    
    def __init__(self, *args):
        super(BEEP, self).__init__(*args)

    # generator function to provide node patch iterator
    @property
    def node_patches(self):
        for x in xrange(self.num_node_patches):
            yield self.get_patch(x)
            
    # generator function to provide mesh
    @property
    def mesh_instances(self):
        for x in xrange(self.num_mesh_instances):
            yield MeshInstance(self.get_mesh_instance(x))

    # generator function to provide mesh
    def get_mesh_inst(self, idx):
        return MeshInstance( self.get_mesh_instance(idx) )

    def kinemage(self, filename, fmax=None, hmax=None):
        """Create a kinemage representation of the current ensemble (with f/h values)"""
        import kintools
        from math import fabs

        print "Writing %s ... " %(filename),
        if fmax is None:
            fmax = max([max([fabs(np.f) for np in minst.node_patches]) for minst in self.mesh_instances])
            print "\nfmax: %f", fmax
        if hmax is None:
            hmax = max([max([fabs(np.h) for np in minst.node_patches]) for minst in self.mesh_instances])
            print "\nhmax: %f", hmax
        
        print fmax, hmax
        self.py_kinemage(filename, fmax, hmax, 100, kintools.hundred_red_blue_colours())
        
        print "done"
                
