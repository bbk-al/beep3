#!/usr/bin/python

import gts_utils
from entities import Charge
from pqrtools import getCharges
from vector import *
from numpy import linalg, eye, zeros, array, ndarray, dot
from _ELLIPSOID import solve_ellipsoid
import sys
from string import atoi

class ellipsoid(object):
    
    def __init__(self, centre, v1, v2, v3, radii):

        self.centre = centre
        self.semi_axis1 = v1
        self.semi_axis2 = v2
        self.semi_axis3 = v3
        self.r1, self.r2, self.r3 = radii

        self.max_radius = max(self.r1, self.r2, self.r3)

    def generate_vertices(self, detail):
        self.vertices, self.triangles = gts_utils.get_spherical_mesh(detail)
        self.vnormals = []
        for v in self.vertices:
            v1 = self.semi_axis1
            v2 = self.semi_axis2
            v3 = self.semi_axis3
            vnew = v1*v.dot(v1)*self.r1 + v2*v.dot(v2)*self.r2 + v3*v.dot(v3)*self.r3
            v.x, v.y, v.z = vnew.x, vnew.y, vnew.z
            
            self.vnormals.append(Vector(v.x,v.y,v.z))
            
            # move the centre of the ellipsoid
            v.x += self.centre.x
            v.y += self.centre.y
            v.z += self.centre.z
    
    def write_gts(self, ellipsoid_filename):
        gts_utils.write_gts(ellipsoid_filename, self.vertices, self.triangles, vertex_normals=self.vnormals)
        
    def pt_in_native_frame(self, pt):
        v = pt - self.centre
        return Vector(v.dot(self.semi_axis1), 
                      v.dot(self.semi_axis2),
                      v.dot(self.semi_axis3))

    def pt_intersection_test(self, pt):
        """Returns if the point pt is inside the ellipsoid."""

        if ((pt - self.centre).length() > self.max_radius):
            return False
        
        # rotate the point into the native coordinate frame of the ellipsoid
        ex, ey, ez = self.pt_in_native_frame(pt)
        
        # see if the scaled vector components add up to less than unity
        if ( ((ex*ex)/(self.r1*self.r1) + 
              (ey*ey)/(self.r2*self.r2) + 
              (ez*ez)/(self.r3*self.r3)) < 1.0):
            return True
        else:
            return False
        
    def kinemage(self, handle=sys.stdout):
        
        import kintools
        print >>handle, kintools.triangleList(self.triangles, self.vertices)
        print >>handle, "@vectorlist {principle axes}"
        radii = Vector(self.r1,self.r2,self.r3)
        for ii,axis_dir in enumerate([self.semi_axis1,self.semi_axis2,self.semi_axis3]):
            axis = self.centre + (axis_dir * radii[ii])
            print "{} X %f %f %f %f %f %f" %(self.centre.x, self.centre.y, self.centre.z, axis.x, axis.y, axis.z)
        print >>handle, "@balllist {charges}"
        for c in self.charges:
            print >>handle, c.kinemage(Vector(0,0,0),None,Vector(0,0,0))

    def dummy_atoms(self):
        """Generates a set of points which are within the ellipse."""
        
        import math
        print "@balllist {dummies}"
        for xctr in range(-int(math.ceil(self.r1)),int(math.ceil(self.r1))):
            for yctr in range(-int(math.ceil(self.r2)),int(math.ceil(self.r2))):
                for zctr in range(-int(math.ceil(self.r3)),int(math.ceil(self.r3))):
                    pt = self.centre + self.semi_axis1*(xctr*1.0) \
                                     + self.semi_axis2*(yctr*1.0) \
                                     + self.semi_axis3*(zctr*1.0)
                    if (self.pt_intersection_test(pt)):
                        #print "r=1.0 {} %f %f %f" %(pt.x, pt.y, pt.z)
                        print "ATOM 0 XXX XXX 0 %f %f %f 0.0 1.0" %(pt.x, pt.y, pt.z)
        
            
def create_bounding_ellipsoid_mesh(pqr_filename):
    
    # get the charges
    charges = getCharges(pqr_filename)
    
    points = [c.position for c in charges]
    centre = Vector(0,0,0)
    v1 = Vector(0,0,0)
    v2 = Vector(0,0,0)
    v3 = Vector(0,0,0)
    radii = Vector(0,0,0)

    # call CGAL to do the heavy lifting
    solve_ellipsoid(points, centre, v1, v2, v3, radii)
    
    # expand by maximum radius of atoms to ensure enclosure
    max_atom_radius = 2.0
    radii += Vector(max_atom_radius, max_atom_radius, max_atom_radius)

    print >>sys.stderr, "%f %f %f" %(centre.x, centre.y, centre.z)
    print >>sys.stderr, "%f %f %f %f" %(v1.x, v1.y, v1.z, radii.x)
    print >>sys.stderr, "%f %f %f %f" %(v2.x, v2.y, v2.z, radii.y)
    print >>sys.stderr, "%f %f %f %f" %(v3.x, v3.y, v3.z, radii.z)

    #print "Semi-axis lengths:", radii
    e = ellipsoid(centre, v1, v2, v3, radii)
    e.charges = charges[:]
    
    return e
    

if (__name__ == "__main__"):
    
    ellipsoid = create_bounding_ellipsoid_mesh(sys.argv[1])
    ellipsoid.generate_vertices(atoi(sys.argv[3]))
    ellipsoid.write_gts(sys.argv[2])
    ellipsoid.kinemage()
    #ellipsoid.dummy_atoms()
