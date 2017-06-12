#!/usr/bin/python
# vim: set fileencoding=UTF8 :
#
# geometry.py
#
# Some useful geometry functions.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

import math
import random
from pybeep import Vector, Quaternion

def normalised(vector):
    """Return normalised (length=1.0) Vector version of the input.

    This function is here to remove some ambiguity in the term 'normal' where
    I sometimes mean perpendicular to a plane, and sometimes I mean,
    normalised to length=1.0 (i.e. a unit vector)."""

    return Vector(vector).normal()

def apply_quaternion_to_vector(quaternion, vector):
    """Return a vector rotated by a quaternion, as a Vector."""
    vnew = Vector(vector)
    vnew.apply_rotation(quaternion)
    return vnew

def local_to_universal(vector, 
                       ref_pt_local, 
                       ref_pt_universe,
                       rotation=None):
    """Convert vector in local co-ordinates to Universe co-ords.

    ref_pt_local and ref_pt_universe are the SAME point in
    each coordinate frame; furthermore the reference point is the centre of
    rotation of the local coordinate frame, so this is the single point in the
    local coordinate frame which is not changed by rotations."""

    # reference points are mandatory!
    assert(ref_pt_local is not None)
    assert(ref_pt_universe is not None)
    
    # get the vector w.r.t the stationary point (common point of reference
    # between local/universe coordinate frames, which is also centre of
    # rotation).
    new_vector = vector - ref_pt_local
    
    # apply rotation (if passed in) (rotation is the rotation of the local
    # coordinate frame w.r.t universe frame)
    if rotation is not None:
        rotated = apply_quaternion_to_vector(rotation, new_vector)
    else:
        rotated = new_vector # no rotation applied
        
    # we know the coordinates of the reference point in universe coords (it
    # gets passed in to this function as reference_point_universe); the
    # universe coordinates of the vector passed in is the universe-frame
    # offset from this point.
    result = rotated + ref_pt_universe
        
    return result

# These are some useful functions derived from Numerical Recipes §21.6
def point_inside_unit_circle():
    """Generate a random point within unit circle.
    
    Pick two randomly uniformly distributed numbers in range [-1,1]; reject
    any where x**2 + y**2 > 1.0.
    
    NB: The random.uniform function returns values up to, but not including,
    the upper boundary, so this is strictly slightly imperfect. Good enough
    for my purposes though."""
    
    while True:
        x = random.uniform(-1.0,1.0)
        y = random.uniform(-1.0,1.0)
        if (x*x + y*y) <= 1.0: break
    
    return (x,y)

def random_point_on_sphere(dimensions=3):
    """Return random point on sphere of unit radius in given number of dimensions."""

    if dimensions==3:

        u0,u1 = point_inside_unit_circle()
        useful = math.sqrt(1.0 - u0*u0 - u1*u1)
        
        x = 2.0*u0*useful
        y = 2.0*u1*useful
        z = 1.0 - 2.0*(u0*u0 + u1*u1)
        
        return (x,y,z)
    
    elif dimensions==4:
        
        u0,u1 = point_inside_unit_circle()
        u2,u3 = point_inside_unit_circle()
        
        x0 = u0
        x1 = u1
        useful = math.sqrt((1.0 - u0*u0 - u1*u1)/(u2*u2 + u3*u3))
        x2 = u2*useful
        x3 = u3*useful        
        
        return (x0,x1,x2,x3)
    
    else:
        # I haven't implemented routines for other dimensions.
        raise ValueError

#def random_rotation_matrix():
#    """Return a random rotation matrix in 3 dimensions."""
#    
#    from Scientific.Geometry import Tensor
#    x0,x1,x2,x3 = random_point_on_sphere(dimensions=4)
#    
#    rot = Tensor([[1.0-2.0*(x1*x1 + x2*x2), 2.0*(x0*x1 - x3*x2), 2.0*(x0*x2 + x3*x1)],
#                 [2.0*(x0*x1 + x3*x2), 1.0-2.0*(x0*x0 + x2*x2), 2.0*(x1*x2 - x3*x0)],
#                 [2.0*(x0*x2 - x3*x1), 2.0*(x1*x2 + x3*x0), 1.0-2.0*(x0*x0 + x1*x1)]])
#    
#    return rot

def cosine_rule(a,b,c):
    """Solves angle A given sides a,b,c (a is opposite side to angle A)."""
    return math.acos((-a*a + b*b + c*c) / (2.0*b*c))

import random
def _rand_xyz(edge_length=100):
    """Generate a random vector within box of given edge_length."""
    xyz = [(random.random() - 0.5)*edge_length for i in range(3)]
    return Vector(xyz[0], xyz[1], xyz[2])

def _rand_rot():
    """Generate a random rotation"""
#    from geometry import random_rotation_matrix as rand_rot_matrix
#    from Scientific.Geometry.Transformation import Rotation 
#    return Rotation(rand_rot_matrix()).asQuaternion()
    x,y,z,w = random_point_on_sphere(dimensions=4)
    return Quaternion(x,y,z,w)

def _rand_one_minus_one():
    
    import math
    n = math.floor(random.random() * 2.0)
    if n == 0: return -1
    else: return +1

def test_rand_rot():
    
    xx = Vector(1.0,0.0,0.0)
    from geometry import apply_quaternion_to_vector
    
    f = open("random_rotation_test.kin", "w")
    print >>f, "@kinemage"
    print >>f, "@dotlist"
    
    for ii in range(10000):
        rand_pt_on_surface_of_sphere = apply_quaternion_to_vector(_rand_rot(), xx)
        print >>f, "%f %f %f" %(rand_pt_on_surface_of_sphere.x,rand_pt_on_surface_of_sphere.y,rand_pt_on_surface_of_sphere.z)

    f.close()

if __name__ == "__main__":

    # TODO: put some test cases here.
    test_rand_rot()
    pass
