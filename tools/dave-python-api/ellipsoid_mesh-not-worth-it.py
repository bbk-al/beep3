#!/usr/bin/python

import gts_utils
from entities import Charge
from pqrtools import getCharges
from vector import *
from numpy import matrix, mat, linalg, eye, zeros, array, ndarray, dot, transpose, multiply
from pulp import *

def create_bounding_ellipsoid_mesh(pqr_filename, ellipsoid_filename):
    
    # get the charges
    charges = getCharges(pqr_filename)
    import random
    random.shuffle(charges)
    
    basis = []
    
    while (len(basis) < 4):
        basis.append(charges.pop())

    A, ellipsoid_centre = find_ellipse_by_LP(init)
    

    
    
    
    print len(charges)
    import sys
    print sys.getrecursionlimit()
    sys.setrecursionlimit(len(charges)*10)
    
    
    for c in charges:
        assert(in_ellipse(A,ellipsoid_centre,c.position))
        
    print ellipsoid_centre
        
    ## create moment of inertia tensor
    #inertia_tensor = zeros([3,3])
    #for c in charges:
        #x,y,z = c.position - cofg
        #inertia_tensor[0,0] += y*y + z*z
        #inertia_tensor[0,1] += -x*y
        #inertia_tensor[0,2] += -x*z
        #inertia_tensor[1,0] += -y*x
        #inertia_tensor[1,1] += x*x + z*z
        #inertia_tensor[1,2] += -y*z
        #inertia_tensor[2,0] += -z*x
        #inertia_tensor[2,1] += -z*y
        #inertia_tensor[2,2] += x*x + y*y

    # get eigenvectors of the inertia matrix -- they will 
    # correspond to axes of the molecule which have greatest
    # spatial extent
    eigvals,eigvecs = linalg.eig(A)

    # get the three principle axes
    v1 = Vector(list(eigvecs[:,0])).normal()
    v2 = Vector(list(eigvecs[:,1])).normal()
    v3 = v1.cross(v2)
    eigvecs[:,2] = v3

    # create rotation matrix from standard xyz to v1v2v3 basis
    rotation_matrix = (eye(3) * matrix(eigvecs).I).T

    # find maximum dimension along each axis
    v1_max, v2_max, v3_max = 0, 0, 0
    for c in charges:
        rel_xyz = c.position - ellipsoid_centre
        v1_max = max(v1_max, abs(rel_xyz.dot(v1))+c.radius)
        v2_max = max(v2_max, abs(rel_xyz.dot(v2))+c.radius)
        v3_max = max(v3_max, abs(rel_xyz.dot(v3))+c.radius)
        
    vertices, triangles = gts_utils.get_spherical_mesh(3)
    vnormals = []
    for v in vertices:
        
        # rotate vertex into new coordinate frame
        newx = (rotation_matrix[0,0] * v.x) + \
               (rotation_matrix[0,1] * v.y) + \
               (rotation_matrix[0,2] * v.z)
        newy = (rotation_matrix[1,0] * v.x) + \
               (rotation_matrix[1,1] * v.y) + \
               (rotation_matrix[1,2] * v.z)
        newz = (rotation_matrix[2,0] * v.x) + \
               (rotation_matrix[2,1] * v.y) + \
               (rotation_matrix[2,2] * v.z)
        v.x = newx
        v.y = newy
        v.z = newz
        
        vnormals.append(v)
        
        # scale axes to form ellipsoid
        v.x *= v1_max
        v.y *= v2_max
        v.z *= v3_max
        
        # move the centre of the ellipsoid to centre of mass
        # of the molecule
        v.x += ellipsoid_centre.x
        v.y += ellipsoid_centre.y
        v.z += ellipsoid_centre.z
    
    gts_utils.write_gts(ellipsoid_filename, vertices, triangles, vertex_normals=vnormals)
    
    import kintools
    print kintools.triangleList(triangles, vertices)
    print "@balllist {centre}\nr=1.0 {color=blue} %f %f %f" \
          %(ellipsoid_centre.x, ellipsoid_centre.y, ellipsoid_centre.z)
    print "@balllist {charges}"
    for c in charges:
        print c.kinemage(Vector(0,0,0),None,Vector(0,0,0))

def smallest_ellipse(interior_Q, support_R):
    
    if (len(interior_Q) == 0 or len(support_R) == 9):
        return find_ellipse_by_LP(support_R[:])
    else:
        q = interior_Q.pop(0)
        A,c = smallest_ellipse(interior_Q[:], support_R[:])
        if (in_ellipse(A,c,q)):
            return A,c

    support_R.append(q)
    return smallest_ellipse(interior_Q[:], support_R[:])
        
def in_ellipse(A,c,p):
    import numpy
    
    if A.all() == 0: return False
    pp = array(p)
    cc = array(c)
    pc = pp - cc
    n = dot(transpose(pc),dot(A, pc))
    return n <= 1.0
        
def find_ellipse_by_LP(constraint_points_raw):

    A = zeros([3,3])
    cx,cy,cz = 0,0,0
    
    if (len(constraint_points_raw) < 7):
        return A,Vector(cx,cy,cz)
    
    constraint_points = constraint_points_raw[:]
    while len(constraint_points) < 9:
        constraint_points.append(constraint_points[-1])
    
    assert(len(constraint_points) <= 9)
    
    # set up problem
    prob = LpProblem("Ellipse", LpMinimize)
    
    a11 = LpVariable("a11", cat='Continuous')
    a12 = LpVariable("a12", cat='Continuous')
    a13 = LpVariable("a13", cat='Continuous')
    a22 = LpVariable("a22", cat='Continuous')
    a23 = LpVariable("a23", cat='Continuous')
    a33 = LpVariable("a33", cat='Continuous')
    mx = LpVariable("mx", cat='Continuous')
    my = LpVariable("my", cat='Continuous')
    mz = LpVariable("mz", cat='Continuous')
    w = LpVariable("w", lowBound=10, cat='Continuous')

    # don't care about objective -- only one compatible solution
    prob += 0, "dummy objective"
    #prob += w > 0, "dummy constraint on w"
    
    for p in constraint_points:
        
        prob += (p.x)*(a11*p.x + a12*p.y + a13*p.z) + \
                (p.y)*(a12*p.x + a22*p.y + a23*p.z) + \
                (p.z)*(a13*p.x + a23*p.y + a33*p.z) \
                 + 2.0*(mx*p.x + my*p.y + mz*p.z) + w == 0, ""
    
    prob.solve()
    #print "Status:", LpStatus[prob.status]
    #for v in prob.variables():
        #print v.name, "=", v.varValue

    A[0,0] = a11.varValue
    A[0,1] = a12.varValue
    A[0,2] = a13.varValue
    A[1,0] = a12.varValue
    A[1,1] = a22.varValue
    A[1,2] = a23.varValue
    A[2,0] = a13.varValue
    A[2,1] = a23.varValue
    A[2,2] = a33.varValue
    
    invA = linalg.inv(A)
    m = array([mx.varValue, my.varValue, mz.varValue])
    c = transpose(m) * invA
    cx = invA[0,0]*mx.varValue + invA[1,0]*my.varValue + invA[2,0]*mz.varValue
    cy = invA[0,1]*mx.varValue + invA[1,1]*my.varValue + invA[2,1]*mz.varValue
    cz = invA[0,2]*mx.varValue + invA[1,2]*my.varValue + invA[2,2]*mz.varValue
    c = array([cx, cy, cz])
    
    p = constraint_points[0]
    pc = array([p.x-cx, p.y-cy, p.z-cz])
    n = dot(transpose(pc),dot(A, pc))
    A /= n

    return A, Vector(cx,cy,cz)
        
if (__name__ == "__main__"):
    
    import sys
    create_bounding_ellipsoid_mesh(sys.argv[1], sys.argv[2])
    
