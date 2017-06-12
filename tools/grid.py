#!/usr/bin/env python3
from vector import *
from geometry import _rand_rot, _rand_xyz

n = 4
radius = 42.0
diameter = radius * 2
origin = Vector(0,0,0)
edge_length = (n-1) * diameter
bottom_left_corner = origin - Vector(edge_length, edge_length, edge_length) / 2
print("edge length:", edge_length + 2*diameter)
print("bottom left:", bottom_left_corner)
#print "@kinemage"
#print "@dotlist"
ctr=0
for xx in range(n):
    xoff = Vector(diameter*xx,0,0)
    for yy in range(n):
        yoff = Vector(0,diameter*yy,0)
        for zz in range(n):
            zoff = Vector(0,0,diameter*zz)
            
            pt = bottom_left_corner + xoff + yoff + zoff
            rot = Quaternion(0,0,0,0)

            print("<instance dielectric=2.0 instance_id=%d mesh_id=0 x=%f y=%f z=%f a=%f b=%f c=%f d=%f/>" %(ctr, pt.x, pt.y, pt.z, rot.x, rot.y, rot.z, rot.w))
            #print "%f %f %f" %(pt.x, pt.y, pt.z)
            ctr += 1
print("")
