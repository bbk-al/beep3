#!/usr/bin/env python3
#from pybeep import Vector
from libBEEP import _Vector, Quaternion

# Simple addition
v1 = _Vector(1,2,3)
v2 = _Vector(4,5,6)
v = v1 + v2
print(v, "= (5,7,9)")

# Null rotation and translation
q = Quaternion(1,0,0,0)
v1 = _Vector(0,0,0)
v2 = _Vector(0,0,0)
v.change_coordinate_frame(v1, q, v2)
print(v, "= (5,7,9)")

# Rotation
q = Quaternion(0,1,0,0)
v1 = _Vector(0,0,0)
v2 = _Vector(0,0,0)
v.change_coordinate_frame(v1, q, v2)
print(v, "= (5,-7,-9)")

# Translation
q = Quaternion(1,0,0,0)
v1 = _Vector(0,0,0)
v2 = _Vector(1,2,3)
v.change_coordinate_frame(v1, q, v2)
print(v, "= (6,-5,-6)")

