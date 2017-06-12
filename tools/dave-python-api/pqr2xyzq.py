#!/usr/bin/python

import pqrtools
import sys

pqr = sys.argv[1]
xyzq = sys.argv[2]
# might be a scale factor passed in to rescale q values
try:
    scale_factor = sys.argv[3]
except IndexError:
    scale_factor = 1.0

charge_list = pqrtools.getCharges(pqr)

f = open(xyzq, 'w')
for c in charge_list:
    x,y,z = c.position
    q = c.charge * scale_factor
    r = c.radius
    print >>f, "%f %f %f %f %f" %(x,y,z,q,r)
f.close()
