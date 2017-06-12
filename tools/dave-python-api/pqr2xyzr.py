#!/usr/bin/python

import pqrtools
import sys

pqr = sys.argv[1]
xyzr = sys.argv[2]
charge_list = pqrtools.getCharges(pqr)

f = open(xyzr, 'w')
for c in charge_list:
    x,y,z = c.position
    r = c.radius
    if (c.charge != 0.0 and r == 0.0):
        print "WARNING: Point charge (%f) defined at (%f,%f,%f) -- setting radius to 1.0" %(c.charge, x,y,z)
        r = 1.0

    print >>f, "%f %f %f %f" %(x,y,z,r)
f.close()
