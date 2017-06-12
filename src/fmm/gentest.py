#!/usr/bin/python

import random
import sys
import math

assert(len(sys.argv) >= 3)

beta = float(sys.argv[1])

def _rand_xyzc(edge_length=1.0):
    """Generate a random vector within box of given edge_length."""
    xyz = [(random.random() - 0.5)*edge_length for i in range(3)]
    xyz.append( 2 * (random.random() - 0.5 ))
    return xyz

for i in range(int(sys.argv[2])):
    pot = 0.0
    x,y,z,c = _rand_xyzc()
    #for ii in range(len(all)):
    #    if i == ii : continue
    #    xx,yy,zz,cc = all[ii]
    #    r =  math.sqrt( (xx-x)**2 + (yy-y)**2 + (zz-z)**2 )
    #    pot += cc * math.exp(-r*beta) / r
    
    print "%f %f %f %f" %(x,y,z,c)
