#!/usr/bin/python

import pdbtools
import sys

pdb = sys.argv[1]
xyzr = sys.argv[2]
atom_list = pdbtools.get_atoms_from_pdb(pdb)

f = open(xyzr, 'w')
for at in atom_list:
    x,y,z = at.position
    r = at.radius
    print >>f, "%f %f %f %f" %(x,y,z,r)
f.close()
