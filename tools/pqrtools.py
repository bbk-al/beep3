#!/usr/bin/python
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
#
# pqrtools.py
#
# A PDBFile-like object to represent PQR files (like PDB format, but with
# charge and radius instead of occupancy and beta factor). Used for electrostatics
# stuff, e.g. APBS, and my BEM implementation.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#
import string
from pybeep import Vector, Charge
from pdb import *

def getCharges(filename):
    """Returns a list of charges and corresponding lines from PQR file."""

    pqr_file = open(filename,'r')
    pqr_data = []
    pqr_atoms, err_lines = readPQR(pqr_file)

    for pqr in pqr_atoms:
        xyz = Vector(pqr.x, pqr.y, pqr.z)
        pqr_data.append(Charge(xyz, pqr.q, pqr.r))
    return pqr_data


def getResidueCharges(filename):
    """Returns a list of residues and total charge located at Calpha from PQR file."""

    pqr_file = open(filename,'r')
    pqr_data = []
    pqr_atoms, err_lines = readPQR(pqr_file)

    atom_set = []
    residue = []
    last_seq = -1
    for pqr in pqr_atoms:
        if (pqr.resSeq == 0): raise Exception # can't identify residues by sequence number
        if (pqr.resSeq != last_seq):
            if (last_seq != -1):
                atom_set.append(residue)
            residue = [pqr]
            last_seq = pqr.resSeq
        else:
            residue.append(pqr)
    if (last_seq != -1 and len(residue) > 0):
        atom_set.append(residue)
        
    for residue in atom_set:
        
        backbone =  [xx for xx in residue if xx.name in ["O","N","C"] ]
        assert(len(backbone) == 3)
        sidechain = [xx for xx in residue if xx.name not in ["CA","O","N","C"] ]
        alpha = [pqr for pqr in residue if pqr.name == "CA"]
        assert(len(alpha)==1)
        
        for pqr in backbone:
            pqr_data.append(Charge(Vector(pqr.x,pqr.y,pqr.z), pqr.q, pqr.r))

        calpha_xyz = Vector(alpha[0].x,alpha[0].y,alpha[0].z)
        total_sidechain_charge = sum([pqr.q for pqr in sidechain])
        pqr_data.append(Charge(calpha_xyz, total_sidechain_charge + alpha[0].q, alpha[0].r))

    return pqr_data

def pqr2calpha(pqr):
    """Convert pqr to xyzr file."""
    
    charge_list = getResidueCharges(pqr)
    
    for c in charge_list:
        xyz = c.position()
        q = c.charge
        r = c.radius
        print "%f %f %f %f %f" %(xyz.x, xyz.y, xyz.z, q, r)

def pqr2xyzr(pqr):
    """Convert pqr to xyzr file."""
    
    charge_list = getCharges(pqr)
    
    for c in charge_list:
        xyz = c.position()
        r = c.radius
        print "%f %f %f %f" %(xyz.x, xyz.y, xyz.z, r)

def pqr2xyzqr(pqr, filename=None, scale_factor=1.0):
    """Convert pqr to xyzqr file."""
    
    charge_list = getCharges(pqr)
    if filename is not None:
        f = open(filename,'w')
        for c in charge_list:
            xyz = c.position()
            q = c.charge * scale_factor
            r = c.radius
            print >>f, "%f %f %f %f %f" %(xyz.x, xyz.y, xyz.z, q,r)
        f.close()
    else:
        for c in charge_list:
            xyz = c.position()
            q = c.charge * scale_factor
            r = c.radius
            print "%f %f %f %f %f" %(xyz.x, xyz.y, xyz.z, q,r)
    return


if __name__ == "__main__":
    import sys
    
    mode = sys.argv[1]
    if (locals().has_key(mode)):
        locals()[mode](*sys.argv[2:])
    else:
        print "Can't run that: try one of these... ", locals()  
       
