#!/usr/bin/env python
#############################################################
#                                                           #
#           THIS FILE IS PART OF PYMACS                     #
#   WRITTEN BY DANIEL SEELIGER (dseelig@gwdg.de)            #
#                   2006-2009                               #
#    COMPUTATIONAL BIOMOLECULAR DYNAMICS GROUP              #
#  MAX-PLANCK-INSTITUTE FOR BIOPHYSICAL CHEMISTRY           #
#                   GOETTINGEN                              #
#                                                           #
#############################################################
# pymacs script to make pdb files ready for use with amber
# dseelig@gwdg.de

import sys,os
from pymacs import *

pdb_format="%6s%5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n"

# describe what the script does

desc=("This script reads a structure file",
      "and outputs a new structure file",
      "with some renamed residues.",
      "The N-terminal and C-terminal residues",
      "are renamed such that we can use the",
      "AMBER force field. We also check for",
      "cysteine residues. So far the script works for"
      "proteins only. \n"
      "Take an all-atom protonated input pdb file (pdb2gmx -ff oplsaa)."
      "Renaming of residues (e.g. HIS) will be done based on"
      "the protons in the structure file.")

# define input/output files

files= [
    (efPDB, "-f", "protein.pdb", ffREAD ),
    (efPDB, "-o", "amber.pdb", ffWRITE ),
    ]

# define options

options=[]

# pass options, files and the command line to pymacs

args=parse_common_args(sys.argv,files,options, desc)


# read structure file
model = Model().read(args['-f'])


# now check cysteines
print '\nChecking cys....'
cysl = model.fetch_residues('CYS')     # get a list with all cysteines

# we do a simple check. If a HG is there it's CYS, else it's CYS2

for res in cysl:
    hg = res.fetch_atoms('HG')        # select HG atom from residue
    if not hg:    # found nothing
        print 'Residue %d-%s (chain %s) will become %s' % (res.id, res.resname,res.chain_id,'CYS2')
        res.set_resname('CYS2')
    else:
        res.set_resname('CYN')
        print 'Residue %d-%s (chain %s) will become %s' % (res.id, res.resname, res.chain_id, res.resname)


# lysine
print 'Checking lys....'
lysl = model.fetch_residues('LYS')
for res in lysl:
    at = res.fetch('HZ3')
    if at:
        res.set_resname('LYP')
        

# histidine
print 'Checking his......'
hisl = model.fetch_residues('HIS')
for res in hisl:
    bHE2 = False
    bHD1 = False
    he2 = res.fetch('HE2')
    if he2: bHE2 = True
    hd1 = res.fetch('HD1')
    if hd1: bHD1 = True
    if hd1 and he2:
        res.set_resname('HIP')
    elif hd1 and not he2:
        res.set_resname('HID')
    elif he2 and not hd1:
        res.set_resname('HIE')
    print 'Residue %d-%s (chain %s) will become %s' % (res.id, 'HIS', res.chain_id, res.resname)
    
        

print 'Checking termini.....'
for chain in model.chains:
    print 'Processing chain %s' % chain.id
    first = chain.residues[0]      # first residue
    last = chain.residues[-1]      # last residue
    if first.resname in library._one_letter.keys() and \
           last.resname in library._one_letter.keys():  # its a protein chain
        first.set_resname('N'+first.resname) # rename e.g. ALA to NALA
        last.set_resname('C'+last.resname)   # rename e.g. ARG to CARG
        try:
            o1,o2 = last.fetchm(['O1','O2'])
            o1.name = 'OC1'
            o2.name = 'OC2'
        except:
            print 'Error: No terminal O1, O2 atoms found in chain %s' % chain.id
            print '       In pdb2gmx generated structures these should be there.'
            print '       Exiting'
            sys.exit(1)

# hack to get pdb file with 4 character residue names
for atom in model.atoms:
    atom.chain_id = ' '
    if len(atom.resname)==4:
        atom.chain_id = atom.resname[-1]+atom.chain_id
    
fp = open(args['-o']['fns'],'w')
print 'Writing ouput to: %s' % args['-o']['fns'] 
for atom in model.atoms:
    print >>fp, atom

thanx()
