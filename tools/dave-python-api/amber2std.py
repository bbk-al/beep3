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
# script to convert amber style pdb/gro files to standard
# dseelig@gwdg.de

import sys,os
from pymacs import *
from pymacs import library

ids = ['A','B','C','D','E','F','G','H','I','J','K','L',\
       'M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def real_resname(r):
    dic = {'LYP':'LYS','LYSH':'LYS','CYM':'CYS',
           'CYS2':'CYS','CYN':'CYS','HIE':'HIS','HIP':'HIS',
           'HID':'HIS','HISA':'HIS','HISB':'HIS',
           'HISH':'HIS'}
    if dic.has_key(r): return dic[r]
    else: return r


desc=("This script reads a structure file",
      "with amber names and outputs",
      "it with standard names.",
      "It also adds chain identifyers.",
      "Works only for proteins....")

# define input/output files

files= [
    (efSTO, "-f", "amber.pdb", ffREAD ),
    (efPDB, "-o", "standard.pdb", ffWRITE ),
    ]

# define options

options=[]

# pass options, files and the command line to pymacs

args=parse_common_args(sys.argv,files,options, desc)



model = Model().read(args['-f'])

# find start and end residues of chains
starts = []
ends = []

for atom in model.atoms:
    if atom.name == 'H1' and atom.resname in library._one_letter.keys():
        starts.append(atom.resnr-1)
    elif atom.name in ['OC1','O1','OT1'] and atom.resname in library._one_letter.keys():
        ends.append(atom.resnr-1)

# add chain identifiers
for i, idx in enumerate(starts):
    resis = model.residues[starts[i]:ends[i]+1]
    for r in resis:
        r.set_chain_id(ids[i])

model.make_chains()
model.make_residues()

# rename ile
for r in model.residues:
    if r.resname == 'ILE':
        try:
            cd, hd1,hd2,hd3 = r.fetchm(['CD','HD1','HD2','HD3'])
            cd.name = 'CD1'
            hd1.name = 'HD11'
            hd2.name = 'HD12'
            hd3.name = 'HD13'
        except:
            pass
# rename terminal ox atoms
for chain in model.chains:
    r = chain.residues[-1]
    try:
        o1,o2 = r.fetchm(['OC1','OC2'])
        o1.name = 'O1'
        o2.name = 'O2'
    except:
        pass

for chain in model.chains:
    nterm = chain.residues[0]
    if nterm.resname in library._protein_residues:
        if len(nterm.resname) == 4:
            nterm.set_resname(nterm.resname[1:])
    cterm = chain.residues[-1]        
    if cterm.resname in library._protein_residues:
        if len(cterm.resname) == 4:
            cterm.set_resname(cterm.resname[1:])


for r in model.residues:
    r.set_resname(real_resname(r.resname))
    

print 'New model:'
print 'Model contains %d chains' % len(model.chains)
for ch in model.chains:
    print '> chain %s: nres = %d, natom = %d' % (ch.id, len(ch.residues), len(ch.atoms))
    print '%s' % ch.get_sequence()
    
print '\nWriting output file: %s' % args['-o']['fns']
model.write(args['-o'])

thanx()
