#!/usr/bin/python

#
# pdbtools.py
#
# Some useful PDB functions (wraps other existing PDB handlers).
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

from Scientific.IO.PDB import PDBFile
from entities import Atom

class SpecialPDBFile(PDBFile):

    def __init__(self, file, subformat = None):
        """
        @param filename: the name of the PDB file
        @type filename: C{str}
        @param mode: the file access mode, 'r' (read) or 'w' (write)
        @type mode: C{str}
        @param subformat: indicates a specific dialect of the PDB format.
                          Subformats are defined in
                          L{Scientific.IO.PDBExportFilters}; they are used
                          only when writing.
        @type subformat: C{str} or C{NoneType}
        """
        self.file = file
        self.output = True
        self.export_filter = None
        if subformat is not None:
            export = export_filters.get(subformat, None)
            if export is not None:
                self.export_filter = export()
        self.open = 1
        if self.output:
            self.data = {'serial_number': 0,
                         'residue_number': 0,
                         'chain_id': '',
                         'segment_id': ''}
            self.het_flag = 0
            self.chain_number = -1

def to_pdb_format(words):
    """Converts a list of PDB line parts into a formatted PDB line."""

    # there should be 10 parts to the words list
    assert(len(words) >= 10)

    output = "%-6s%5i %-4s %3s %1s %3i    %8.3f%8.3f%8.3f%6.2f%6.2f" % (
        words[0],        # "ATOM"
        int(words[1]),   # atom number
        words[2],        # element
        words[3],        # amino acid
        ' ',             # chain
        int(words[4]),   # residue number
        float(words[5]), # x
        float(words[6]), # y
        float(words[7]), # z
        float(words[8]), # 
        float(words[9])  # 
        )

    return output

def get_atoms_from_pdb(pdb_filename):
    
    from Scientific.IO.PDB import Structure
    struct = Structure(pdb_filename)
    atoms = get_atoms_from_pdb_structure(struct)
    
    return atoms
    
def get_atoms_from_pdb_structure(struct):
    """Allow re-cycling of structure object."""
    
    atoms = []
    for r in struct.residues:
        for atom in r:
            element = atom.properties['element']
            atoms.append(Atom(atom.position, name=element))

    return atoms

if __name__ == "__main__":

    import sys

    try:
        input_file = sys.argv[1]
    except IndexError:
        print "Please specify a (badly formed) PDB file."
        sys.exit(1)

    try:
        for line in open(input_file):
            words = line.split()
            print to_pdb_format(words)

    except:
        print "Failed."
        sys.exit(1)

    # if we get here then everything's ok.
    sys.exit(0)

