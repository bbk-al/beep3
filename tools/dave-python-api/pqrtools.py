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
from _BEM import Vector
from Scientific.IO.PDB import PDBFile
from Scientific.IO.FortranFormat import FortranFormat
from entities import Charge

class PQRFile(PDBFile):

    atom_format = FortranFormat('A6,I5,1X,A4,A1,A4,A1,I4,A1,3X,3F8.2,2F6.2,' +
                            '6X,A4,2A2')

    def readLine(self):
        """
        Override of Scientific.IO.PDB.PDBFile.readLine function.

        @returns: the contents of one PQR record
        @rtype: C{tuple}
        """
        while True:
            line = self.file.readline()
            if not line: return ('END','')
            if line[-1] == '\n': line = line[:-1]
            line = string.strip(line)
            if line: break
        line = string.ljust(line, 80)
        type = string.strip(line[:6])
        if type == 'ATOM' or type == 'HETATM':
            line = [l.strip() for l in line.split()]
            #line = FortranLine(line, PQRFile.atom_format)
            data = {'serial_number': line[1],
                    'name': line[2],
                    'alternate': string.strip(line[3]),
                    'residue_name': string.strip(line[4]),
                    'position': Vector([string.atof(xyz) for xyz in line[5:8]]),
                    'charge': string.strip(line[8]),
                    'radius': string.strip(line[9])}
            return type, data
        else:
            return type, line[6:]

    
def getCharges(filename):
    """Returns a list of charges and corresponding lines from PQR file."""
    
    pqr = PQRFile(filename)
    pqr_data = []
    while True:
        type, data = pqr.readLine()
        if type == 'END': break
        if type == 'ATOM' or type == 'HETATM':
            num = string.atoi(data['serial_number'])
            xyz = data['position']
            charge = string.atof(data['charge'])
            radius = string.atof(data['radius'])
            pqr_data.append(Charge(xyz, charge, radius))
    
    return pqr_data

