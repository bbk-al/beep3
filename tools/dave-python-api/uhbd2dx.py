#!/usr/bin/python

from _BEM import GridParms
import opendx_utils
import math
from collections import deque

# UHBD values are in kCal/mol.e -- conv_factor is kT/e in those units (APBS
# format is dimensionless in units of kT/e, so divide UHBD vals by this
# constant to get consistency with APBS)
conv_factor = 0.2388*6.02e23*1.38e-23*293.0*1e-3

def read_uhbd(uhbd_grd_filename, nx, ny, nz, h, ox, oy, oz):

    grid = GridParms()
    grid.x_pts = nx
    grid.y_pts = ny
    grid.z_pts = nz
    grid.x_angstroms = (nx-1) * h
    grid.y_angstroms = (ny-1) * h
    grid.z_angstroms = (nz-1) * h
    grid.origin_x = ox 
    grid.origin_y = oy
    grid.origin_z = oz

    from numpy import zeros
    grd_data = zeros((grid.x_pts,grid.y_pts,grid.z_pts),'d')
    
    f = open(uhbd_grd_filename, 'r')

    for preamble in range(5):
        f.readline()
    
    for kk in range(grid.z_pts):
        
        # sanity check
        chk_line = [int(val) for val in f.readline().split()]
        assert(len(chk_line) == 3)
        kn, itot, jtot = chk_line
        assert(kn == kk+1 and itot == grid.x_pts and jtot == grid.y_pts)
        
        # six numbers per line
        expected_num_lines = int(math.ceil(itot*jtot / 6.0))

        # put data into a fifo (deque is double ended list)
        buf = deque()
        for line_ctr in range(expected_num_lines):
            buf.extend([float(val) for val in f.readline().split()])

        # check we got as much data as we were expecting
        assert(len(buf) == itot*jtot)
            
        # popleft from deque to use as a fifo
        for jj in range(jtot):
            for ii in range(itot):
                grd_data[ii,jj,kk] = buf.popleft() / conv_factor
                
        # sanity check
        assert(len(buf) == 0)

    f.close()
    
    return grid, grd_data

if __name__ == "__main__":
    
    import sys
    filename = sys.argv[1]
    nx = int(sys.argv[2])
    ny = int(sys.argv[3])
    nz = int(sys.argv[4])
    h = float(sys.argv[5])
    ox = float(sys.argv[6])
    oy = float(sys.argv[7])
    oz = float(sys.argv[8])
    
    grid, data = read_uhbd(filename, nx, ny, nz, h, ox, oy, oz)
    opendx_utils.write_opendx_file(filename+".dx", grid, data)
    
    

    
