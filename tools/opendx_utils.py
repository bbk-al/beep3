# -*- coding: utf-8 -*-
from libBEEP import GridParms
import math

# exception class used for opendx grids
class BadGridPoint(Exception):
    def __init__(self, *args):
        Exception.__init__(self, *args)

class PotentialGrid(object):
        
    def __init__(self, filename):
        
        self.dx_filename = filename
        self.grid, self.matrix = PotentialGrid.read_opendx_file(filename)
        
    def get_potential(self, xyz):
        """Get the grid potential at an xyz location"""

        try:
            x,y,z = self.real_xyz_to_grid_idx(self.grid, xyz.x, xyz.y, xyz.z)
        
        except BadGridPoint:
            #print "position passed in (%s) is outside of the grid defined by this opendx file" %(xyz)
            return 0.0
        
        return self.matrix[x,y,z]
    
    @staticmethod
    def copy_GridParms_from_opendx_file(filename):
        """Copy the grid layout from an existing OpenDX file."""

        # open the openDX file
        opendx_file = open(filename,'r')

        # get the preamble lines
        preamble = PotentialGrid.get_preamble(opendx_file)
        
        # close the dx file
        opendx_file.close()

        #
        # extract the gridparms
        #
        
        # extract the number of grid points in xyz directions
        gridx, gridy, gridz = [int(griddim) 
                            for griddim in preamble[0].split()[-3:]]
        # bottom left corner
        or_x, or_y, or_z    = [float(griddim) 
                            for griddim in preamble[1].split()[-3:]]
        delta_x = float(preamble[2].split()[1])
        delta_y = float(preamble[3].split()[2])
        delta_z = float(preamble[4].split()[3])

        grid = GridParms()
        grid.x_pts = gridx
        grid.y_pts = gridy
        grid.z_pts = gridz
        grid.x_angstroms = (gridx-1) * delta_x
        grid.y_angstroms = (gridy-1) * delta_y
        grid.z_angstroms = (gridz-1) * delta_z
        grid.origin_x = or_x
        grid.origin_y = or_y
        grid.origin_z = or_z
        
        return grid
        
    @staticmethod
    def get_preamble(opendx_file_handle):
        
        preamble = []

        while len(preamble) < 7:
            line = opendx_file_handle.readline()
            if line[0] != '#':
                preamble.append(line.strip())

        return preamble
    
    @staticmethod
    def read_opendx_file(opendx_filename):
        
        grid = PotentialGrid.copy_GridParms_from_opendx_file(opendx_filename)
        
        from numpy import zeros
        results = zeros((grid.x_pts,grid.y_pts,grid.z_pts),'d')
        num_lines = math.ceil(grid.num_grid_points() / 3.0)
        
        opendx = open(opendx_filename)
        PotentialGrid.get_preamble(opendx) # skip preamble lines

        ctr = 0
        from collections import deque
        all_results = deque()

        # helper function
        def refill(opendx, all_results, ctr):
            
            while ctr < num_lines and len(all_results) < 1000:
                
                line_raw = opendx.readline().strip()
                line = [float(val) for val in line_raw.split()]
                all_results.extend(line)
                ctr += 1
            return ctr

        # put some results into the all_results deque
        ctr = refill(opendx, all_results, ctr)

        for x in range(grid.x_pts):
            for y in range(grid.y_pts):
                for z in range(grid.z_pts):

                    # stash the result from opendx file into matrix
                    results[x,y,z] = all_results.popleft()
                    
                    # refill our results deque from the file if it's running low
                    if (len(all_results) == 0):
                        ctr = refill(opendx, all_results, ctr)

        opendx.close()

        return grid, results
    
    @staticmethod
    def real_xyz_to_grid_idx(gridparms, realx, realy, realz):
        """Real x,y,z positions are in Angstroms, not metres."""
        
        x = math.floor((realx - gridparms.origin_x) / gridparms.deltax())
        y = math.floor((realy - gridparms.origin_y) / gridparms.deltay())
        z = math.floor((realz - gridparms.origin_z) / gridparms.deltaz())
        
        
        if (x < 0 or x > gridparms.x_pts or 
            y < 0 or y > gridparms.y_pts or 
            z < 0 or z > gridparms.z_pts):
            raise BadGridPoint
        
        return x,y,z

    @staticmethod
    def get_potential_from_file(opendx_filename, realx, realy, realz, gridparms=None):
        
        if gridparms is None:
            gridparms = copy_GridParms_from_opendx_file(opendx_filename)

        try:
            x,y,z = real_xyz_to_grid_idx(gridparms, realx, realy, realz)
        
        except BadGridPoint:
            
            print "position passed in is outside of the grid defined by this opendx file"
            return 0.0

        idx_in_file = int(x * (gridparms.y_pts*gridparms.z_pts) + \
                        y * (gridparms.z_pts) + \
                        z)
        
        line = idx_in_file / 3
        idx_within_line = idx_in_file - (line * 3)
        
        opendx = open(opendx_filename,'r')
        
        # get past the preamble
        PotentialGrid.get_preamble(opendx)

        ctr = 0
        while ctr < line:
            opendx.readline()
            ctr += 1
        line_components = opendx.readline().split()
        
        result = float(line_components[idx_within_line])
        opendx.close()

        return result

if __name__ == "__main__":
    
    opendx_filename = "barnase-barstar/barnase.dx"
    grd = copy_GridParms_from_opendx_file(opendx_filename)
    print grd.OpenDX_Preamble()

    print PotentialGrid.get_potential_from_file(opendx_filename, 0.0, 0.0, 0.0)
    print PotentialGrid.get_potential_from_file(opendx_filename, 2.0e-9, -1.0e-8, 0.0)
    
    #grid, results = read_opendx_file(opendx_filename)
    #print results[0,0,0]
    #print get_potential(results, grid, 2.0e-9, -1.0e-8, 0.0)
    
