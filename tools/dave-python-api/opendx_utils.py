from _BEM import GridParms
import math

def copy_GridParms_from_opendx_file(filename):
    """Copy the grid layout from an existing OpenDX file."""

    # open the openDX file
    opendx_file = open(filename,'r')

    # get the preamble lines
    preamble = get_preamble(opendx_file)
    
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
    
def get_preamble(opendx_file_handle):
    
    preamble = []

    while len(preamble) < 7:
        line = opendx_file_handle.readline()
        if line[0] != '#':
            preamble.append(line.strip())

    return preamble

def read_opendx_file(opendx_filename):
    
    grid = copy_GridParms_from_opendx_file(opendx_filename)
    
    from numpy import zeros
    results = zeros((grid.x_pts,grid.y_pts,grid.z_pts),'d')
    num_lines = math.ceil(grid.num_grid_points() / 3.0)
    
    opendx = open(opendx_filename)
    get_preamble(opendx) # skip preamble lines

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
            for z in range(grid.x_pts):

                # stash the result from opendx file into matrix
                results[x,y,z] = all_results.popleft()
                
                # refill our results deque from the file if it's running low
                if (len(all_results) == 0):
                    ctr = refill(opendx, all_results, ctr)

    opendx.close()

    return grid, results

def write_opendx_file(filename, grid, data):
    
    opendx_file = open(filename, 'w')
    
    opendx_file.write(grid.OpenDX_Preamble())
    
    def print_out_buffer(buf):

        out_str = ""
        for val in buf:
            out_str += "%01.6e " %(val)
        
        print >>opendx_file, out_str.strip()
        return
        
    out_buffer = []
    for ii in range(grid.x_pts):
        for jj in range(grid.y_pts):
            for kk in range(grid.z_pts):
                
                out_buffer.append(data[ii,jj,kk])
                
                if len(out_buffer) == 3:
                    print_out_buffer(out_buffer)
                    out_buffer = []

    # print any last remaining vals in the buffer!
    print_out_buffer(out_buffer)
    
    opendx_file.write(grid.OpenDX_Suffix())
    opendx_file.close()
    
    return

class BadGridPoint(Exception):
    def __init__(self, *args):
        Exception.__init__(self, *args)
    
def real_xyz_to_grid_idx(gridparms, realx, realy, realz):
    """Real x,y,z positions are in Angstroms, not metres."""
    
    x = round((realx - gridparms.origin_x) / gridparms.deltax())
    y = round((realy - gridparms.origin_y) / gridparms.deltay())
    z = round((realz - gridparms.origin_z) / gridparms.deltaz())
    
    if (x < 0 or x >= gridparms.x_pts or 
        y < 0 or y >= gridparms.y_pts or 
        z < 0 or z >= gridparms.z_pts):
        
        raise BadGridPoint
    
    return x,y,z

def grid_idx_to_real_xyz(gridparms, gx, gy, gz):

    realx = gridparms.origin_x + gx*gridparms.deltax()
    realy = gridparms.origin_y + gy*gridparms.deltay()
    realz = gridparms.origin_z + gz*gridparms.deltaz()
    
    return realx, realy, realz

def get_field(opendx_matrix, opendx_grid, realx, realy, realz):
    """Returns potential and field at the given point."""
    
    from Scientific.Geometry import Vector
    
    try:
        x,y,z = real_xyz_to_grid_idx(opendx_grid, realx, realy, realz)
    
        if (x+2 >= opendx_grid.x_pts or
            y+2 >= opendx_grid.y_pts or
            z+2 >= opendx_grid.z_pts):
            raise BadGridPoint
    except BadGridPoint:
        
        #print "position passed in is outside of the grid defined by this opendx file"
        return 0.0, Vector(0.0, 0.0, 0.0)
    
    locx, locy, locz = grid_idx_to_real_xyz(opendx_grid,x,y,z)
    deltx = 2.0*opendx_grid.deltax()
    delty = 2.0*opendx_grid.deltay()
    deltz = 2.0*opendx_grid.deltaz()
    
    def dxdydz(_x,_y,_z):
        dx = (opendx_matrix[_x+1,_y,_z] - opendx_matrix[_x-1,_y,_z]) / deltx
        dy = (opendx_matrix[_x,_y+1,_z] - opendx_matrix[_x,_y-1,_z]) / delty
        dz = (opendx_matrix[_x,_y,_z+1] - opendx_matrix[_x,_y,_z-1]) / deltz
        return Vector(dx, dy, dz)
    
    d2x = (dxdydz(x+1,y,z) - dxdydz(x-1,y,z)) / deltx
    d2y = (dxdydz(x,y+1,z) - dxdydz(x,y-1,z)) / delty
    d2z = (dxdydz(x,y,z+1) - dxdydz(x,y,z-1)) / deltz

    dx, dy, dz = closest_field = dxdydz(x,y,z)
    
    closest_pot = opendx_matrix[x,y,z]
    interp_pot = closest_pot + (realx-locx)*dx \
                             + (realy-locy)*dy \
                             + (realz-locz)*dz
    
    interp_field = closest_field + (realx-locx)*d2x \
                                 + (realy-locy)*d2y \
                                 + (realz-locz)*d2z
    
    return interp_pot, interp_field
        
def get_potential(opendx_matrix, opendx_grid, realx, realy, realz):
    """Linearly interpolates around the nearest grid point."""
    
    try:
        x,y,z = real_xyz_to_grid_idx(opendx_grid, realx, realy, realz)
    
    except BadGridPoint:
        
        print "position passed in is outside of the grid defined by this opendx file"
        return 0.0

    locx, locy, locz = grid_idx_to_real_xyz(opendx_grid,x,y,z)
    
    dx = (opendx_matrix[x+1,y,z] - opendx_matrix[x-1,y,z]) / (2.0*opendx_grid.deltax())
    dy = (opendx_matrix[x,y+1,z] - opendx_matrix[x,y-1,z]) / (2.0*opendx_grid.deltay())
    dz = (opendx_matrix[x,y,z+1] - opendx_matrix[x,y,z-1]) / (2.0*opendx_grid.deltaz())

    closest_pt = opendx_matrix[x,y,z]
    interpolated = closest_pt + (realx-locx)*dx + (realy-locy)*dy + (realz-locz)*dz
    #print closest_pt, interpolated
    
    return interpolated

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
    get_preamble(opendx)

    ctr = 0
    while ctr < line:
        opendx.readline()
        ctr += 1
    line_components = opendx.readline().split()
    
    result = float(line_components[idx_within_line])
    opendx.close()

    return result

class OpenDXGrid(object):
    
    __dx_cache = {}

    class OpenDXGrid_Item(object):
        
        def __init__(self, filename):
            
            print "DEBUG: reading %s" %(filename+".dx")
            gridparms, opendx = read_opendx_file(filename+".dx")
            self.__gridparms = gridparms
            self.__opendx = opendx
    
        def get_potential_at(self, xyz):
            """Get the potential at the xyz location within the grid."""
            
            realx, realy, realz = xyz
            return get_potential(self.__opendx, 
                                 self.__gridparms, 
                                 realx, realy, realz)
    
        def get_field_at(self, xyz):
            
            realx, realy, realz = xyz
            return get_field(self.__opendx, 
                             self.__gridparms,
                             realx, realy, realz)
    
    @classmethod
    def get_ecm_model(cls, filename):
        
        try:
            return cls.__dx_cache[filename]
        except KeyError:
            new_dx_grid = cls.OpenDXGrid_Item(filename)
            cls.__dx_cache[filename] = new_dx_grid
        
        return cls.__dx_cache[filename]

if __name__ == "__main__":
    
    opendx_filename = "ecm_model.dx"
    grd = copy_GridParms_from_opendx_file(opendx_filename)
    print grd.OpenDX_Preamble()

    print get_potential_from_file(opendx_filename, 10.0, 10.0, 0.0)
    print get_potential_from_file(opendx_filename, 20, -10, 0.0)
    
    #grid, results = read_opendx_file(opendx_filename)
    #print results[0,0,0]
    #print get_potential(results, grid, 2.0e-9, -1.0e-8, 0.0)
    