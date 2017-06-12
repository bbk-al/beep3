from vector import Vector
from string import atof, atoi

class MeshDef(object):
    
    """Holds the filenames for defining a mesh."""
    
    def __init__(self,
                 gts_file,
                 xyzq_file,
                 apbs_dx_file,
                 centre_file,
                 fh_vals_file,
                 shape_file):
        self.gts_file = gts_file
        self.xyzq_file = xyzq_file
        self.apbs_dx_file = apbs_dx_file
        self.centre_file = centre_file
        self.fh_vals_file = fh_vals_file
        self.shape_file = shape_file

class SimDef(object):
    
    def __init__(self, mesh_def, xyz_offset):
        
        self.mesh_def = mesh_def
        self.xyz_offset = xyz_offset

    def addToGlobalMesh(self, global_mesh_obj):
        
        # open the gts file

        # get the global_mesh_obj current vertex counter
        
        # add the vertices to the global mesh
        
        # add the triangles
        
        # add the charges
        
        return
    
class ConfigFile(object):
    """A class to contain the data within a config file.
    
    TODO: convert this to an xml format rather than shonky flat file."""
    
    def __init__(self, filename):
        
        config = open(filename, 'r')
        lines = config.readlines()
        config.close()
        
        from string import find, split
        new_lines = []
        for line in lines:
            idx = find(line, "//")
            if (idx != -1):
                line = line[:idx]
            new_lines.append(line.strip())
        lines = new_lines[:]
        
        self.edge_length = atof(lines.pop(0))
        self.maxsize = atof(lines.pop(0))
        self.bottom_chare_level = atoi(lines.pop(0))
        self.kappa, self.Dint, self.Dext = [atof(xx) 
                                            for xx in lines.pop(0).split()]
        self.num_mesh_defs = atoi(lines.pop(0))
        self.mesh_defs = []
        
        for mesh_def_ctr in range(self.num_mesh_defs):
            
            self.mesh_defs.append( MeshDef(lines.pop(0), #gts
                                           lines.pop(0), #xyzq
                                           lines.pop(0), #apbs dx grid
                                           lines.pop(0), #centre
                                           lines.pop(0), #fh vals
                                           lines.pop(0)) ) #shape definition
            
        self.num_molecules_in_sim = atoi(lines.pop(0))
        self.sim_defs = []
        for sim_ctr in range(self.num_molecules_in_sim):
            line_def = lines.pop(0).split()
            mesh_id = atoi(line_def[0])
            assert(mesh_id < self.num_mesh_defs) # assert mesh_id is ok
            xyz = Vector( atof(line_def[1]), 
                          atof(line_def[2]), 
                          atof(line_def[3]) )
            self.sim_defs.append( SimDef(self.mesh_defs[mesh_id], xyz) ) 
        
        assert(lines.pop(0).strip() == "END")
        
