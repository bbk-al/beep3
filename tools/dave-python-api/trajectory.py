import os
import shelve # persistent object type thing
from Scientific.Geometry import Vector, Quaternion
from geometry import local_to_universal

class TrajectoryPoint(list):
    
    """This represents a point in trajectory space -- the collection of
    rotations and translations which define the locations of all molecules in
    the system at a given timestep."""

    def __init__(self, time, entities):
        
        self.time = time
        for de in entities:
            self.append( (de.id, de.xyz, de.rotation) )
 
        return

class TrajectoryBase(object):
    
    def __init__(self, trajectory_foldername):
        
        self.trajectory_foldername = trajectory_foldername
        assert(os.path.exists(self.trajectory_foldername))

        # persistent objects to represent the entities in the simulation and
        # the trajectories of all diffusing entities
        self.sim_metadata = shelve.open(self.mkpath("sim_metadata"))        
        self.entity_dictionary = shelve.open(self.mkpath("entity_dictionary"))
        self.trajectory_points = shelve.open(self.mkpath("trajectory_points"), protocol=2, writeback=True)

    def mkpath(self, some_path):
        """Utility function to prepend the trajectory foldername onto paths."""
        return "%s/%s" %(self.trajectory_foldername, some_path)
        
class TrajectoryWriter(TrajectoryBase):
    
    def __init__(self, trajectory_foldername, Dext, Dint, kappa):
        
        # ensure that the output folder doesn't already exist
        if (os.path.exists(trajectory_foldername)):
            raise Exception
        os.mkdir(trajectory_foldername)

        super(TrajectoryWriter, self).__init__(trajectory_foldername)        
        
        # store the simulation metadata
        self.sim_metadata["Dext"] = Dext
        self.sim_metadata["Dint"] = Dint
        self.sim_metadata["kappa"] = kappa
        
        if not self.trajectory_points.has_key("trajectory"):
            self.trajectory_points["trajectory"] = []

    def add_trajectory_point(self, time, diffusing_entities):

        # check that the entities are defined in the entity dictionary
        for de in diffusing_entities:
            if (de.id not in self.entity_dictionary.keys()):
                self.entity_dictionary["%d" %(de.id)] = (de.centre_of_diffusion, de.raw_vis_points, de.bounding_radius)
                
                #pqr_archived_path = os.path.basename(os.path.normpath(self.mkpath(de.pqr_filename)))
                #if not os.path.exists(pqr_archived_path):
                    #import shutil
                    #print de.pqr_filename, pqr_archived_path
                    #shutil.copyfile(de.pqr_filename, pqr_archived_path)

        #traj_key = "%09d" %(len(self.trajectory_points.keys())+1)
        self.trajectory_points["trajectory"].append(TrajectoryPoint(time, diffusing_entities))
        self.trajectory_points.sync()
        
        return

class TrajectoryReader(TrajectoryBase):
   
    def __init__(self, trajectory_foldername):

        super(TrajectoryReader, self).__init__(trajectory_foldername)        
        os.chdir(trajectory_foldername)
        
        from vtk import VTK_XML_Serial_Unstructured
        vtk_writer = VTK_XML_Serial_Unstructured()
        
        for particle_list in self.trajectory_points["trajectory"]:
            #print particle_list.time
            xx,yy,zz,radii = [], [], [], []
            for (id, xyz, rot) in particle_list:

                # extract visualisation points
                centre, vis_pts, max_radius = self.entity_dictionary["%d" %(id)]

                #x,y,z = xyz
                #xx.append(x)
                #yy.append(y)
                #zz.append(z)
                #radii.append(max_radius)
                
                #ave = 0.0
                #for v in vis_pts:
                    #ave += (v.position-centre)
                #ave /= len(vis_pts)
                for v in vis_pts:
                    x,y,z = local_to_universal(v.position, 
                                               ref_pt_local=centre,
                                               ref_pt_universe=xyz,
                                               rotation=rot)
                    xx.append(x)
                    yy.append(y)
                    zz.append(z)
                    radii.append(v.radius)
                
            vtk_writer.snapshot("%f.vtu" %(particle_list.time), xx, yy, zz, radii=radii)
   
        vtk_writer.writePVD("trajectory.pvd")
            
if __name__=="__main__":
    
    import sys
    
    TrajectoryReader(sys.argv[1])
    
