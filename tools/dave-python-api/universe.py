# vim: set fileencoding=UTF8 :
#
# universe.py
#
# A Brownian Dynamics universe object.
# Based on octree for fast collision checking.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

import octree
import entities
from vector import *
import numpy
from constants import *
import trajectory

# Monte Carlo based on energy change
def monte_carlo(energy_change):
    
    # negative changes
    if energy_change < 0.0:
        # always accept moves which 
        # lower the overall energy
        return True
    
    import random
    import math
    p = math.exp(-energy_change) # acceptance probability
    num = random.random() # roll the dice
    return (num <= p) # accept or reject

class Universe(octree.Octree):

    """Class representing the entire Brownian Dynamics universe.
    
    This is really just an Octree root node, with convenient accessor
    functions."""

    def __init__(self, Dprotein, Dsolvent, screening_length=None, size=None, output_filepath=None):

        """Create a universe with given ionic screening and dielectrics.
        
        Ionic screening length is in Angstroms; this is converted into the
        inverse Angstrom quantity kappa."""

        # 1 object per node since this is a collision detection tree.
        if size is not None:
            super(Universe, self).__init__(worldSize=size, max_objects_per_node=1)
        else:
            super(Universe, self).__init__(max_objects_per_node=1)
        
        # store constants
        self.Dext = float(Dsolvent)
        self.Dint = float(Dprotein)
        self.epsilon = float(self.Dext / self.Dint)
        if screening_length is not None:
            self.kappa = 1.0 / screening_length
        else:
            self.kappa = 0.0 # infinite screening length

        # prepare output folder for trajectories
        if output_filepath is None:
            import tempfile
            output_filepath = tempfile.mktemp(prefix="uclbd_output_", dir=".")
            print "Output will be written to %s" %(output_filepath)
        self.trajectory_writer = trajectory.TrajectoryWriter(output_filepath, 
                                                             self.Dext, 
                                                             self.Dint, 
                                                             self.kappa)

    def addDiffusingEntity(self, de):

        xyz = de.xyz
        
        # assign diffusing entity to a volume
        try:
            new_obj = self.insertObject(de, xyz)
        except ValueError:
            # this will happen if the universe size is not big enough...
            offset = xyz - self.position
            
            cx,cy,cz = self.position
            xlim_pos = max([cx + self.edge_length, xyz[0]])
            xlim_neg = min([cx - self.edge_length, xyz[0]])
            ylim_pos = max([cy + self.edge_length, xyz[1]])
            ylim_neg = min([cy - self.edge_length, xyz[1]])
            zlim_pos = max([cz + self.edge_length, xyz[2]])
            zlim_neg = min([cz - self.edge_length, xyz[2]])
            
            new_edge_length = max([xlim_pos - xlim_neg, ylim_pos - ylim_neg, zlim_pos - zlim_neg])
            
            new_centre = Vector( (xlim_neg + xlim_pos),
                                 (ylim_neg + ylim_pos),
                                 (zlim_neg + zlim_pos) ) / 2.0
            
            copy_items = self.master_list[:]
            
            if len(copy_items) > 0:
                self.clear()
                self.edge_length = new_edge_length * 1.1 # add 10% extra length
                self.position = new_centre
            
                for item in copy_items:
                    self.insertObject(item.obj, item.position)
            else:
                # special case of an empty universe -- just shift the centre
                # to the position of the new object to be added
                self.position = xyz
            
            new_obj = self.insertObject(de, xyz)

        # resolve any clashes
        #if len(self.diffusing_entities) > 1:
        #    self.check_for_clashes([new_obj])

        # TODO: reinit octree as new_obj may have moved a bit

        return

    def reinit(self):
        """Reinitialise the octree."""
        
        contents = self._data[:]
        
        self.clear()
        for wrapped_obj in contents:
            try:
                self.insertObject(wrapped_obj.obj, wrapped_obj.obj.xyz)
            except ValueError:
                x,y,z = wrapped_obj.position - self.position
                half_width = self.edge_length / 2.0
                dx = abs(x) - half_width
                dy = abs(y) - half_width
                dz = abs(z) - half_width
                if x < -half_width:
                    x += dx
                elif x > half_width:
                    x -= dx
                if y < -half_width:
                    y += dy
                elif y > half_width:
                    y -= dy
                if z < -half_width:
                    z += dz
                elif z > half_width:
                    z -= dz
                wrapped_obj.obj._xyz = self.position + Vector(x,y,z)
                
                self.insertObject(wrapped_obj.obj, wrapped_obj.obj.xyz)
                
        return
    
    def solve_electrostatics(self):
        """Solve the electrostatics for all entities in the Universe."""
        
        # calculate electrostatics
        print "calculating boundary element method electrostatics"
        import BEM
        BEM.solveElectrostatics(self.diffusing_entities, 
                                self.kappa, 
                                self.Dext, 
                                self.Dint)
        print "done BEM calculations"
        
        return
    
    def runTimeStep(self, t):
        """t is the target timestep in nanoseconds.
        
        The actual timestep used will be adjusted to ensure force move is not too much greater than the diffusion step."""

        target_time = t
        t_accum = 0.0

        # init diffusing entity forces
        for de in self.diffusing_entities:
            try:
                de.f_trans
                de.f_rot
            except AttributeError:
                de.f_trans, de.f_rot, e = de.getForcedMove(t, 
                                                           self.kappa, 
                                                           self.Dext, 
                                                           self.diffusing_entities)
        
        # TODO: calculate the energy and decide whether to accept or reject this
        # move
        while (t_accum < target_time):

            # make sure t can't exceed next checkpoint
            t = min([t, target_time - t_accum])
            max_move = 0.0
            
            # loop over all entities
            for de in self.diffusing_entities:
    
                # diffusion move
                diff_trans, diff_rot = de.getDiffusionMove(t) # in Angstroms, in Universe coords
                
                # now get a set of force moves which, individually, are not
                # massively greater than the diffusion motion
                f_trans = de.f_trans # already calculated previously at same
                                     # time as energy eval for previous step
                f_rot = de.f_rot
                
                assert(diff_trans is not None and f_trans is not None)
                trans = diff_trans + f_trans
                
                if f_rot is None:
                    rot = diff_rot
                elif diff_rot is None:
                    rot = f_rot
                else:
                    # assume rotation due to forces, followed by diffusive
                    # rotation
                    rot = f_rot * diff_rot
                
                # commit the move to the diffusing entity
                de.apply_movement(trans, rot)
                if trans.length() > max_move:
                    max_move = trans.length()

            #MAXMOVE = 5.0
            #if max_move > MAXMOVE:
               
                #print "too fast (%f)" %(max_move)
                #for de in self.diffusing_entities:
                    #de.undo_last_movement()
                #t /= (max_move / MAXMOVE)
                #continue

            #t /= (max_move / MAXMOVE)
            
            if self.reject_timestep():
                
                print "steric clashes"
                for de in self.diffusing_entities:
                    de.undo_last_movement()
                continue
                
            # automatically reject any timestep which moves anything by more
            # than 5 Angstroms
            #if max_move > 5.0:
            #    print "rejected timestep size %f due to large force motion %f"\
            #          %(t, max_move)
            #    for de in self.diffusing_entities:
            #        de.undo_last_movement()
            #    t /= max_move
            #    continue
                
            # fix steric clashes
            #self.check_for_clashes()

            #if self.reject_timestep():
            #    print "rejected due to steric clash at stepsize=%f" %(t)
            #    for de in self.diffusing_entities:
            #        de.undo_last_movement()
            #    if delta_e > 1.0:
            #        t /= delta_e
            #    else:
            #        t *= 0.9
            #    continue

            # join anything within 3 Angstroms (1 water molecule) of any other
            # object -- desolvate the contact by merging the meshes at this
            # point
            
            # now re-solve electrostatics to get new electrostatic energies
            #self.solve_electrostatics()
            
            delta_e = 0.0
            for de in self.diffusing_entities:
                de.f_trans, de.f_rot, e_new = de.getForcedMove(t, self.kappa, self.Dext, self.diffusing_entities)
                delta_e += e_new

            print "delta_e: %f" %(delta_e)
            
            if delta_e < 0.0:
                
                print "energy lowered by %f: timer is at %f (target: %f)"\
                      "(current stepsize=%f)" %(delta_e, t_accum, target_time, t)
                t_accum += t
                #t *= 1.1

                # reset object positions in octree
                for wrapped in self.master_list:
                    wrapped._position = wrapped.obj.xyz
                self.reinit()
            
            elif monte_carlo(delta_e):
                
                print "Accepted unfavourable MC step (dE=%f): timer is at %f (target: %f)"\
                      "(current stepsize=%f)" %(delta_e, t_accum, target_time, t)
                
                t_accum += t
                # reset object positions in octree
                for wrapped in self.master_list:
                    wrapped._position = wrapped.obj.xyz
                self.reinit()
                
            else:
                
                print "Monte-Carlo rejected unfavourable energy %f (stepsize: %f)" %(delta_e, t)
                for de in self.diffusing_entities:
                    de.undo_last_movement()
                    
        return

    def reject_timestep(self):
        """Return true if configuration is bad (due to steric clashes)."""
        
        for wrapped_de in self.master_list:
    
            # the actual diffusing entity object is inside a wrapper
            de = wrapped_de.obj
            
            # find the node which is appropriate size for this entity
            node = self.find_leaf_at_position(wrapped_de.position)
            radius = de.bounding_radius
            while (node.edge_length < radius and node.parent is not None):
                node = node.parent
            
            # get all other entities within neighbourhood
            others = []
            for neighbour in node.colleagues:
                others.extend(neighbour._data)
            
            others.remove(wrapped_de) # remove this entity
            
            for chk in others:

                p1 = de.diffusion_centre_universe_coordinates
                p2 = chk.obj.diffusion_centre_universe_coordinates
                
                # if they overlap then reject this timestep
                if entities.DiffusingEntity.is_overlapping(de, chk.obj):
                    return True
        return False
            
    def check_for_clashes(self, specific_obj=None):
        """Check that clashing entities within the Universe octree."""
        ctr = 0

        print "Resolving steric clashes"
        dirty_octree = False
        
        while True:
            ctr += 1
            
            # loop until all clashes resolved

            # this variable gets set to True if we find any clashes
            # if it's still False by the end of this loop, we break
            # (because our job is done!)
            clashes = False
            
            # init clash lists within wrapped objects (bit hacky this)
            for wrapped_de in self.master_list:
                wrapped_de.check_list = []
                wrapped_de.corrections = []

            if specific_obj is None:
                list_of_objs_to_check = self.master_list
            else:
                list_of_objs_to_check = specific_obj[:]

            print "Clash iteration %d (%d entities to check)" \
                  %(ctr, len(list_of_objs_to_check))
                
            for wrapped_de in list_of_objs_to_check:
                
                # the actual diffusing entity object is inside a wrapper
                de = wrapped_de.obj
                
                # find the node which is appropriate size for this entity
                node = self.find_leaf_at_position(wrapped_de.position)
                radius = de.bounding_radius
                while (node.edge_length < radius and node.parent is not None):
                    node = node.parent
                
                # get all other entities within neighbourhood
                others = []
                for neighbour in node.colleagues:
                    others.extend(neighbour._data)
                
                others.remove(wrapped_de) # remove this entity
                
                for chk in others:
                    if chk not in wrapped_de.check_list:

                        p1 = de.diffusion_centre_universe_coordinates
                        p2 = chk.obj.diffusion_centre_universe_coordinates
                        
                        # get overlapping atoms
                        #print "get overlaps..."
                        overlaps=entities.DiffusingEntity.get_overlaps(de, 
                                                                       chk.obj)
                        
                        ## check the overlaps algorithm is working!
                        #if (len(overlaps) == 0):
                            #assert((p1 - p2).length() > de.bounding_radius + chk.obj.bounding_radius)
                        
                        #print "got %d overlaps" %(len(overlaps))
                        
                        # apply correction forces at the overlaps
                        for clash_de, clash_chk, mag in overlaps:
                            wrapped_de.corrections.append((clash_de, clash_chk, chk, mag))
                            chk.corrections.append( (clash_chk, clash_de, wrapped_de, mag) )
                            clashes = True
                        
                        # record the fact that we have done the overlap check
                        # between these two entities, so we don't need to do it
                        # again
                        chk.check_list.append(wrapped_de)
                        wrapped_de.check_list.append(chk)
            
            specific_obj = []
            # loop over all the wrapped entities which have corrections to apply
            for wrapped_de in [w_e for w_e in list_of_objs_to_check 
                               if len(w_e.corrections) > 0]:

                dirty_octree = True
                
                # unwrap the diffusing entity from octree data wrapper
                de = wrapped_de.obj
                
                resultant_force = Vector(0.0,0.0,0.0)
                resultant_axis = Vector(0.0,0.0,0.0)

                mags = []
                for loc, other, ent, mag in wrapped_de.corrections:

                    mags.append(mag)
                    
                    # calculate force on clashing atom
                    if (loc-other).length() > 0.0:
                        force = (loc - other).normal() * mag
                    else:
                        force = _rand_xyz().normal() * mag
                    
                    # make sure the force is going to actually improve the
                    # inter-body distance
                    #r = (ent.diffusion_centre_universe_coordinates - de.diffusion_centre_universe_coordinates)
                    #if force * r > 0:
                    #    force = -force
                    #    other, loc = loc, other
                    
                    c_to_f = (loc - de.diffusion_centre_universe_coordinates)
                    if c_to_f.length() > 0:
                        resultant_axis += c_to_f.normal().cross(force.normal())
                    
                    # add this elemental force to the resultant force vector
                    resultant_force += force
                
                # target_displacement is the max overlap distance required
                target_displacement = min(mags)
                #target_displacement = len(wrapped_de.corrections) ** 0.333
                #target_displacement = 0.1 * (1+2*(ctr/10)) # slowly increase the shimmy
                
                # apply translation and rotation composed of the possibly
                # multiple correction force/torques. Translate by 1
                # Angstrom; ensure that the rotation doesn't cause any gross
                # movement of more than 1 Anstrom.
                if resultant_force.length() > 0:
                    resultant_force = resultant_force.normal()                

                trans = resultant_force * target_displacement
                if resultant_axis.length() > 0.0:
                    rot = Rotation(resultant_axis.normal(), 
                                   target_displacement/de.bounding_radius).asQuaternion()
                    rot = None
                else:
                    rot = None
                
                de.apply_movement(trans, rot)
                if wrapped_de not in specific_obj:
                    specific_obj.append(wrapped_de)
                for other_checks in wrapped_de.check_list:
                    if other_checks not in specific_obj:
                        specific_obj.append(other_checks)
                
            if not clashes: break
        
        # clear the temporary holders in wrapped objects
        for wrapped_de in self.master_list:
            del(wrapped_de.check_list)
            del(wrapped_de.corrections)

        # octree may have been disrupted by moving stuff about
        #if dirty_octree:
        #    self.reinit()

        # TODO: Debug logging
        print "done clash checks"
        return
    
    @property
    def diffusing_entities(self):
        """Return a list of all diffusing entities in the universe.
        
        Get this from the master list of stuff in the tree."""

        return [wrapped_obj.obj for wrapped_obj in self.master_list]
    
    def snapshot_kinemage(self, filename):
        """Output snapshot of the Universe as a kinemage."""
        
        kin = open(filename, 'w')
        print >>kin, "@kinemage"
        for de in self.diffusing_entities:
            de.mesh.kinemage(f=kin, normals=True, charges=True,
                             centre_of_rotation=de.centre_of_diffusion,
                             rotation=de.rotation, 
                             translation=de.xyz)
        kin.close()
        
        return
    
    def snapshot_pdb(self, pdb_filename):
        """Output a snapshot of the Universe in PDB format."""

        pdb_file = open(pdb_filename, 'w')

        for de in self.diffusing_entities:
            de.writePDB(pdb_file)

        pdb_file.close()
        
        return

    def open_dcdfile(self, dcd_filename="trajectory"):
        
        # number of atoms is the total number of visualisation points in all
        # entities
        num_atoms = sum([len(de.vis_points) for de in self.diffusing_entities])
        from MDAnalysis import DCD
        self.writer = DCD.DCDWriter(dcd_filename + ".dcd", num_atoms)

        f = open(dcd_filename+".xyz", 'w')
        f.write("%d\n" %(num_atoms))
        f.write("Simulation geometry.\n")
        i = 0
        for de in self.diffusing_entities:
            for atom in de.atoms:
                x,y,z = de.convert_local_coordinate(atom.position)
                f.write("%s %f %f %f\n" %(atom.name, x, y, z))
                i += 1
        f.close()
    
    def close_dcdfile(self):
        try:
            self.writer.close_trajectory()
            self.writer = None
        except AttributeError:
            pass

    def snapshot(self, time):
        
        self.trajectory_writer.add_trajectory_point(time, self.diffusing_entities)
        
    def snapshot_DCD(self):

        from MDAnalysis.DCD import Timestep

        # get full list of visualisation points (representative points of
        # diffusing entities; could be atoms, or just surface points
        vis_points = []
        for de in self.diffusing_entities:
            vis_points.extend(de.vis_points)
            
        pos = numpy.array(vis_points,'f')
        ts = Timestep(pos)
        ts.numatoms = len(vis_points)
        self.writer.write_next_timestep(ts)
        return

#
# TEST FUNCTIONS
#

# helper function for test purposes
import random
def _rand_xyz(de_list=[],edge_length=100):
    """Generate a random vector within box of given edge_length."""
    while True:
        
        clean = True
        v = Vector([(random.random() - 0.5)*edge_length for i in range(3)])
        for de in de_list:
            if (v - de.xyz).length() < 5.0:
                clean = False
                
        if clean:
            break
    return v

from geometry import _rand_rot

def _rand_charge():
    """Generate random charge (+1 or -1)."""
    if random.random() > 0.5:
        return +1
    else:
        return -1

def _salt_test(conc, num_timesteps, timestep=0.01, edge_length=100.0):
    """Generate a test set of Na+ and Cl- ions at given concentration (mM)."""
    
    avogadro = 6.022e23
    vol = (edge_length * 1e-10) ** 3
    num_entities = int(conc * 1e-3 * avogadro * 1e3 * vol)

    print "%d Na+ and Cl- ions in %dx%dx%d (Angstroms) box, concentration=%dmM"\
          %(num_entities, edge_length, edge_length, edge_length, conc)
    
    from vector import Vector, Rotation, Quaternion

    uni = Universe(1.0, 80.0, screening_length=None, size=edge_length)
    origin = Vector(0.0, 0.0, 0.0) # the origin
    rot = Rotation(Vector(1,1,1), 0.0) # no rotation

    from mesh import MeshedSphere
    
    for i in range(num_entities):
        
        # create diffusing entity
        de = entities.DiffusingEntity()
        de.set_isotropic_diffusion(70.0, 0.0) # ions diffuse at approx. 70 Angstroms^2/ns
        de.set_xyz(_rand_xyz(de_list=uni.diffusing_entities,edge_length=edge_length))
        chloride_ion = entities.Charge(origin, 
                                       charge=-1, 
                                       name="Cl")
        de.atoms = de.charges = [chloride_ion]
        de.mesh = MeshedSphere(radius=chloride_ion.radius, num_subdivides=3)
        de.mesh._charges = [chloride_ion]
        
        # add ion to universe
        uni.addDiffusingEntity(de)

    for i in range(num_entities):
    
        # create diffusing entity
        de = entities.DiffusingEntity()
        de.set_isotropic_diffusion(70.0, 0.0) # ions diffuse at approx. 70 Angstroms^2/ns
        de.set_xyz(_rand_xyz(de_list=uni.diffusing_entities,edge_length=edge_length))
        sodium_ion = entities.Charge(origin, 
                                     charge=+1, 
                                     name="Na")
        de.atoms = de.charges = [sodium_ion]
        de.mesh = MeshedSphere(radius=sodium_ion.radius, num_subdivides=3)
        de.mesh._charges = [sodium_ion]
    
        # add ion to universe
        uni.addDiffusingEntity(de)        
        
    # resolve clashes
    uni.check_for_clashes()
#    uni.solve_electrostatics()
    
    #uni.snapshot_kinemage("NaCl_%dmM_mesh.kin" %(conc))
    
    # open trajectory file for snapshots
    #uni.open_dcdfile("NaCl_%dmM_trajectory" %(conc))
    sim_clock = 0.0 # simulation clock
    uni.snapshot(sim_clock)  # init configuration
    for i in range(num_timesteps):
        print "Running timestep %d (of %d)" %(i+1, num_timesteps)
        uni.runTimeStep(timestep) # timestep is in nanoseconds
        sim_clock += timestep
        print "Taking snapshot"
        uni.snapshot(sim_clock)        
    #uni.close_dcdfile()
    
    return

def _ion_pair(num_timesteps, timestep=0.1):
    
    from vector import Vector, Rotation, Quaternion

    size = 10.0
    rad = 1.0
    detail = 1
    
    uni = Universe(2.0, 80.0, screening_length=None, size=size)
    origin = Vector(0.0, 0.0, 0.0) # the origin

    from mesh import MeshedSphere
    from math import pi

    # ion 1
    # create diffusing entity
    de = entities.DiffusingEntity()
    de.set_isotropic_diffusion(70.0, 0.0) # ions diffuse at approx. 70 Angstroms^2/ns
    de.set_xyz(_rand_xyz(edge_length=size))
    de.atoms = [entities.Atom(Vector(0.0,0.0,0.0), rad)] 
    de.mesh = MeshedSphere(radius=rad, num_subdivides=detail)
    
    # correct mesh areas to better approximate sphere
    #correction = (4.0*pi*rad*rad) / sum([t.area for t in de.mesh.triangles])
    #for t in de.mesh._triangles:
    #    t._area *= correction
        
    de.mesh._charges = [entities.Charge(Vector(0.0,0.0,0.0), +1, radius=rad)]
    uni.addDiffusingEntity(de)
    de1 = de

    # ion 2
    # create diffusing entity
    de = entities.DiffusingEntity()
    de.set_isotropic_diffusion(70.0, 0.0) # ions diffuse at approx. 70 Angstroms^2/ns
    de.set_xyz(_rand_xyz(edge_length=size))
    de.atoms = [entities.Atom(Vector(0.0,0.0,0.0), rad)] 
    de.mesh = MeshedSphere(radius=rad, num_subdivides=detail)
    
    # correct mesh areas to better approximate sphere
    #correction = (4.0*pi*rad*rad) / sum([t.area for t in de.mesh.triangles])
    #for t in de.mesh._triangles:
    #    t._area *= correction
        
    de.mesh._charges = [entities.Charge(Vector(0.0,0.0,0.0), -1, radius=rad)]
    uni.addDiffusingEntity(de)
    de2 = de
    
    # resolve clashes
    uni.check_for_clashes()
    
    uni.snapshot_kinemage("ion_pair_mesh.kin")
    
    # open trajectory file for snapshots
    uni.open_dcdfile("ion_pair_trajectory")        

    from math import exp
    from constants import force_conversion_Newtons, force_conversion_kT_per_Angstrom
    
    for i in range(num_timesteps):
        print "Running timestep %d (of %d)" %(i+1, num_timesteps)

        uni.solve_electrostatics()
        
        # get distance between ions at start of timestep
        dist = (de1.xyz - de2.xyz).length()
        
        # print out the electrostatic forces
        f1 = de1.mesh.calculate_force(uni.kappa, uni.Dext, si_units=False)
        f2 = de2.mesh.calculate_force(uni.kappa, uni.Dext, si_units=False)
        actual = force_conversion_kT_per_Angstrom * exp(-uni.kappa * dist) / \
                   (4.0*pi*uni.Dext*dist*dist)
        print "FORCES: ", actual, f1.length(), f2.length(), abs(actual - f1.length())*100.0 / actual, f1.length() / actual

        # now run a timestep
        uni.runTimeStep(timestep) # timestep is in nanoseconds
        
        print "Taking snapshot"
        uni.snapshot(clock)
        
    uni.close_dcdfile()
    
    return
    
def _ubiquitin_test(num_entities, num_timesteps, timestep=0.1):
    
    from vector import Vector, Rotation, Quaternion
    import random

    #filename = "ubiquitin/ubiquitin"
    filename = "ecm_model"

    uni = Universe(2.0, 80.0, screening_length=None, size=243)
    origin = Vector(0.0, 0.0, 0.0) # the origin

    #from Scientific.IO.PDB import Structure
    #struct = Structure(filename + ".pdb")
    
    for i in range(num_entities):
        #sphere = MeshedSphere(radius=1.0, xyz=origin, num_subdivides=1)
        
        # create diffusing entity
        de = entities.DiffusingEntity()
        de.set_xyz(_rand_xyz(edge_length=243)) # 243 for 10 entities @ 10mg/ml; 534 for 100
        de.set_rotation(_rand_rot())
        #de.set_hydropro_diffusion_tensor(filename + ".res")
        de.set_isotropic_diffusion(25,0.0)
        de.set_charges_from_pqr(filename + ".pqr", set_atoms_from_pqr=True)
        de.create_effective_charges(filename)
        #de.set_atoms_from_pdb_structure(struct)
        #de.create_mesh()
        
        # centre of diffusion not correct from hydropro???
        centre = Vector(0.0,0.0,0.0)
        for atom in de.atoms:
            centre += atom.position
        de.set_centre_of_diffusion(centre / len(de.atoms))

        uni.addDiffusingEntity(de)
        print "Added entity %d" %(len(uni.diffusing_entities))
        #print i, len(uni.diffusing_entities)

    #uni.snapshot_kinemage("mesh.kin")

    #uni.open_dcdfile("ubiquitin_clash_trajectory")                            
    # resolve clashes
    #uni.open_dcdfile("ubiquitin_mesh_test")
    uni.check_for_clashes()    
    #uni.open_dcdfile("ubiquitin")
    #uni.close_dcdfile()
    
    sim_clock = 0.0 # simulation clock
    for i in range(num_timesteps):
        print "Running timestep %d" %(i+1)
        uni.runTimeStep(timestep) # timestep is in nanoseconds
        sim_clock += timestep
        print "Taking snapshot"
        uni.snapshot(sim_clock)
    #uni.close_dcdfile()    
    
#
# RUN TESTS
#
if __name__ == "__main__":

    #_salt_test(100,100) # 200mM salt, 100 timesteps at 0.05ns per step
    #_ion_pair(num_timesteps=100, timestep=0.01)
    _ubiquitin_test(num_entities=50, num_timesteps=100, timestep=0.1)
