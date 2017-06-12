#
# entities.py
#
# This file is for entities within the Brownian Dynamics universe -- such as
# diffusing entities, individual charges ... whatever.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

# TODO: put in some test cases!

# standard python libs
import math
import random
import string
import StringIO

# some numpy stuff which we need for constructing a DiffusingEntity
from numpy import matrix, mat, linalg, eye, zeros, array

# Scientific library imports
#from Scientific.Geometry.Quaternion import Quaternion
#from Scientific.Geometry.Transformation import Rotation
from vector import *

# required for coordinate conversions from local to universe frame
from geometry import apply_quaternion_to_vector, local_to_universal
import octree

class Atom(object):

    """Class representing a point charge (with a radius property)."""
    
    # dictionary of atomic radii
    atomic_radii = { 'C':1.9,
                     'N':1.6,
                     'O':1.5,
                     'S':1.9,
                     'H':1.0,
                     'Na':1.02,
                     'Cl':1.81 }
    
    def __init__(self, position, radius=None, name=None):
        """Construct a charge object from given position, charge and radius."""
        
        # convert position to Vector object, if it's not one already
        if isinstance(position, Vector):
            self._position = position
        else:
            self._position = Vector(position)

        # set the atom name
        if name is None:
            self._name = "atom"
        else:
            self._name = name
        
        # atom radius - get from atomic_radii dictionary if the name is valid;
        # if radius is supplied as an argument then that is used
        if radius is None:
            try:
                self._radius = Atom.atomic_radii[name]
            except KeyError:
                self._radius = 0.0 # default zero radius
        else:
            self._radius = radius

    @property
    def name(self):
        """Return name of the atom."""
        return self._name
        
    @property
    def position(self):
        """Return the position of this charge."""
        return self._position
    
    @property
    def radius(self):
        """Return the radius of this charge."""
        return self._radius

    def kinemage(self, centre_of_rotation=None, rotation=None, translation=None):
        """Return kinemage stylee representation of this Atom.
        
        e.g. 'r=radius {} x y z'
        
        """
        
        from geometry import local_to_universal
        x,y,z = local_to_universal(self.position, centre_of_rotation, translation, rotation)
        return "r=%f {} %f %f %f" %(self.radius,x,y,z)
        
    
    def __str__(self):
        """String representation of Atom."""
        return "Atom '%s' at position %s, radius %f" %(self.name, self.position, self.radius)
    
class Charge(Atom):

    """Class representing a point charge (with a radius property)."""
    
    def __init__(self, position, charge, radius=None, name=None):
        """Construct a charge object from given position, charge and radius."""
        
        super(Charge, self).__init__(position, radius=radius, name=name)
        self._charge = charge

    @property
    def charge(self):
        """Return the magnitude of this charge."""
        return self._charge
    
    def set_charge(self, val):
        self._charge = val
        return
    
    def __repr__(self):
        return "Charge magnitude %f at %s" %(self.charge, self.position)

    @staticmethod
    def writeChargesAsPQR(charges, filename=None,):
        """Write a set of charges in pseudo-PQR format.
        
        If filename is not passed in, returns a string."""

        from pdbtools import to_pdb_format
        from string import join
        
        lines = []
        for num, ch in enumerate(charges):

            x,y,z = ch.position
            q = ch.charge
            rad = ch.radius
            pdb_line = ["ATOM",num+1,"X","XXX",1,x,y,z,q,rad]
            
            lines.append(to_pdb_format(pdb_line))
        
        # turn list of lines into a string
        strung = join(lines, "\n")
        
        if filename is not None:
            open(filename, "w").write(strung)
            return
        else:
            return strung
    
class DiffusingEntity(object):

    """Class to represent a diffusing entity within a Brownian Dynamics
    simulation.
    
    A DiffusingEntity has some diffusion characteristics (diffusion tensor)
    and possibly some atom coordinates and a representation of the surface.
    
    The DiffusingEntity should have some electrostatics attribute registered
    which will hook into the electrostatics code (which should be independent
    of the implementation details of the DiffusingEntity and it's internal
    representation."""
    
    # total number of translative and rotational dimensions
    dimensions = 6
    
    # unique id generator
    __uid_ctr = 0
    
    @classmethod
    def get_uid(cls):
        cls.__uid_ctr += 1
        return cls.__uid_ctr

    def __init__(self):

        #filename = None
        #surface_mesh = None
        self.mesh = None
        self._xyz = Vector(0.0,0.0,0.0)
        self._rotation = Rotation(Vector(1,0,0), 0.0).asQuaternion()
        self.id = self.get_uid() # get an id unique for this simulation
        
        # the centre of diffusion of the entity is the point about which
        # rotations occur
        self._centre_of_diffusion = Vector(0.0,0.0,0.0)
        
        # atoms is the set of 'atoms' which make up this protein -- used for
        # surface meshing, collision detection, and trajectory visualisation
        self.atoms = []
        self.charges = []

        # this is the spherical-octree for multi-level collision checking
        self._collision_tree = None
        self._centre = None
        
        # guranteed clash radius -- anything within this radius must be in
        # contact
        self._guaranteed_clash_radius = None
        
    @property
    def centre_of_diffusion(self):
        return self._centre_of_diffusion
    
    def create_mesh(self, filename=None):
        """Create a molecular surface mesh from GTS file."""
        
        import mesh
        
        self.mesh = mesh.GTS_Mesh(filename+".gts")
        #self.mesh.enforce_no_charge_clashes(self.atoms)
        ##self.mesh.gts_cleanup()
        #while True:
            #self.mesh.gts_coarsen(0)
            #if not self.mesh.enforce_no_charge_clashes(self.atoms): break
            
        #self.mesh = mesh.SphericalMesh(atom_list=self.atoms)
        #self.mesh = mesh.TriangularPyramidMesh(atom_list=self.atoms)
        #self.mesh.coalesce(0.1)
        #self.mesh.expand(2.0)
        self.mesh._charges = self.charges[:]
        self.mesh.init_fh(Vector(0.0,0.0,0.0), Vector(0.0,0.0,0.0), None)
        return 
        
    @property
    def centre(self):
        """Return geometric centre of the atoms collection, in local coords."""

        # centre is cached
        if self._centre is None:
            cent = Vector(0.0,0.0,0.0)
            for at in self.atoms:
                cent += at.position
                
            self._centre = cent / len(self.atoms)
            
        return self._centre
        
    @property
    def guaranteed_clash_radius(self):
        
        # cached value -- derived from the minimum atom distance in the
        # collision tree (which only contains surface-like atoms)
        if self._guaranteed_clash_radius is None:
            r = min([(atom.position - self.centre_of_diffusion).length() + atom.obj.radius 
                     for atom in self.collision_tree.master_list])
            self._guaranteed_clash_radius = r
        
        return self._guaranteed_clash_radius
    
    @property
    def collision_tree(self):
        """Octree containing all component atoms of the entity.
        
        Used as a hierarchical method for checking overlap between entities."""
        
        # collision tree is in local coords, and is cached
        if self._collision_tree is None:
            
            self._collision_tree = octree.Octree(worldSize=2.0*self.bounding_radius,
                                                 worldCentre=self.centre_of_diffusion,
                                                 max_objects_per_node=8)
            for atom in self.atoms:
                self._collision_tree.insertObject(atom, atom.position)
                
            ## now do a little housekeeping on the collision tree
            ## No need to have empty nodes for a start
            #self._collision_tree.prune()
            
            ## also don't need any nodes which don't define surface parts of
            ## the entity. Cull them by checking for any nodes where all
            ## neighbours are non-empty. Since we just pruned all empty nodes
            ## from the tree the quickest way to do this is to mark all nodes
            ## with exactly 27 items in the neighbours list for deletion.
            ## NB: This is only really worthwhile for very large proteins
            #ctr = 0
            #for node in self._collision_tree.all_nodes:
                #if len(node.colleagues) == 27:
                    #ctr += 1
                    #node.delete_me = True
            
            #data_content = len(self._collision_tree.master_list)
            #num_nodes = len(self._collision_tree.all_nodes)

            ## the above loop marked all interior nodes with 'delete_me'. Now
            ## delete them.
            #for node in self._collision_tree.all_nodes[:]:

                ## if the delete_me attribute is set, delete node
                #try:
                    #if node.delete_me:
                        #self._collision_tree.deleteNode(node)
                #except AttributeError: # no delete_me attribute
                    #pass
                #except ValueError: # no such node in tree (already deleted)
                    #pass 
                
            ## tree should now only contain non-empty nodes of collections of
            ## surface atoms
            #self._collision_tree.prune()
            #print "%d of %d nodes deleted; content before: %d after: %d" %(ctr, num_nodes, data_content, len(self._collision_tree.master_list))

        return self._collision_tree
        
    # utility function for coordinate transformations
    def convert_local_coordinate(self, pt):
        """Convert a point to Universe coordinate frame."""
        return local_to_universal(pt,
                                  ref_pt_local=self.centre_of_diffusion,
                                  ref_pt_universe=self.xyz,
                                  rotation=self.rotation)

    @property
    def diffusion_centre_universe_coordinates(self):
        """Return the centre of diffusion of the entity in Universe coordinate
        frame."""
        
        # by definition this is the xyz position w.r.t the universe coordinate frame
        return self.xyz
    
    @property
    def vis_points(self):
        """Return representative points for this diffusing entity, in Universe
        coordinates.
        
        The vis points are intended to be the set of points which represent
        this entity within a trajectory file. In the simplest case they are
        just the atoms of the molecule; for more coarse representations
        (blob-like proteins) we will need some reduced set of vis_points --
        perhaps a set of intersecting spheres which more-or-less encompass the
        blob volume."""

        # create list of visualisation points
        pts = [self.convert_local_coordinate(wrapped_atom.obj.position) 
               for wrapped_atom in self.collision_tree.master_list]
        
        # if pts is an empty list then just use the position of the entity
        if len(pts) == 0:
            # this should just return the centre of diffusion of the molecule
            # in the Universe coordinate frame -- which is the same as the xyz
            pts = [self.xyz]

        return pts
        
    @property
    def raw_vis_points(self):
        return [w.obj for w in self.collision_tree.master_list]
    
    @property
    def xyz(self):
        """Get xyz position of this diffusing entity with respect to Universe
        origin."""
        return self._xyz

    def set_xyz(self, new_xyz):
        """Set the position of this diffusing entity with respect to Universe
        origin."""
        self._xyz = Vector(new_xyz)
        return

    @property
    def rotation(self):
        """Get rotation of this diffusing entity with respect to Universe
        origin (as Quaternion object)."""
        return self._rotation

    def set_rotation(self, new_rotation):
        """Set rotation of this diffusing entity with respect to Universe
        origin (as Quaternion object)."""
        from copy import copy
        assert(isinstance(new_rotation, Quaternion))
        self._rotation = copy(new_rotation)
        return

    def apply_movement(self, translation=None, rotation=None):
        """Apply given translation and rotation to diffusing entity."""

        # sanity check the input parms
        if translation is None and rotation is None:
            raise ValueError

        # stash previous state in case we need to 'undo'
        self._undo_translation = self._xyz
        self._undo_rotation = self._rotation
        
        # only apply translation / rotation if they've actually been passed
        # in! (for example, for test purposes we may turn off rotational
        # diffusion)
        if translation is not None:
            self._xyz += translation
        if rotation is not None:
            self._rotation = rotation * self._rotation 

        # delete the cached universe coordinates
        try:
            del(self._cached_centre_universe_coords)
        except AttributeError:
            pass
        try:
            del(self._universe_atoms)
        except AttributeError:
            pass

        # clear the cached collision tree positions
        if self._collision_tree is not None:
            self._collision_tree.clear_universal_position_cache()
        
        return

    def undo_last_movement(self):
        """Undo the last movement."""
        
        self._xyz = self._undo_translation
        self._rotation = self._undo_rotation
        
        # delete the cached universe coordinates
        try:
            del(self._cached_centre_universe_coords)
        except AttributeError:
            pass
        try:
            del(self._universe_atoms)
        except AttributeError:
            pass

        # clear the cached collision tree positions
        if self._collision_tree is not None:
            self._collision_tree.clear_universal_position_cache()
        
        return
    
    @property
    def universe_atoms(self):
        """Return a list of all atom positions in universe coordinates."""
        
        try:
            return self._universe_atoms
        except AttributeError:
            uni_atoms = [Atom(self.convert_local_coordinate(at.position), 
                              at.radius) for at in self.atoms]
            self._universe_atoms = uni_atoms
        
        return self._universe_atoms

    @staticmethod
    def is_overlapping(entity1, entity2):
        """Return true if there is any overlap between entities 1 and 2."""
        
        # distance between two points
        r = (entity1.diffusion_centre_universe_coordinates - \
             entity2.diffusion_centre_universe_coordinates).length()
        
        # if this condition is met then the entities are not clashing at all
        if (r > (entity1.bounding_radius + entity2.bounding_radius)):
            return False
        
        # if this condition is met then the entities are seriously clashing
        if r < (entity1.guaranteed_clash_radius + \
                entity2.guaranteed_clash_radius):
            return True

        # otherwise we need to do something a little more subtle to find
        # whether the entities are indeed overlapping
        tree1 = entity1.collision_tree
        tree2 = entity2.collision_tree
                        
        top_level_clashes = clash_list(entity1, entity2, tree1, tree2)
        real_clashes = []
        for tlc in top_level_clashes:
            some_real_clashes = clash_list(entity2, entity1, tlc, tree1)
            for src in some_real_clashes:
                real_clashes.append( (src, tlc) )

        # now real_clashes contains a list of octnodes which overlap
        for node1, node2 in real_clashes:
            for atom1 in [n1.obj for n1 in node1._data]:
                for atom2 in [n2.obj for n2 in node2._data]:
                    
                    # don't forget to use universe coords :-)
                    p1 = entity1.convert_local_coordinate(atom1.position)
                    p2 = entity2.convert_local_coordinate(atom2.position)
                    
                    overlap = atom1.radius + atom2.radius - (p1 - p2).length()
                    if overlap > 0:
                        # ok *genuinely* these do overlap
                        return True
            
        # if get here then cannot be any overlapping!
        return False
    
    @staticmethod
    def get_overlaps(entity1, entity2, max_overlaps=None):
        """Return list of paired contacts between two entities."""

        # distance between two points
        r = (entity1.diffusion_centre_universe_coordinates - \
             entity2.diffusion_centre_universe_coordinates).length()
        
        # if this condition is met then the entities are not clashing at all
        if (r > (entity1.bounding_radius + entity2.bounding_radius)):
            return []
        
        # if this condition is met then the entities are seriously clashing
        if r < (entity1.guaranteed_clash_radius + \
                entity2.guaranteed_clash_radius):
            return [(entity1.diffusion_centre_universe_coordinates,
                     entity2.diffusion_centre_universe_coordinates,
                     entity1.guaranteed_clash_radius + \
                     entity2.guaranteed_clash_radius - r)]

        # otherwise we need to do something a little more subtle to find
        # whether the entities are indeed overlapping, and if so then where
        tree1 = entity1.collision_tree
        tree2 = entity2.collision_tree
                        
        top_level_clashes = clash_list(entity1, entity2, tree1, tree2)
        real_clashes = []
        for tlc in top_level_clashes:
            some_real_clashes = clash_list(entity2, entity1, tlc, tree1)
            for src in some_real_clashes:
                real_clashes.append( (src, tlc) )

        # now real_clashes contains a list of octnodes which overlap
        final_clash_list = []
        for node1, node2 in real_clashes:
            for atom1 in [n1.obj for n1 in node1._data]:
                for atom2 in [n2.obj for n2 in node2._data]:
                    
                    # don't forget to use universe coords :-)
                    p1 = entity1.convert_local_coordinate(atom1.position)
                    p2 = entity2.convert_local_coordinate(atom2.position)
                    
                    overlap = atom1.radius + atom2.radius - (p1 - p2).length()
                    if overlap > 0:
                        final_clash_list.append( (p1, p2, overlap) )
            
        return final_clash_list
    
    @staticmethod
    def get_overlaps_slow(entity, other, max_overlaps=None):
        """Return list of paired contacts between two entities."""

        # for quick collision check, the minimum separation between the two
        # objects is the sum of their bounding radii
        min_separation = entity.bounding_radius + other.bounding_radius
        
        # if the distance is greater than the min_separation, then it's not
        # possible for the two entities to overlap, so don't bother with an
        # atomic-level check
        r = entity.diffusion_centre_universe_coordinates - \
            other.diffusion_centre_universe_coordinates
        if (r.length() > min_separation):
            return []
        
        # ok so the bounding radii overlap, possibly there are atomic clashes,
        # possibly not. Since entities can be made up of a large number of
        # atoms, create an octree to hold the atoms of both entities and cross
        # check.
        overlaps = []
        import octree
        centre_for_octree = (entity.diffusion_centre_universe_coordinates + \
                             other.diffusion_centre_universe_coordinates) / 2.0
        root = octree.Octree(worldSize=2.0*min_separation,
                             worldCentre=centre_for_octree, 
                             max_objects_per_node=8)
        
        # entity 1 atoms into octree
        for atom in entity.universe_atoms:
            pos = atom.position
            #if (pos - other.diffusion_centre_universe_coordinates).length() < atom.radius + other.bounding_radius:
            wrapped = root.insertObject(atom, pos)
            wrapped.entity = entity
        
        # entity 2 atoms into octree
        for atom in other.universe_atoms:
            pos = atom.position
            #if (pos - entity.diffusion_centre_universe_coordinates).length() < atom.radius + entity.bounding_radius:
            wrapped = root.insertObject(atom, pos)
            wrapped.entity = other

        # optimize the tree by pruning empty nodes
        root.prune()
            
        # now iterate over all atoms and check for overlaps
        for wrapped_atom in root.master_list:
            
            atom = wrapped_atom.obj
            node = root.find_leaf_at_position(wrapped_atom.position)
            radius = atom.radius
            while (node.edge_length < radius and node.parent is not None):
                node = node.parent
            
            # get all atoms in neighbourhood which belong to the other entity
            other_atoms = []
            for neighbour in node.colleagues:
                other_atoms.extend(neighbour._data)
            other_atoms = [i for i in other_atoms 
                           if i.entity != wrapped_atom.entity]
            
            for other_atom in other_atoms:
                
                dist = other_atom.obj.radius + atom.radius
                r = (other_atom.position - wrapped_atom.position).length()
                if r < dist:
                    overlaps.append( (wrapped_atom.position, other_atom.position) )

        return overlaps
    
    @property
    def bounding_radius(self):
        """Return the bounding radius of this diffusing entity."""
        #return self.mesh.bounding_radius
        return max([(atom.position - self.centre_of_diffusion).length() + atom.radius
                    for atom in self.atoms] )
        
    @property
    def universe_bounding_box(self):
        """Return the bounding box (lower left, upper right corners).
        
        Based on the maximum radius of this molecule."""
        
        r = self.bounding_radius        
        upper_right = self.mesh._xyz + Vector(r,r,r)
        lower_left = self.mesh._xyz - Vector(r,r,r)
        return lower_left, upper_right

    @staticmethod
    def readHydroproResults(hydropro_results):

        # read the hydropro results file into an array
        h = open(hydropro_results, 'r')
        hydro = h.readlines()
        h.close()

        # lines 44-46 (inclusive) contain the centre of diffusion in cm
        # mult by 1e8 to convert cm to Angstroms
        cod = [string.atof(hydro[line].split()[4]) * 1e8 for line in [43,44,45]]
        centre_of_diffusion = Vector(cod)

        # lines 52-54,57-59 (inclusive) contain the 6x6 diffusion tensor
        diffusion_tensor = hydro[51:54]
        diffusion_tensor.extend(hydro[56:59])
        diffusion_tensor = [ [string.atof(x) for x in line.split()] \
                             for line in diffusion_tensor ]

        diff_tensor_as_matrix = matrix(diffusion_tensor)

        # convert units for the 3x3 translation component
        for i in range(3):
            for j in range(3):
                # units in hydropro are cm2/s -- we want Angstroms2/ns
                # which works out to a multiplication by 1e7:
                # (1e-2)cm*(1e-2)cm*(1e10)Angs*(1e10)Angs*(1e-9)ns
                diff_tensor_as_matrix[i,j] *= 1e7 

        # convert units for the translation-rotation coupling 3x3 matrices
        for i in range(3):
            for j in range(3,6):
                # units in hydropro are cm.rad/s we want Angstroms.rad/ns
                # which works out to multiplication by 1e-1
                diff_tensor_as_matrix[i,j] *= (1e-1)
                diff_tensor_as_matrix[j,i] *= (1e-1)

        # convert units for the rotation component 
        # from rad^2 per s to rad^2 per ns
        for i in range(3,6):
            for j in range(3,6):
                diff_tensor_as_matrix[i,j] *= 1e-9
                
        return diff_tensor_as_matrix, centre_of_diffusion    
    
    def set_hydropro_diffusion_tensor(self, hydropro_filename):
        """Extract and store diffusion characteristics from a hydropro results
        file.
        
        Hydropro is Garcia de la Torre's program for estimating the diffusion
        tensor of molecule from a PDB file."""
        
        diff, centre = DiffusingEntity.readHydroproResults(hydropro_filename)
        self.diffusion_tensor = diff
        self.set_centre_of_diffusion(centre)
        self.sqrt_diffusion_tensor = self.sqrt_matrix(diff)
        return

    def set_centre_of_diffusion(self, centre):
        """Set the centre of diffusion of the entity.
        
        All rotation of the entity occurs around this point. Also the xyz
        vector (defining the position of the diffusing entity coordinate frame
        w.r.t universe coordinate frame) is the vector from the universe
        origin to this point."""
        
        self._centre_of_diffusion = centre
        return
        
    def set_isotropic_diffusion(self, Dtr, Drot):
        """Set an isotropic diffusion tensor using the given coefficients.
        
        Dtr should be in Angstroms squared per nanosecond.
        
        Is a matrix which looks like:
        
        [[ Dt,  0.,  0., 0., 0., 0.],
         [ 0.,  Dt,  0., 0., 0., 0.],
         [ 0.,  0.,  Dt, 0., 0., 0.],
         [ 0.,  0.,  0., Dr, 0., 0.],
         [ 0.,  0.,  0., 0., Dr, 0.],
         [ 0.,  0.,  0., 0., 0., Dr]]
        """

        tensor = eye(DiffusingEntity.dimensions)

        assert(DiffusingEntity.dimensions == 6)
        D = array([Dtr, Dtr, Dtr, Drot, Drot, Drot])
        
        self.diffusion_tensor = mat(D * tensor)
        self.sqrt_diffusion_tensor = self.sqrt_matrix(mat(D * tensor))
        return
    
    @staticmethod
    def sqrt_matrix(matrix):
        """Returns the 'square root' of a matrix.
        
        i.e. solves for S in Q=S.S, given Q by diagonalization."""
        
        w, v = linalg.eig(matrix)
        w = [math.sqrt(eigval) for eigval in w]
        ww = mat(eye(matrix.shape[0]))
        for i,eigval in zip(range(matrix.shape[0]), w):
            ww[i,i] = eigval
        return mat(v) * (ww * mat(v).I)        


    def set_atoms_from_pdb(self, pdb_filename):
        """Warning: deprecated function. Use set_charges_from_pqr method
        instead.
        
        Set the atoms of the entity from a PDB file. Note radii all set to
        zero, since they are not defined in a PDB file!
        
        NB: The 'atoms' don't have to be true atoms -- can be large blob-like
        things to more coarsely approximate the diffusing entity. These atoms
        can be used for surface meshing and for trajectory visualisation. Also
        they may be used as the basis of collision-detection for volume
        exclusion (unless a surface collision algorithm is implemented
        instead).
        
        In practice the atoms should be set from a PQR file using
        set_charges_from_pqr with the set_atoms_from_pqr flag set the True.
        """
        
        # load the PDB file into this object
        # we should use biopython instead of MMTK (which is frankly a bit unstable...)
        import pdbtools
        self.atoms = pdbtools.get_atoms_from_pdb(pdb_filename)
    
        return
    
    def set_charges_from_pqr(self, pqr_filename, set_atoms_from_pqr=True):
        """Set the charges within the diffusing entity from a PQR file.
        
        NB: set_atoms_from_pqr controls whether to use the Charge objects as
        atoms in the atom_list container. This has implications for surface
        mesh generation and for what is used for visualisation points. For
        atomically detailed representations, this can be set to True without
        any problem. However if the representation is more 'blob-like' then
        the charges may not reside at the atom locations, so charges and atoms
        are not interchangeable."""

        # DAVE HACK
        self.pqr_filename = pqr_filename
        
        # get the PQR data (Electrostatic charge assignments) the pqrtools
        # 'module' provides a parse for PQR files, based on underlying class
        # within Scientific library (Scientific.IO.PDB). (Could use BioPython
        # instead?)
        from pqrtools import getCharges
        self.charges = getCharges(pqr_filename)
        
        if set_atoms_from_pqr:
            self.atoms = [Atom(c.position, c.radius) for c in self.charges]

        # remove zero-charges from charge list
        self.charges = [c for c in self.charges if c.charge != 0.0]
        
        return

    def create_effective_charges(self, filename):
        """Calculates and stores effective charges to reproduce the potential
        pattern outside the diffusing entity."""

        from pqrtools import getCharges
        self.__ecm = getCharges(filename + ".ecm")

        from opendx_utils import OpenDXGrid
        self.opendx = OpenDXGrid.get_ecm_model(filename)

        return

    @property
    def ecm(self):
        
        ecm_copy = self.__ecm[:]
        for ecm in ecm_copy:
            ecm._position = local_to_universal(ecm.position, 
                                               self.centre_of_diffusion, 
                                               self.xyz, 
                                               self.rotation)
        return ecm_copy
    
    def universe_charges_and_positions(self):
        """Return lists of charges and charge positions."""
    
        charges = []
        charge_positions = []
    
        # get the electrostatic mesh object from the diffusing entity
        mesh = self.mesh
    
        # this will get the charges (in the Universe co-ordinate frame)
        from geometry import local_to_universal
        from BEM import VectorPoint
        cvt = lambda x: local_to_universal(x, 
                                           ref_pt_local=self.centre_of_diffusion,
                                           ref_pt_universe=self.xyz, 
                                           rotation=self.rotation)
        for charge_obj in mesh.charges:
            position = cvt(charge_obj.position)
            charge = charge_obj.charge
            charges.append(charge)
            charge_positions.append( VectorPoint(position) )
        
        return charges, charge_positions
    
    def calculateElectrostaticForce(self, kappa, Dext):
        return self.mesh.calculate_force(kappa, Dext)

    def calculateElectrostaticMoment(self, kappa, Dext):
        return self.mesh.calculate_moment(kappa, Dext)

    def ecm_energy_force_moment(self, others):
        
        total_energy = 0.0
        total_force = Vector(0.0,0.0,0.0)
        total_moment = Vector(0.0,0.0,0.0)
        for other_entity in others:
            if self is other_entity: continue # don't do self-self charge/field interaction
            if ((self.xyz - other_entity.xyz).length() < 100.0):
            
                rot = other_entity.rotation.inverse()
                for ch in self.ecm:
                    xyz_in_opendx_grid = apply_quaternion_to_vector(rot, 
                                                                    ch.position - other_entity.xyz)
                    potential, field = self.opendx.get_field_at(xyz_in_opendx_grid)
                    force = field*ch.charge
                
                    total_energy += potential*ch.charge
                    total_force += force
                    total_moment += (ch.position - self.xyz).cross(force)
        
        force = apply_quaternion_to_vector(self.rotation.inverse(), total_force)
        mom = apply_quaternion_to_vector(self.rotation.inverse(), total_moment)
        
        return total_energy, force, mom
    
    def getForcedMove(self, t, kappa, Dext, others):

        #F = self.calculateElectrostaticForce(kappa, Dext)
        #M = self.calculateElectrostaticMoment(kappa, Dext)
        e,F,M = self.ecm_energy_force_moment(others)

        # put F and M together into one vector
        FM = [F.x(), F.y(), F.z(), M.x(), M.y(), M.z()]
        
        # NB: at this point the forces are in units of kT per Angstrom (where
        # kT is boltzmann constant multiplied by temperature)

        # Return Diffusion tensor x force/torque
        deltas = t * self.diffusion_tensor * matrix(FM).T
        
        # deltas is dx,dy,dz (translations) and drx,dry,drz (rotations)
        [dx, dy, dz, drx, dry, drz] = [float(v) for v in deltas]
        trans_local = Vector(dx, dy, dz)
        trans_universe = apply_quaternion_to_vector(self.rotation, trans_local)
        
        # use unbiased rotation operator
        # see Beard and Schlick Biophys Journal, 2003
        rot_vector = Vector(drx, dry, drz)
        if rot_vector.length() > 0.0:
            rot_vector_uni = apply_quaternion_to_vector(self.rotation, rot_vector)
            rot = Rotation(rot_vector_uni.normal(), rot_vector_uni.length()).asQuaternion()
        else:
            rot = None
            
        return trans_universe, rot, e

    def getDiffusionMove(self, t):
        """Return diffusional translation vector and rotation quaternion,
        given a timestep (in nanoseconds)."""

        # this diffusing entity has a diffusion tensor; use it to iterate the
        # position via Ermak and McCammon algorithm.
        rand = [random.gauss(0.0, math.sqrt(2.0*t))
                for i in range(DiffusingEntity.dimensions)]

        diffusion = self.sqrt_diffusion_tensor * matrix(rand).T
        [dx, dy, dz, drx, dry, drz] = [float(v) for v in diffusion]
        trans_local = Vector(dx, dy, dz)
        trans_universe = apply_quaternion_to_vector(self.rotation, trans_local)
        
        # use unbiased rotation operator
        # see Beard and Schlick Biophys Journal, 2003
        rot_vector = Vector(drx, dry, drz)
        if rot_vector.length() > 0.0:
            rot_vector_uni = apply_quaternion_to_vector(self.rotation, rot_vector)
            rot = Rotation(rot_vector_uni.normal(), rot_vector_uni.length()).asQuaternion()
        else:
            rot = None
        
        return trans_universe, rot
    
def old_code():

    # the surface_mesh passed in is for electrostatics stuff
    # it also keeps track of where this diffusing_entity is in relation
    # to the rest of the Universe.
    if surface_mesh is not None:
        self.mesh = surface_mesh
    else:
        self.mesh = MSMS_Mesh(self.atom_list)

    # this imports existing solutions to the electrostatic mesh
    # to reduce gmres iterations to convergence
    pre_calculated_pots = []
    try:
        pre_calc_filename = filename + ".surfpot"
        pot_file = open(pre_calc_filename, 'r')
        for line in pot_file:
            if line.strip() != "":
                solution = (f,h) = [float(x) for x in line.split()]
                pre_calculated_pots.append(solution)
    except TypeError:
        pass
    except IOError:
        pass
    
    if len(pre_calculated_pots) == len(self.mesh.vertices):
        for solpt, vert in zip(pre_calculated_pots, self.mesh.vertices):
            vert.f, vert.h = solpt

    #self.mesh.resetOrigin(diff_centre)
    self.mesh.setRotationFromUniverse(rot.asQuaternion())
    self.mesh.setXYZFromUniverse(xyz)

    # HACK: figure out whether charges have already been added to the mesh!
    # For MSMS they will have; for ObjMesh they won't.
    if len(self.mesh.charges) == 0:

        # add the charges to the mesh
        for atom in self.atom_list:
            if atom.properties['charge'] != 0.0:
                self.mesh.addCharge(atom.position, atom.properties['charge'] * 1.6e-19) # note conversion to Coulombs

    # move the centre of diffusion to the origin of the molecule
    # (simplifies the geometry)
    ## HACK!! ## 
    #for atom in self.atom_list:
    #    atom.position -= diff_centre
    
           # assign charge and radius to atom properties
    for r in struct.residues:

        DiffusingEntity.res_ctr += 1

        # write the residue as PDB to dummy fileobject
        dummy = StringIO()
        dummy_pdbfile = SpecialPDBFile(dummy)
        r.writeToFile(dummy_pdbfile)
        pdb_lines = dummy.getvalue().split('\n')
        dummy_pdbfile.close()

        for atom, pdb_line in zip(r, pdb_lines):

            DiffusingEntity.atom_ctr += 1

            # implicit unit is angstroms
            ##atom.position *= 1e-10

            self.atom_list.append(atom)
            num = atom['serial_number']
            atom.properties['charge'], atom.properties['radius'] = pqr_data[num]
            atom.properties['occupancy'], atom.properties['temperature_factor'] = pqr_data[num]

            line = pdb_line.split()
            line[1] = "%d" %(DiffusingEntity.atom_ctr)
            line[4] = "%d" %(DiffusingEntity.res_ctr)
            line[8] = "%f" %(pqr_data[num][0])
            line[9] = "%f" %(pqr_data[num][1])
            atom.pdb_line = " ".join(line) 

def clash(entity1, entity2, node1, node2):
    """Utility function to reveal if two spherical nodes in octree
    overlap."""
    
    # 4.0 angstroms is the approximate radius of 2 carbon atoms --
    # note that the octree assumes point entities
    r = node1.radius + node2.radius + 4.0 
    dist = (node1.uni_position(entity1) - node2.uni_position(entity2)).length()
    return dist < r

def clash_list(entity1, entity2, node1, node2):

    clashes = []
    check_list = [node2]
    
    # loop until no more nodes to check
    while len(check_list) > 0:
        check_me = check_list.pop()
        
        # if this node clashes, then if it's a leaf node then record
        # it, otherwise check all child nodes of this node
        if clash(entity1, entity2, node1, check_me):
            if check_me.isLeafNode:
                clashes.append(check_me)
            else:
                check_list.extend(check_me.children)
        else:
            # check_me does not clash with node1, so ignore it
            pass
                
    return clashes
