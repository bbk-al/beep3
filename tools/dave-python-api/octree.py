#
# octree.py
#
# Octree is a very useful structure for hierarchical storage of objects by xyz
# co-ordinates.
#
# This forms the basis for some Fast Multipole Method stuff (see BEM
# electrostatics), as well as fast-collision detection for diffusing entities
# within the 'Universe' of a Brownian Dynamics simulation.
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

from Scientific.Geometry import Vector
import math

sq_root_three = math.sqrt(3)

class OctNode(object):

    """Node class for Octree."""
    
    def __init__(self, parent, position, edge_length, data=[]):
        """Constructs an Octnode (octree element)."""
        
        # OctNode cubes; position is centre of cube with edges of edge_length.
        self.position = position
        self.edge_length = edge_length
        self.parent = parent
        self.radius = math.sqrt(3.0*0.25*edge_length*edge_length)
        
        # set caches for direction list etc.
        self.clear_lists()
        
        if parent is None:
            self.root = self
            self.level = 1
        else:
            self.root = parent.root
            self.level = parent.level + 1
        
        # refer to parent for maximum_objects_per_node; for top level node this
        # will be set in the Octree c'tor.
        try:
            self.maximum_objects_per_node = self.parent.maximum_objects_per_node
        except AttributeError:
            pass

        # This is where we will contain the objects belonging to this node.
        self._data = data[:]
        
        # no children by default
        self._children = []

        # this is used for collision checking in universe coordinates
        # TODO: this should be part of a specific collision tree subclass
        self._uni_pos = None

    def __getattr__(self,name):
        
        if name == 'mpole':
            import numpy
            nterms = 9
            return(numpy.zeros((nterms+1, nterms+1), 'complex'))
        
    def clear_lists(self):
        """Reset the cached interaction/direction/colleague lists."""
        
        self._colleagues = None
        self._interaction_list = None
        self._neighbourhood = None
        self._outer_edges = None

        return

    def _calculate_interaction_lists(self):    
        
        # clear lists
        self.clear_lists()

        # force recalc
        self.calc_direction_lists()
        #self.calc_sub_direction_lists()

        # recurse
        for child in self.children:
            child._calculate_interaction_lists()
        
        return
    
    @property
    def _decendents(self):
        """Return list of all decendents of this node (incl. self).
        
        Recursively decends tree."""
        
        all_decendents = [self]
        
        for child in self.children:
            all_decendents.extend(child._decendents)
        
        return all_decendents
    
    @property
    def isLeafNode(self):
        """Return true if is a leaf node (i.e. no children)."""
        return len(self.children) == 0
        
    def _insertObject(self, wrapped_object):
        """Insert an object into the octree.
        
        This should be a private function -- only insert objects into the tree
        by calling the Octree.insertObject function, which wraps this."""
        
        # add the object to the data contained within this node.
        self._data.append(wrapped_object)
        
        # if we now have too many items in this node then subdivide
        if (self.isLeafNode and 
            len(self._data) > self.maximum_objects_per_node and 
            self.root.max_depth < self.root._max_allowed_depth): 
            self._subdivide() # subdivide() will also pass the obj to children
            
        elif not self.isLeafNode:
            
            # if this is not a leaf node, then figure out which child node
            # of this node the data also belongs to.
            next_level_node = self._findChild(wrapped_object.position)
            next_level_node._insertObject(wrapped_object)
            
        return

    def _prune(self):
        """Recursively cull empty nodes of the tree (i.e. nodes which contain
        no data)."""
        
        if len(self._data) == 0:
            self._destroy()
        else:
            for c in self.children:
                c._prune()
        
        return
    
    def _destroy_data(self, data_list):
        """Removes data from the tree, upwardly recursive."""
        
        # recurse upwards through the tree
        if self.parent is not None:
            self.parent._destroy_data(data_list)

        for item in data_list:
            self._data.remove(item)
        
        return
    
    def _destroy(self):
        """Destroy this node (and all it's decendents).
        
        Remove all children and self from the master leaf list."""
        
        # detroy all children
        for child in self.children:
            child._destroy()
        
        # remove from the 'children' list of the parent
        self.parent._children.remove(self)

        return

    def _subdivide(self):
        """Subdivide this node into 8 children."""

        # can't subdivide an already subdivided node!
        assert(self.isLeafNode)
            
        # it's a leaf node, time to subdivide it.
        new_edge_length = self.edge_length / 2.0

        # this is how far the new centres will be from the current
        # node centre after subdivision... (in each dimension...)
        offset = self.edge_length / 4.0
        pos = self.position
        
        # create 8 children
        self.wsd = OctNode(self, pos + Vector(-offset, -offset, -offset),
                           new_edge_length, [])
        self.wsu = OctNode(self, pos + Vector(-offset, -offset, +offset),
                           new_edge_length, [])
        self.wnd = OctNode(self, pos + Vector(-offset, +offset, -offset),
                           new_edge_length, [])
        self.wnu = OctNode(self, pos + Vector(-offset, +offset, +offset),
                           new_edge_length, [])
        self.esd = OctNode(self, pos + Vector(+offset, -offset, -offset),
                           new_edge_length, [])
        self.esu = OctNode(self, pos + Vector(+offset, -offset, +offset),
                           new_edge_length, [])
        self.end = OctNode(self, pos + Vector(+offset, +offset, -offset),
                           new_edge_length, [])
        self.enu = OctNode(self, pos + Vector(+offset, +offset, +offset),
                           new_edge_length, [])

        # list of children according to Huang ordering
        self._children = [self.end, self.wnd, self.wsd, self.esd,
                          self.enu, self.wnu, self.wsu, self.esu]
        
        # set Huang octant id (for use with Fortran FMM code)
        for idx, child in enumerate(self._children):
            child.Huang_octant_id = idx+1

        # update the maximum recursion depth counter in the tree
        self.root._max_depth = max(self.root._max_depth, self.level + 1)

        # deal out the data objects amongst the new children
        for wrapped_object in self._data:
            next_level_node = self._findChild(wrapped_object.position)
            next_level_node._insertObject(wrapped_object)
           
        return
 
    def _subdivide_to_level(self, max_depth):
        """Ensure subdivision is at max_depth level (ignore if empty node)."""
        
        # if this is a leaf node and not at max depth, and contains some stuff
        if self.isLeafNode and self.level < max_depth:# and len(self._data) > 0:
            self._subdivide()

        # pass this call onto the children (if any exist)
        for child in self.children:
            child._subdivide_to_level(max_depth)
        
        return
    
    @property
    def children(self):
        """Return the children of this box."""
        return self._children

    @property
    def colleagues(self):
        """Return the colleagues of this box.
        
        Colleagues are boxes at the same refinement level which share a
        boundary point."""

        if self._colleagues is None:
            
            # if this box is the root node then it is it's only colleague
            if self.parent is None: return [self]
            
            # colleagues are all the siblings of this box...
            colleagues = self.parent.children[:]
            #assert(self in colleagues)
            
            # ...plus the children of parents colleagues which adjoin
            uncles_aunts = [c for c in self.parent.colleagues 
                            if c is not self.parent]
    
            for uncle_aunt in uncles_aunts:
                for cousin in uncle_aunt.children:
                    deltas = (cousin.position - self.position) \
                                    / self.edge_length
                    [dx, dy, dz] = [abs(int(round(dd))) for dd in deltas]
                    if dx <= 1 and dy <= 1 and dz <= 1:
                        colleagues.append(cousin)
            
            self._colleagues = colleagues
            assert( len(self._colleagues) <= 27 )
            
        return self._colleagues[:]

    def isWellSeparatedFrom(self, other):
        """Returns whether or not the other node is well separated from this
        node."""

        # ensure we always compare at same level in the tree
        assert( self.level == other.level )
        
        # check that we're at the same refinement level
        #if self.level != other.level:
            # if not at the same level then assume not well separated
        #    return False
        
        return other not in self.colleagues
    
    @property
    def interaction_list(self):
        """Return the interaction list of this cube."""

        # cache once calculated
        if self._interaction_list is None:
            
            interaction_list = []
            
            if self.root is not self:
                uncles_aunts = [c for c in self.parent.colleagues 
                                if c is not self.parent]

                for aunt_uncle in uncles_aunts:
                    for cousin in aunt_uncle.children:
                        if self.isWellSeparatedFrom(cousin):
                            interaction_list.append(cousin)
            
            self._interaction_list = interaction_list

        return self._interaction_list[:]

    def calc_direction_lists(self):
        """Precalc the up- down- north- south- east- west- lists."""

        # working copy of the interaction list
        interaction_list = self.interaction_list[:]

        buf = self.edge_length

        self.uplist = [box for box in interaction_list 
                       if box.position.z() > (self.position.z() + buf)]
        self.downlist = [box for box in interaction_list 
                         if box.position.z() < (self.position.z() - buf)]
        
        # remove the uplist and downlist from the interaction list
        # components left in play
        interaction_list = [box for box in interaction_list 
                            if (box not in self.uplist) and 
                            (box not in self.downlist)]
        
        self.northlist = [box for box in interaction_list 
                          if box.position.y() > (self.position.y() + buf)]
        self.southlist = [box for box in interaction_list 
                          if box.position.y() < (self.position.y() - buf)]
        
        # remove the north- and south- lists from the interaction list
        # components left in play
        interaction_list = [box for box in interaction_list 
                           if (box not in self.northlist) and 
                           (box not in self.southlist)]
        
        self.eastlist = [box for box in interaction_list 
                         if box.position.x() > (self.position.x() + buf)]
        self.westlist = [box for box in interaction_list 
                         if box.position.x() < (self.position.x() - buf)]

        # check that we've accounted for all the boxes in the interaction list
        interaction_list = [box for box in interaction_list 
                           if (box not in self.eastlist) and 
                           (box not in self.westlist)]
        assert(len(interaction_list) == 0) 
        
        return

    @property
    def neighbourhood(self):
        """The 3x3x3 cladding of child-sized cubes around this cube (octet)."""
        
        if self._neighbourhood is None:
            
            neighbourhood = []
            for neighbour in self.colleagues:
                for child in neighbour.children:
                    deltas = (child.position - self.position) \
                           / (child.edge_length / 2.0)
                    [dx, dy, dz] = [abs(int(delt)) for delt in deltas]
                    if dx <= 3 and dy <= 3 and dz <=3:
                        neighbourhood.append(child)
                        
            neighbourhood = [n for n in neighbourhood if n not in self.children]
            self._neighbourhood = neighbourhood
        return self._neighbourhood

    @property
    def outer_edges(self):
        
        if self._outer_edges is None:
            
            # all outer edge boxes
            outer_edges = []
            for neighbour in self.colleagues:
                for child in neighbour.children:
                    if child not in self.neighbourhood:
                        outer_edges.append(child)
            self._outer_edges = outer_edges
            
        return self._outer_edges

    def calc_sub_direction_lists(self):
        
        edges = self.outer_edges[:]
        buf = self.edge_length/2.0
        self.far_uplist = [outer for outer in edges
                           if outer.position.z() > \
                            (self.position.z() + self.edge_length)]
        self.near_uplist = [n for n in self.neighbourhood 
                            if n.position.z() > self.position.z()+buf]
        
        self.far_downlist = [outer for outer in edges
                             if outer.position.z() < \
                             (self.position.z() - self.edge_length)]
        self.near_downlist = [n for n in self.neighbourhood 
                              if n.position.z() < self.position.z()-buf]

        edges = [edge for edge in edges 
                 if (edge not in self.far_uplist) and 
                 (edge not in self.far_downlist)]
        self.far_northlist = [outer for outer in edges
                           if outer.position.y() > \
                            (self.position.y() + self.edge_length)]
        self.near_northlist = [n for n in self.neighbourhood 
                            if n.position.y() > self.position.y()+buf]
        
        self.far_southlist = [outer for outer in edges
                             if outer.position.y() < \
                             (self.position.y() - self.edge_length)]
        self.near_southlist = [n for n in self.neighbourhood 
                              if n.position.y() < self.position.y()-buf]

        edges = [edge for edge in edges 
                 if (edge not in self.far_northlist) and 
                 (edge not in self.far_southlist)]
        self.far_eastlist = [outer for outer in edges
                             if outer.position.x() > \
                            (self.position.x() + self.edge_length)]
        self.near_eastlist = [n for n in self.neighbourhood 
                              if n.position.x() > self.position.x()+buf]
        
        self.far_westlist = [outer for outer in edges
                             if outer.position.x() < \
                             (self.position.x() - self.edge_length)]
        self.near_westlist = [n for n in self.neighbourhood 
                              if n.position.x() < self.position.x()-buf]
        
        return

    def _findPosition(self, target_position):
        """Return the leaf node corresponding to position given."""
        
        # TODO: Can be optimized by storing a list of all leaf (bottom level)
        # nodes and going straight to the correct leaf node by doing a modulo
        # divide on the target position (requires that the tree is totally
        # uniform with equally sized leaf nodes spanning the whole volume --
        # which is the case for the current version of the code anyway.

        # at the moment we find the leaf node recursively. It's not *that*
        # slow.
        
        if self.isLeafNode: 
            return self
        else:
            child = self._findChild(target_position)
            return child._findPosition(target_position)

    def toSpherical(self, target_position):
        """Returns the target_position in spherical co-ords from centre of
        cube.

        Returns rho, theta, phi.
        """

        # get the cartesian vector from centre of cube to target position
        x,y,z = vec = target_position - self.position

        # delegate the calculation to the spherical.f f2py module for speed
        import spherical        
        return spherical.spherical_coordinates(x,y,z)
            
    def _findChild(self, target_position):
        """Find the child octant corresponding to the object target
        position."""

        # follow the definitions that 
        # +ve x is east, -ve x is west;
        # +ve y is north, -ve y is south;
        # +ve z is up, -ve z is down.

        # this can probably be optimized but I can't be bothered.
        
        if target_position.x() < self.position.x():
            # westerly
            if target_position.y() < self.position.y():
                # southerly
                if target_position.z() < self.position.z():
                    # down
                    return self.wsd
                else:
                    # up
                    return self.wsu
            else:
                # northerly                
                if target_position.z() < self.position.z():
                    # down
                    return self.wnd
                else:
                    # up
                    return self.wnu                
        else:
            # easterly
            if target_position.y() < self.position.y():
                # southerly
                if target_position.z() < self.position.z():
                    # down
                    return self.esd
                else:
                    # up
                    return self.esu
            else:
                # northerly                
                if target_position.z() < self.position.z():
                    # down
                    return self.end
                else:
                    # up
                    return self.enu                

        # shouldn't get here.
        raise Exception

    def clear_universal_position_cache(self):
        
        self._uni_pos = None
        for c in self.children:
            c.clear_universal_position_cache()
            
        return
    
    def uni_position(self, reference_entity):
        
        if self._uni_pos is None:
            self._uni_pos = reference_entity.convert_local_coordinate(self.position)
            
        return self._uni_pos
    
    
    
class OctreeDataObject(object):
    
    """Represents an object in the octree. This is an almost transparent
    read-only wrapping about the original object in that any get-attribute
    request on this wrapping passes straight through to the __getattribute__
    accessor on the object; i.e. it's as though you were directly accessing
    that object. There is no pass-through wrapping for setattr to avoid
    accidental rebinding of important attributes in the wrapped object when
    what is really intended is setting a temporary attribute on the wrapper
    object, as a workspace.
    
    The only exception is the position attribute which returns the position of
    this object within the Octree data structure, which for various reasons
    may not be the same as the position according to the object itself, e.g.
    due to different coordinate system between octree and object. TODO: put
    some example syntax here."""
    
    def __init__(self, wrapped_obj, position):
        """Wrap the object at the given position."""
        
        self._position = position
        self._obj = wrapped_obj
        
    @property
    def position(self):
        """Return the xyz position in the tree of this data object."""
        return self._position
        
    @property
    def obj(self):
        """Return the wrapped object."""
        return self._obj

    # all other attribute requests should be addressed directly to the wrapped
    # object
    def __getattr__(self, name):
        return self._obj.__getattribute__(name)
    
class Octree(OctNode):
    
    def __init__(self, 
                 worldSize=0.0,
                 worldCentre=Vector([0.0,0.0,0.0]), 
                 max_objects_per_node=100):
        """Constructor for an octree. All we need is a worldSize for starters.
        
        Add stuff to the tree with insertObject."""
        
        super(Octree, self).__init__(None, worldCentre, worldSize, [])
        self.maximum_objects_per_node = max_objects_per_node
        
        # keep track of how many levels there are in the tree in total
        self._max_depth = 1
        
        # don't allow more than a certain number of levels
        self._max_allowed_depth = 20
        
    def clear(self):
        """Clear the octree."""
        
        self._data = []
        self._children = []
        self.clear_lists()
        
        return

    @property
    def all_nodes(self):
        """All nodes in the tree are decendents of the root node."""
        
        return self._decendents
    
    @property
    def master_list(self):
        """Return copy of list of all data in the tree."""
        return self._data[:]

    @property
    def total_leaves(self):
        """Return the total number of leaf nodes."""
        return len(self.leaves)
    
    @property
    def max_depth(self):
        """Return the highest number of subdivisions in the tree."""
        return self._max_depth
        
    def deleteNode(self, node):
        """Remove a node (and it's children, and data) from the tree."""

        # root node is indelible!
        if node is self:
            return
        
        # check the node is actually in the tree
        if node not in self.all_nodes:
            raise ValueError
        
        # delete the node data, then the node
        node._destroy_data(node._data[:])
        node._destroy()
        
        return
    
    def insertObject(self, objData, position=None):
        """Insert an object into the tree at supplied position.
        
        If position is not specified, will check position attribute of
        objData.
        
        Returns the wrapped object which has been inserted into the tree."""

        # if no position supplied, try to find a position attribute
        if position is None:
                position = objData.position

        # check for this object falling outside the size of this octree
        for i, xyz in enumerate(position):
            if xyz > self.position[i] + (self.edge_length / 2.0): raise ValueError
            if xyz < self.position[i] - (self.edge_length / 2.0): raise ValueError
                
        # wrap the object, attach position to wrapped object
        wrapped_object = OctreeDataObject(objData, position)
            
        self._insertObject(wrapped_object)
        
        return wrapped_object
            
    @property
    def leaves(self):
        """Return a list of all leaf nodes of this tree."""
        return [n for n in self.all_nodes if n.isLeafNode]

    def optimize(self):
        """Homogenize (enforce constant level depth) then prune (remove empty
        nodes) then pre-calculate the interaction lists.
        
        This is called to slightly optimize the FMM methods."""
        
        self.homogenize()
        #self.prune()
        self.calculate_interaction_lists()
        
        return    
    
    def homogenize(self):
        """Homogenize the tree -- force equal refinement at lowest level."""
        self._subdivide_to_level(self.max_depth)
        
    def prune(self):
        """Remove all empty leaf nodes. 
        
        Only do this on a homegenized tree when you're sure you're not going
        to add any more objects (as they might want to be added to
        non-existant leaf nodes!)."""
        
        self._prune() # see implementation in OctNode
        return
    
    def calculate_interaction_lists(self):

        self.clear_lists()
        
        # call underlying node-level method which traverses tree
        self._calculate_interaction_lists()
        
        return

    def find_leaf_at_position(self, position):
        """Return the bottom level node corresponding to this position."""
        
        return self._findPosition(position)
    
# put some good tests here
if __name__ == "__main__":

    myTree = Octree(1.0)
    for i in range(10000):
        import random
        v = [random.random() - 0.5,random.random() - 0.5,random.random() - 0.5]
        myTree.insertObject(None, Vector(v))
        
    print myTree.max_depth
