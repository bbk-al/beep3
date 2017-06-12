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
# Uses C++ module _OCTREE.so
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

import _OCTREE
from _BEM import Vector
    
class Octree(_OCTREE.Octree, object):
    
    def __init__(self, 
                 worldSize=0.0,
                 worldCentre=Vector([0.0,0.0,0.0]),
                 num_levels=2):
        """Constructor for an octree. All we need is a worldSize for starters.
        
        Add stuff to the tree with insertObject."""
        
        super(Octree, self).__init__(worldCentre, worldSize, num_levels)
        return
    
# put some good tests here
if __name__ == "__main__":

    myTree = Octree(1.0)
    for i in range(100):
        import random
        v = [random.random() - 0.5,random.random() - 0.5,random.random() - 0.5]
        node = myTree.get_node_at_position(Vector(v))
        print node.get_centre()
        
    