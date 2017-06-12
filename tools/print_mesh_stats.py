#!/usr/bin/python
# -*- coding: utf-8 -*-
import pybeep

if __name__=="__main__":

    import sys
    from math import sqrt
    from random import uniform
    import constants

    mesh_filename = sys.argv[1]
    m = pybeep.Mesh(mesh_filename)
    total_area = sum([n.bezier_area() for n in m.node_patches])
    area_per_vertex = total_area / m.num_node_patches

    print("%s np=%d area=%f rad=%f vol=%f apv=%f" %(mesh_filename, m.num_node_patches, total_area, m.get_radius(), m.calculate_volume() / 6., area_per_vertex))

