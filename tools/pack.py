#!/usr/bin/python
# -*- coding: utf-8 -*-
from pybeep import Vector
from geometry import _rand_rot, _rand_xyz
import sys
from math import pi, ceil

def pack(crowdant_radius, volume_fraction, distance_limit, forbidden, minimum_separation):
    """pack a load of spherical crowders around a set of existing spherical objects, to a 
       specified distance_limit with specified total volume_fraction occluded (approx.)"""
    
    # centroid of the forbidden set
    centroid = Vector(0,0,0)
    for xx,rad in forbidden:
        centroid += xx
    centroid /= len(forbidden)
    
    # semi-axis limits -- centroid to furthest point in forbidden set
    dists = [(xx-centroid).length()+radx for (xx,radx) in forbidden]
    edge_length = max(dists)*2 + distance_limit*2 + crowdant_radius*4
    volume = edge_length ** 3
    
    # calculate number required to reach specified volume fraction
    individual_vol = (4. / 3.) * pi * crowdant_radius ** 3
    target_vol = volume * volume_fraction
    n = int(ceil(target_vol / individual_vol))
    
    # now loop 
    centres = forbidden[:]
    num_added = 0
    iterations_since_successful_addition = 0
    while (num_added < n and iterations_since_successful_addition < 10000):
        
        new_pt = _rand_xyz(edge_length) + centroid
        ok = True
        for xx,r in centres:
            if ((new_pt - xx).length() < (minimum_separation+crowdant_radius+r)):
                ok = False
                break
        if ok:
            iterations_since_successful_addition = 0
            num_added += 1
            centres.append((new_pt,crowdant_radius))
        else:
            iterations_since_successful_addition += 1

    ##rot = Quaternion(0,0,0,0)
    #for ctr,(pt,d) in enumerate(centres):
        ##rot = _rand_rot()
        #rot = Quaternion(0,0,0,0)
        #print "<instance instance_id=%d mesh_id=0 conformation=0 x=%f y=%f z=%f a=%f b=%f c=%f d=%f/>" %(ctr, pt.x, pt.y, pt.z, rot.x, rot.y, rot.z, rot.w)
                ##print "%f %f %f" %(pt.x, pt.y, pt.z)
    #print

    final_vol_fraction = num_added * individual_vol / volume
    print("Added %d crowders, reached volume fraction of %f" %(num_added, final_vol_fraction))
    
    centres = centres[len(forbidden):]
    
    print("Culling anything further than %f from the forbidden range" %(distance_limit))
    
    # utility function- returns true if the point/radius is further than
    # distance_limit from any object in forbidden (NB: could use lambda
    # function here)
    def cull_function(x,r):
        dists = [(x - fx).length()-fr-r for fx,fr in forbidden]
        return (min(dists) < distance_limit)

    # cull using list comprehension
    centres = [x for x,r in centres if cull_function(x,r)]
    
    print("%d crowders remain" %(len(centres)))
    
    return centres
