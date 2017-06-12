#!/usr/bin/python
# vim: set fileencoding=UTF8 :

#
# coalesce.py
#
# Utilities for cleaning up MSMS-based triangulations. Based on the algorithms
# contained in COALESCE by Sergio Aragon (see Aragon2003).
#
# Author: David Fallaize, University College London, 2009
# E-mail: drf33@cantab.net
#
from vector import *
from random import shuffle

def coalesce(vertex_list, triangle_defs, normals, minimum_length):

    # sanity check the input parameter
    #assert(short_fraction > 0 and short_fraction < 1)
    
    # figure out what is 'short'
    #average_length = average_side_length(vertex_list, triangle_defs)
    #too_short = average_length * 0.1
    too_short = minimum_length

    init_num_vertices = len(vertex_list)
    init_num_triangles = len(triangle_defs)
    
    print "before: %d verts, %d triangles" %(init_num_vertices, init_num_triangles)
    
    merge_short_triangles(vertex_list, triangle_defs, normals, too_short=minimum_length)
    rationalise(vertex_list, triangle_defs, normals)
    #fix_slender_triangles(vertex_list, triangle_defs)
    #rationalise(vertex_list, triangle_defs, normals)
    #remove_non_triplet_triangles(vertex_list, triangle_defs)
    #rationalise(vertex_list, triangle_defs, normals)
    
    final_num_vertices = len(vertex_list)    
    final_num_triangles = len(triangle_defs)
    print "after: %d verts, %d triangles" %(final_num_vertices, final_num_triangles)
    
    return

def find_close_vertices(vertex_list):
    
    mappings = {}
    processed = []
    while True:
        print processed
        quit_loop = True
        for i,v in enumerate(vertex_list):
            if i in processed: continue

            nearby = []
            idxs = [i]
            for j, other_v in enumerate(vertex_list[i+1:]):
                if (v - other_v).length() < 0.1:
                    nearby.append(other_v)
                    idxs.append(j+i)

            # if there are nearby vertices, then fix them
            if len(nearby):
                centre = Vector(0.0,0.0,0.0)
                for n in nearby:
                    centre += n
                centre /= len(nearby)
                new_vertex_idx = len(vertex_list)
                vertex_list.append(centre)
                for idx in idxs:
                    mappings[i] = new_vertex_idx
                processed.extend(idxs)
                quit_loop = False
                break
            
        # only stop merging vertices when no vertex is within 0.1 of another
        if quit_loop: break
        
    return

def merge_short_triangles(vertex_list, triangle_defs, normals, too_short):

    new_triangle_defs = []
    vertex_conversions = {}
    ctr = 0
    while True:
        start_ctr = len(triangle_defs)
        for t_vidx1, t_vidx2, t_vidx3 in triangle_defs:
    
            # renumber vertices if necessary...
            try:
                while(True):
                    t_vidx1 = vertex_conversions[t_vidx1]
            except KeyError:
                pass
            try:
                while(True):
                    t_vidx2 = vertex_conversions[t_vidx2]
            except KeyError:
                pass
            try:
                while(True):
                    t_vidx3 = vertex_conversions[t_vidx3]
            except KeyError:
                pass
    
            # if this triangle is already collapsed then don't bother
            if count_common_vertices([t_vidx1, t_vidx2, t_vidx3],
                                     [t_vidx1, t_vidx2, t_vidx3]) != 3:
                continue
            
            # these are the vertices
            v1 = vertex_list[t_vidx1]
            v2 = vertex_list[t_vidx2]
            v3 = vertex_list[t_vidx3]
    
            vertex_idx_pairs = [(t_vidx2,t_vidx1),
                                (t_vidx3,t_vidx2),
                                (t_vidx1,t_vidx3)]
            vertex_pairs = [(v2,v1),(v3,v2),(v1,v3)]
            lens = [(vv1 - vv2).length() for vv1, vv2 in vertex_pairs]
            min_len = min(lens)
    
            if min_len < too_short:
                ctr += 1
                idx = lens.index(min_len)
                merge_1, merge_2 = vertex_pairs[idx]
                
                new_vertex = (merge_1 + merge_2) / 2.0
                
                merged_vertex_idx = len(vertex_list)
                vertex_list.append(new_vertex)
                
                # generate new normal vector for the new vertex
                n1, n2 = vertex_idx_pairs[idx]
                merged_normal = ((normals[n1] + normals[n2]) / 2.0).normal()
                normals.append(merged_normal)
                
                for obsolete_vert_idx in vertex_idx_pairs[idx]:
                    vertex_conversions[obsolete_vert_idx] = merged_vertex_idx
    
            else:
                new_triangle_defs.append((t_vidx1,t_vidx2,t_vidx3))

        triangle_defs[:] = new_triangle_defs[:]                
        new_triangle_defs = []
        if len(triangle_defs) == start_ctr:
            break

    #print start_ctr
    
    for t_vidx1, t_vidx2, t_vidx3 in triangle_defs:

        try:
            while(True):
                t_vidx1 = vertex_conversions[t_vidx1]
        except KeyError:
            pass
        try:
            while(True):
                t_vidx2 = vertex_conversions[t_vidx2]
        except KeyError:
            pass
        try:
            while(True):
                t_vidx3 = vertex_conversions[t_vidx3]
        except KeyError:
            pass

        # this will abort any triangles which have collapsed to lines or
        # points
        if count_common_vertices([t_vidx1, t_vidx2, t_vidx3],
                                 [t_vidx1, t_vidx2, t_vidx3]) != 3:
            continue 
        
        new_triangle_defs.append((t_vidx1, t_vidx2, t_vidx3))
    
    triangle_defs[:] = new_triangle_defs[:]
        
    return

def average_side_length(vertex_list, triangle_defs):
    """Given vertex and triangle lists, calculate average side length.
    
    NB: This is inefficient as each side is shared between two triangles and
    so is calculated twice.  Oh well."""

    accum = 0.0
    for t_vidx1, t_vidx2, t_vidx3 in triangle_defs:
        
        v1 = vertex_list[t_vidx1]
        v2 = vertex_list[t_vidx2]
        v3 = vertex_list[t_vidx3]
        
        accum += (v2 - v1).length()
        accum += (v3 - v2).length()
        accum += (v1 - v3).length()
        
    accum /= (len(triangle_defs) * 3)
    
    return accum

def cull_vertex_list(vertex_list, triangle_defs, normals):

    # sanity check the vertex/normal lists
    assert(len(vertex_list) == len(normals))
    
    used = [False for n in range(len(vertex_list))]
    
    # record list of the vertices which are in use
    for v1,v2,v3 in triangle_defs:
        used[v1] = True
        used[v2] = True
        used[v3] = True

    # create culled vertex list (don't include unused vertices)
    vertex_list[:] = [vert for vert, in_use in zip(vertex_list,used) if in_use]
    normals[:] = [norm for norm, in_use in zip(normals,used) if in_use]

    return

def rationalise(vertex_list, triangle_defs, normals):

    # first cleanup useless triangles -- collapsed to lines/points or
    # duplicate entries
    cleanup_triangle_list(triangle_defs)
    
    mapping = {}
    for i,vertex in enumerate(vertex_list):
        mapping[i] = vertex

    cull_vertex_list(vertex_list, triangle_defs, normals)        
        
    # match vertices to new indices
    for k in mapping.keys():
        vert = mapping[k]
        try:
            i = vertex_list.index(vert)
        except ValueError:
            i = -1
        mapping[k] = i

    new_triangle_defs = []
    def renumber(vert):
        """Utility function for renumbering vertices."""
        new_num = mapping[vert]
        assert(new_num != -1)
        return new_num

    for v1,v2,v3 in triangle_defs:
        
        [new_v1, new_v2, new_v3] = [renumber(vv) for vv in [v1,v2,v3]]
        new_triangle_defs.append((new_v1, new_v2, new_v3))
        
    triangle_defs[:] = new_triangle_defs[:]
    return

# trivial helper function
def count_common_vertices(tri1, tri2):
    """Count common vertex indices -- tri, tr2 should be tuples of vertex
    indices such as: [2,51,231] and [1243,24,231] which would give a reult of
    1 (vertex 231 is shared).
    
    NB: Not necessarily equal length tuples, or any specific length ... counts
    repeats i.e. inputs [1,2,3] and [1, 1, 1] will return 3 not 1.
    """

    # this is what we're doing, in legible form...
    #ctr = 0
    #for v in tri1:
    #    if v in tri2:
    #        ctr += 1
    #return ctr
            
    # but we can use a list comprehension just to show how clever we are
    return sum([tri2.count(v) for v in tri1])

def fix_slender_triangles(vertex_list, triangle_defs):
    """Fix slender triangles -- if used after the merge_short_triangles
    function then this should deal with the slender-obtuse cases which
    remain.
    
    Find any triangle where the longest side is not very much longer than the
    sum of the other two sides; find the other triangle which shares that
    edge, then re-shuffle the vertices to produce 2 sensibly shaped triangles.
    It is possible that the resulting triangles will then be ill-formed in
    that one edge could be much shorter than the others. In which case,
    merge_short_triangles needs running again."""
    
    #triangle_vertex_sum = [(sum(list(verts)), i) for i, verts in enumerate(triangle_defs)]
    #triangle_vertex_sum.sort()
    vertex_conversions = {}
    processed = []

    ctr=0
    removals = []
    for i, tri in enumerate(triangle_defs):
        
        if (count_common_vertices(list(tri),list(tri)) != 3): 
            removals.append(tri)
            continue
        
        if ctr==1000:
            print i
            ctr=0
        ctr += 1
        
        if i in processed:
            continue

        # these are the vertex indices in vertex_list
        t_vidx1, t_vidx2, t_vidx3 = tri[:]

        # these are the vertices themselves (Vector-like objects)
        v1 = vertex_list[t_vidx1]
        v2 = vertex_list[t_vidx2]
        v3 = vertex_list[t_vidx3]

        vertex_idx_pairs = [(t_vidx2,t_vidx1),
                            (t_vidx3,t_vidx2),
                            (t_vidx1,t_vidx3)]
        vertex_pairs = [(v2,v1),(v3,v2),(v1,v3)]
        lens = [( (vv1 - vv2).length(), idx) 
                for idx, (vv1, vv2) in enumerate(vertex_pairs)]
        lens.sort(reverse=True)
        hypotenuse_verts = list(vertex_idx_pairs[lens[0][1]])
        remaining_vertex = [v for v in tri if v not in hypotenuse_verts][0]
        
        aspect_ratio = (lens[1][0]+lens[2][0]) / lens[0][0]
        if abs(lens[1][0] - lens[2][0]) < 0.5 and aspect_ratio < 1.1: # slender

            #print lens
            for j, other_tri in enumerate(triangle_defs):
                
                if i != j and count_common_vertices(hypotenuse_verts, list(other_tri)) == 2:
                    other_remaining_vertex = [v for v in other_tri 
                                              if v not in hypotenuse_verts][0]
                    
                    #print "before:", triangle_defs[i], triangle_defs[j]
                    triangle_defs[i] = (remaining_vertex, 
                                        hypotenuse_verts[0],
                                        other_remaining_vertex)
                    triangle_defs[j] = (remaining_vertex,
                                        hypotenuse_verts[1],
                                        other_remaining_vertex)
                    #print "after:", triangle_defs[i], triangle_defs[j]
                    processed.append(j)
                    break
                    #processed.sort()
                    
        
    print "Flipped %d triangles" %(len(processed))

    for r in removals:
        triangle_defs.remove(r)
    
    return

def cleanup_triangle_list(triangle_defs):
    
    # remove collapsed triangles
    new_triangle_defs = []
    for i, tri in enumerate(triangle_defs):
        if count_common_vertices(list(tri),list(tri)) == 3:
            new_triangle_defs.append(tri)
    triangle_defs[:] = new_triangle_defs[:]
    
    # remove repeated triangles
    new_triangle_defs = []
    for tri in triangle_defs:
        if tri not in new_triangle_defs:
            new_triangle_defs.append(tri)
    triangle_defs[:] = new_triangle_defs[:]

    return

def polyfilla(vertex_list, triangle_defs, min_tris_per_vertex=3):
    """Ensure that an umbrella of triangles surround all vertices."""

    from collections import deque
    for central_vertex_idx,v in enumerate(vertex_list):
        
        triangles = [(list(t),tri_idx) for tri_idx,t in enumerate(triangle_defs)
                     if central_vertex_idx in t]
        for t,idx in triangles:
            t.remove(central_vertex_idx) # remove central vertex from each triangle

        if len(triangles) == 0:
            # broken vertex, no triangles point to it
            # use the rationalise function to remove dead vertices
            continue
            
        expected_length = len(triangles)
           
        umbrella = deque()
        umbrella.append(triangles.pop()[0])

        forward = True
        back = True
        while forward or back:
            
            # grow forwards
            if forward:
                t = umbrella[-1]
                next_t = [(tri,idx) for tri,idx in triangles 
                          if count_common_vertices(t, tri)]
                if len(next_t) == 0:
                    forward = False
                elif len(umbrella)==1 and len(next_t)==2:
                    triangles.remove(next_t[0])
                    umbrella.append(next_t[0][0])
                elif len(next_t) > 1:
                    # bad geometry
                    forward = False
                else:
                    triangles.remove(next_t[0])
                    umbrella.append(next_t[0][0])

            # grow backwards
            if back:
                t = umbrella[0]
                next_t = [(tri,idx) for tri,idx in triangles 
                          if count_common_vertices(t, tri)]
                if len(next_t) == 0:
                    back = False
                elif len(next_t) > 1:
                    # bad geometry
                    back = False
                else:
                    triangles.remove(next_t[0])
                    umbrella.appendleft(next_t[0][0])

        # unused triangle definitions should be culled
        remove = []
        for tri,idx in triangles:
            remove.append(triangle_defs[idx])
        for r in remove:
            triangle_defs.remove(r)
                    
        # ok if the umbrella is linked at front and back then the vertex is
        # ok, otherwise we need a patch triangle to fill the hole.
        shared_verts = count_common_vertices(umbrella[0], umbrella[-1])        
        if len(umbrella) < min_tris_per_vertex or shared_verts > 1:
            # hmm that's very odd indeed
            # destroy the entire vertex, it's malformed
            remove = [t for t in triangle_defs if central_vertex_idx in t]
            for r in remove:
                triangle_defs.remove(r)
        elif shared_verts == 1:
            # no probs, all fine
            pass
            
        else:
            # no shared vertices, so we need a patch.
            
            # figure out the two orphaned vertex indices
            if umbrella[0][0] in umbrella[1]:
                lone_vert_a = umbrella[0][1]
            else:
                lone_vert_a = umbrella[0][0]
            assert(lone_vert_a not in umbrella[1])
            if umbrella[-1][0] in umbrella[-2]:
                lone_vert_b = umbrella[-1][1]
            else:
                lone_vert_b = umbrella[-1][0]
            assert(lone_vert_b not in umbrella[-2])
            
            # append the patch
            triangle_defs.append( (lone_vert_a, lone_vert_b, central_vertex_idx) )

    return       
                
            
def remove_non_triplet_triangles(vertex_list, triangle_defs):
    """Remove any triangles which do not have exactly three other triangles
    adjoining."""

    while True:
        removals = []
        for t in triangle_defs:
            
            other_triangles = [other for other in triangle_defs 
                               if count_common_vertices(list(t),list(other))==2]
    
            if len(other_triangles) != 3:
                removals.append(t)

        # keep looping until nothing to fix
        if len(removals) == 0: break                
                
        # remove the highly shonky looking triangles
        for r in removals:
            triangle_defs.remove(r)
        
        print "removed %d triangles" %(len(removals))
        
        print "polyfilling"
        polyfilla(vertex_list, triangle_defs)
        print "done polyfilling"

    return
