#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pybeep import Vector, Octree, BasicTriangle

def average_normal(n1, n2, n3):
    """Returns the 'average' of three vectors, normalised to unit length."""
    return (n1 + n2 + n3).normalised()

def tnormal(a,b,c,direction=None):
    """Returns the normal vector given three vertices of a triangle.

    Ensures the normal is in the correct direction if supplied. Assumes a,b,c
    are Vector objects (with cross and normal methods defined correctly)."""

    # calculate normal vector from cross product
    #try:
    n = ((b - a).cross(c - a))
    n.normalise()
    #except:
    #    print a,b,c
    #    raise Exception

    # ensure correct direction
    if direction is not None and n.dot(direction) < 0:
        n *= -1

    return n

def tarea(a,b,c):
    cross_prod = (b - a).cross(c - a)
    return cross_prod.length() * 0.5

# utility functions for GTS file format
def isGTSComment(line):
    try:
        if (line[0] in ["#","!"]): return True
    except IndexError:
        pass
    return False

def get_next_line_not_comment(gts_file):
    while (True):
        line = gts_file.readline()
        if not isGTSComment(line): break
    return line

def get_gts_info(filename, return_handle=False):
    """Returns a tuple containing number of vertices, edges, faces in gts file.

    If return_handle is True then will return an open file handle too.
    """

    f = open(filename,"r")

    try:

        # get the first line, ignoring comment lines (starting with # or !)
        first_line = get_next_line_not_comment(f)

        # number of vertices, edges and faces
        nv, ne, nf = [int(xx) for xx in first_line.split()[:3]]

    except IOError:
        print("Bad GTS file.")
        return None
    finally:

        if not return_handle:
            f.close()

    if return_handle:
        return (nv, ne, nf), f
    else:
        return (nv, ne, nf)

def get_vertices_from_gts(filename):
    """returns a list of vertices and the number of faces from a gts file."""

    # first get the number of vertices, edges, faces
    (nv, ne, nf), f = get_gts_info(filename, return_handle=True)

    vertices = []
    try:

        while (len(vertices) < nv):

            # convert to Vertex object
            vx, vy, vz = [float(xx) for xx in get_next_line_not_comment(f).split()]
            vertices.append(Vector(vx, vy, vz))

    except IOError:
        print("Bad GTS file.")
        vertices = []

    finally:
        f.close()

    return vertices

def gts_overmesh(filename_in, filename_out):
    """Overmesh the mesh by subividing all triangles into 6 new triangles."""

    # first get the number of vertices, edges, faces
    (nv, ne, nf), f = get_gts_info(filename_in, return_handle=True)

    vertices, triangles = [], []

    mid_edge_vertices = {}
    tri_norms = []

    try:

        while (len(vertices) < nv):

            # convert to Vertex object
            vx, vy, vz = [float(xx) for xx in get_next_line_not_comment(f).split()]
            vertices.append(Vector(vx, vy, vz))

        edges = []
        # edges are defined as an ordered pair of vertices
        while (len(edges) < ne):
            edge_verts = get_next_line_not_comment(f).split()
            edges.append( set([int(edge_verts[0])-1, int(edge_verts[1])-1]) )

        # triangle defs come after the edge defs
        while (len(triangles) < nf*6):
            next_line = get_next_line_not_comment(f).split()
            e1_idx, e2_idx, e3_idx = [int(xx)-1 for xx in next_line]
            e1, e2, e3 = [edges[int(xx)-1] for xx in next_line]

            # get the 3 vertices, in order
            v1_idx = e1.intersection(e2).pop()
            v2_idx = e2.intersection(e3).pop()
            v3_idx = e3.intersection(e1).pop()

            #new_t = (v1_idx,v2_idx,v3_idx)
            t_centre = (vertices[v1_idx] + vertices[v2_idx] + vertices[v3_idx]) / 3.0
            t_centre_idx = len(vertices)
            vertices.append(t_centre)

            # add new vertices to vertex list if necessary
            def get_mid_edge_vertex(edge_idx):
                try:
                    e_mid = mid_edge_vertices[edge_idx]
                except KeyError:
                    e_mid = len(vertices)
                    v1_idx, v2_idx = edges[edge_idx]
                    vertices.append( (vertices[v1_idx] + vertices[v2_idx]) / 2.0)
                    mid_edge_vertices[edge_idx] = e_mid
                return e_mid

            def add_triangle(v1,v2,v3):
                triangles.append( (v1,v2,v3) )
                tri_norms.append( tnormal(vertices[v1],vertices[v2],vertices[v3]) )

            add_triangle(v1_idx, t_centre_idx, get_mid_edge_vertex(e1_idx))
            add_triangle(get_mid_edge_vertex(e1_idx), t_centre_idx, v3_idx)
            add_triangle(v3_idx, t_centre_idx, get_mid_edge_vertex(e3_idx))
            add_triangle(get_mid_edge_vertex(e3_idx), t_centre_idx, v2_idx)
            add_triangle(v2_idx, t_centre_idx, get_mid_edge_vertex(e2_idx))
            add_triangle(get_mid_edge_vertex(e2_idx), t_centre_idx, v1_idx)

    except IOError:
        print("Bad GTS file.")

    finally:
        f.close()

    write_gts(filename_out,vertices, triangles)

    return

def get_vertices_triangles_from_gts(filename):

    # first get the number of vertices, edges, faces
    (nv, ne, nf), f = get_gts_info(filename, return_handle=True)

    vertices, triangles = [], []

    try:

        while (len(vertices) < nv):

            # convert to Vertex object
            vx, vy, vz = [float(xx) for xx in get_next_line_not_comment(f).split()]
            vertices.append(Vector(vx, vy, vz))

        # edges are defined as an ordered pair of vertices
        edges = []
        while (len(edges) < ne):
            edge_verts = get_next_line_not_comment(f).split()
            edges.append( set([int(edge_verts[0])-1, int(edge_verts[1])-1]) )

        # triangle defs come after the edge defs
        while (len(triangles) < nf):
            e1, e2, e3 = [edges[int(xx)-1] for xx in get_next_line_not_comment(f).split()]

            # get the 3 vertices, in order
            v1_idx = e1.intersection(e2).pop()
            v2_idx = e2.intersection(e3).pop()
            v3_idx = e3.intersection(e1).pop()

            new_t = (v1_idx,v2_idx,v3_idx)
            triangles.append(new_t)

    except IOError:
        print("Bad GTS file.")
        vertices, triangles = [], []

    finally:
        f.close()

    return vertices, triangles

def jagged_sphere(num_subdivides, radius, filename):

    from random import gauss
    from math import pi, sqrt

    vertices, triangles = get_spherical_mesh(num_subdivides)
    randscale = sqrt(4*pi / len(vertices)) / 8.0

    for v in vertices:
        v.x += gauss(0,randscale)
        v.y += gauss(0,randscale)
        v.z += gauss(0,randscale)

    vnormals = vertices[:]
    vertices = [v*radius for v in vertices]

    write_gts(filename,
              vertices,
              triangles)
    print("Wrote jagged sphere with %d vertices" %(len(vertices)))
    return

def noisy_sphere(num_subdivides, radius, filename):

    from random import gauss
    from math import pi, sqrt

    vertices, triangles = get_spherical_mesh(num_subdivides)
    randscale = sqrt(4*pi / len(vertices)) / 8.0

    for v in vertices:
        v.x += gauss(0,randscale)
        v.y += gauss(0,randscale)
        v.z += gauss(0,randscale)
    vertices = [v/v.length() for v in vertices]

    for v in vertices:
        v.x = round(v.x * 1e4) / 1e4;
        v.y = round(v.y * 1e4) / 1e4;
        v.z = round(v.z * 1e4) / 1e4;

    vnormals = vertices[:]
    vertices = [v*radius for v in vertices]

    write_gts(filename,
              vertices,
              triangles)
    print("Wrote noisy sphere with %d vertices" %(len(vertices)))
    return

def reflect(gts_filename, output_filename):

    vertices, triangles = get_vertices_triangles_from_gts(gts_filename)
    vertices = [Vector(-v.x, v.y, v.z) for v in vertices]

    tnormals = [tnormal(vertices[t[0]],vertices[t[1]],vertices[t[2]])
                for t in triangles]

    write_gts(output_filename, vertices, triangles)
    return;

def hemisphere(num_subdivides, radius, filename):

    from random import random

    vertices, triangles = get_spherical_mesh(num_subdivides)

    # copy normals (same as vertices at this point)
    vnormals = vertices[:]

    # flatten negative z points into z=0 plane
    for v,vn in zip(vertices,vnormals):
        if v.z < 0.0:
            v.z = 0.0
            vn = (0,0,-1)

    # fix radius
    vertices = [v*radius for v in vertices]

    # write output
    write_gts(filename,
              vertices,
              triangles)
    print("Wrote hemisphere with %d vertices" %(len(vertices)))
    return

def get_spherical_mesh(num_subdivides, radius=1.0, spherical_file=""):
    """Get approximate spherical mesh from subdivided icosahedron (uses GTS).

    Returns vertices and triangles in simple list format."""
    #TODO: fix this
    raise Exception
    #import tempfile
    #tmp_file = tempfile.mktemp(suffix=".gts",
                               #prefix="gts_sphere_",
                               #dir=".")

    #gts.write_spherical_mesh(num_subdivides, tmp_file)

    #vertices, triangles = get_vertices_triangles_from_gts(tmp_file)
    #vertices = [v*radius for v in vertices]

    ## cleanup
    #import os
    #os.unlink(tmp_file)

    #if (spherical_file != ""):
        #write_gts(spherical_file, vertices, triangles, vertex_normals=vertices)

    #return vertices, triangles

def ellipsoid(filename, a, b, c):
    """Write an ellipsoid GTS file, based on a stretched sphere.

    NB:  Should be post-processed using meshlab to create uniform sampling"""

    # write noisy sphere
    noisy_sphere(4, 1.0, filename)

    verts, tris = get_vertices_triangles_from_gts(filename)
    for v in verts:
        v.x *= a
        v.y *= b
        v.z *= c

    write_gts(filename, verts, tris, verts)
    return


def clean_gts(filename):
    
    vertices, triangles = get_vertices_triangles_from_gts(filename)
    
    # get centre of set of points
    centre = Vector(0,0,0)
    for v in vertices:
        centre += v
    centre = centre / len(vertices)

    # get maximum xyz extent of set of points
    maxdim = max([ (v - centre).length() for v in vertices ])

    # insert vertices into an Octree
    max_neighbours = 27
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)

    duplicates = []
    ctr = 0
    mappingA = {}
    mappingB = {}
    for i,v in enumerate(vertices):

        # check for duplicate already in tree
        pre_existing = adaptive_tree.check_for_duplicate(v)
        if (pre_existing == -1):
            mappingA[i] = ctr
            mappingB[ctr] = i
            adaptive_tree.insert(v)
            ctr += 1
        else:
            pre_existing_idx = mappingB[pre_existing]
            print("Clash between %d %s and %d %s" %(i, v, pre_existing_idx, vertices[pre_existing_idx]))
            duplicates.append(i)
            mappingA[i] = pre_existing

    vertices = [v for i,v in enumerate(vertices) if i not in duplicates]

    print("Number of vertices: ", len(vertices))
    new_tri = []
    for i,t in enumerate(triangles):
        #print t[0]
        new_v1 = mappingA[t[0]]
        new_v2 = mappingA[t[1]]
        new_v3 = mappingA[t[2]]
        if (new_v1 in (new_v2, new_v3) or
            new_v2 in (new_v1, new_v3) or
            new_v3 in (new_v1, new_v2)
            ):
            print("skipped a dud triangle (%d) " %(i), t)
            continue
        new_tri.append([new_v1, new_v2, new_v3])

    # write the gts file
    write_gts(filename, vertices, new_tri)

def write_gts(filename, vertices, triangles):
    """Write the mesh as a GTS (Gnu Triangulated Surface) format file.

    vertices should be a list of Vector's describing vertex pts;

    triangles should be list of triples with vertex indices."""

    eout = []
    edges = {}

    def edge_find(va, vb):

        assert(va != vb)
        if (vb < va):
            vxx = va
            va = vb
            vb = vxx

        if (va not in edges.keys()): edges[va] = []
        stuff = [eidx for (eidx, vidx) in edges[va] if vidx==vb]
        if (len(stuff) > 0):
            assert(len(stuff) == 1)
            return stuff[0]
    
        # append edge to text output
        eout.append("%d %d" %(va, vb))

        edges[va].append( (len(eout),vb) )
        
        return len(eout)
    
    tout = []
    for t in triangles:
        tri = BasicTriangle(vertices[t[0]],vertices[t[1]],vertices[t[2]])
        #if (tri.area() == 0): continue
        v1_idx = t[0] + 1
        v2_idx = t[1] + 1
        v3_idx = t[2] + 1
        
        e1_idx = edge_find( v1_idx, v2_idx)
        e2_idx = edge_find( v2_idx, v3_idx)
        e3_idx = edge_find( v3_idx, v1_idx)
        tout.append("%d %d %d" %(e1_idx, e2_idx, e3_idx))
    
    gts = open(filename, 'w')
    print("%d %d %d GtsSurface GtsFace GtsEdge GtsVertex" %(len(vertices), len(eout), len(triangles)), file=gts)
    for v in vertices:
        print("%f %f %f" %(v.x, v.y, v.z), file=gts)
    print("\n".join(eout), file=gts)
    print("\n".join(tout), file=gts)
        
    gts.close()

    return
        

def get_clean_msms_vertices(msms_vertex_filename):

    vfloat = [line.split() for line in open(msms_vertex_filename,'r') if line[0] != '#']

    # delete the first lines of each list, which are column totals
    vfloat.remove(vfloat[0])

    vertices = []
    normals = []
    for vertex_descriptor in vfloat:
        x,y,z,xn,yn,zn = [float(v) for v in vertex_descriptor[:6]]
        vertices.append( Vector(x,y,z) )
        normals.append(  Vector(xn,yn,zn) )

    # get centre of set of points
    centre = Vector(0,0,0)
    for v in vertices:
        centre += v
    centre = centre / len(vertices)

    # get maximum xyz extent of set of points
    maxdim = max([ (v - centre).length() for v in vertices ])

    # insert vertices into an Octree
    max_neighbours = 27
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)


    remap = {}
    for i,(v,n) in enumerate(zip(vertices,normals)):

        # check for duplicate already in tree
        pre_existing = adaptive_tree.check_for_duplicate(v)
        if (pre_existing == -1):
            adaptive_tree.insert(v)
        else:
            pre_existing = new_vertex_map[pre_existing]
            print("Clash between %d %s %s and %d %s %s" %(i, v, n, pre_existing, vertices[pre_existing], normals[pre_existing]))
            remap[i] = pre_existing
            normals[pre_existing] += n

    vertices = [v for i,v in enumerate(vertices) if i not in remap.keys()]
    normals = [v.normalised() for i,v in enumerate(normals) if i not in remap.keys()]

    print("done cleaning duplicate vertices")

    return vertices, normals

def msms_remesh(msms_vertex_filename, gts_filename):
    """Triangulate an MSMS vertex list to GTS format."""

    vertices, normals = get_clean_msms_vertices(msms_vertex_filename)

    # get centre of set of points
    centre = Vector(0,0,0)
    for v in vertices:
        centre += v
    centre = centre / len(vertices)

    # get maximum xyz extent of set of points
    maxdim = max([ (v - centre).length() for v in vertices ])

    patch_vertex_list = [ [] for v in vertices ]

    # insert vertices into an Octree
    max_neighbours = 100
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)

    for i,(v,n) in enumerate(zip(vertices,normals)):
        adaptive_tree.insert(Vector(v.x, v.y, v.z))


    def get_neighbourlist(idx):
        near_pts = [0]*(max_neighbours)
        try:
            num = adaptive_tree.get(idx,near_pts)
        except IndexError:
            print("IndexError: %d (%d vertices)" %(idx, len(vertices)))
            raise IndexError
        #print near_pts
        near_pts = [pt for pt in near_pts[:num] if pt != idx]
        return near_pts

    for working_idx,(v,vn) in enumerate(zip(vertices, normals)):

        near_pts = get_neighbourlist(working_idx)

        # select set of working points which define this patch
        proximity_scores = []
        for pt_idx in near_pts:
            pt = vertices[pt_idx]
            r = (pt - v)
            score = (r.length()**2) / (1.0 - r.normalised().dot(vn))
            proximity_scores.append( (score, pt_idx) )
        proximity_scores.sort()
        closest = proximity_scores[0][1]

        points = [closest]
        near_pts.remove(closest)

        while True:

            if len(points) == 3:
                near_pts.append(points[0])


            # now find the next point which best forms a triangle
            for pt_idx in near_pts:
                pt = vertices[pt_idx]
                r = (pt - v)
                score = (r.length()**2) / (1.0 - r.normalised().dot(vn))
                proximity_scores.append( (score, pt_idx) )
            proximity_scores.sort()
            closest = proximity_scores[0][1]









    return

def upsample_fh_vals(low_res_mesh, high_res_mesh, low_res_results, high_res_results):
    """Forces all vertices in the mesh_to_realign to lie on the surface defined by reference_mesh."""

    low_res_results = [(float(x.split()[0]), float(x.split()[1])) for x in open(low_res_results,'r').readlines()]
    high_res_out = open(high_res_results, 'w')

    low_res_vertices, low_res_triangles = get_vertices_triangles_from_gts(low_res_mesh)
    high_res_vertices, high_res_triangles = get_vertices_triangles_from_gts(high_res_mesh)

    # get centre of set of points
    centre = Vector(0,0,0)
    for v in low_res_vertices:
        centre += v
    centre = centre / len(low_res_vertices)

    # get maximum xyz extent of set of points
    maxdim = max(max([ (v - centre).length() for v in high_res_vertices ]),
                 max([ (v - centre).length() for v in low_res_vertices ])
                 )

    # insert reference vertices into an Octree
    max_neighbours = 27
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)

    for i,v in enumerate(low_res_vertices):
        v.uid = i # Vector is python object so can add an attribute
        adaptive_tree.insert(v)

    new_vertices = []
    for i,v in enumerate(high_res_vertices):
        
        closest_vertex = adaptive_tree.nearest(v)
        closest_original_vertex = -1
        for ii,vv in enumerate(low_res_vertices):
            if ((closest_vertex - vv).length() < 1e-6): 
                closest_original_vertex = ii
                break
        if (closest_original_vertex == -1):
            print("Failed to find closest vertex.")
            raise Exception

        
        
        #print "dist to closest vertex: ", (low_res_vertices[closest_original_vertex] - v).length()
        possible_triangles = [t for t in low_res_triangles if closest_original_vertex in t]

        if len(possible_triangles) == 0:
            print(v, centre, maxdim*2)
            print(v, closest_original_vertex, low_res_vertices[closest_original_vertex], (low_res_vertices[closest_original_vertex] - v).length())
            raise Exception

        vert_checklist = []
        for t in possible_triangles:
            for v_idx in t:
                if v_idx == closest_original_vertex: continue
                if v_idx not in vert_checklist:
                    vert_checklist.append(v_idx)

        def get_closest_vertex_from_list(v, vlist):
            vert_dists = [ ( (v - low_res_vertices[vert_chk]).length(), vert_chk)
                           for vert_chk in vlist]
            vert_dists.sort()
            #print "sorted list: ",  vert_dists
            return vert_dists[0][1]

        second_closest_vertex = get_closest_vertex_from_list(v, vert_checklist)

        possible_triangles = [t for t in possible_triangles
                              if second_closest_vertex in t]
        vert_checklist = []
        for t in possible_triangles:
            for v_idx in t:
                if v_idx == closest_original_vertex: continue
                if v_idx == second_closest_vertex: continue
                if v_idx not in vert_checklist:
                    vert_checklist.append(v_idx)

        third_closest_vertex = get_closest_vertex_from_list(v, vert_checklist)

        closest_triangles = [t for t in possible_triangles
                             if third_closest_vertex in t]
        ct = closest_triangles[0]
        #print ct


        # at this point we have a vertex v which needs projecting onto a point
        # within the triangle defined by "ct"
        v1 = low_res_vertices[ct[0]]
        v2 = low_res_vertices[ct[1]]
        v3 = low_res_vertices[ct[2]]
        v1f = low_res_results[ct[0]][0]
        v1h = low_res_results[ct[0]][1]
        v2f = low_res_results[ct[1]][0]
        v2h = low_res_results[ct[1]][1]
        v3f = low_res_results[ct[2]][0]
        v3h = low_res_results[ct[2]][1]

        new_x_axis = (v2-v1).normalised()
        new_y_axis = (v3-v1).normalised()

        pt_from_v1 = v - v1
        new_f = v1f + pt_from_v1.dot(new_x_axis)*(v2f-v1f) + pt_from_v1.dot(new_y_axis)*(v3f-v1f);
        new_h = v1h + pt_from_v1.dot(new_x_axis)*(v2h-v1h) + pt_from_v1.dot(new_y_axis)*(v3h-v1h);

        print(new_f, new_h, file=high_res_out)

    high_res_out.close();


    return

def gts2off(gts_filename, off_filename):

    vertices, triangles = get_vertices_triangles_from_gts(gts_filename)

    vnormals = [Vector(0,0,0) for ii in range(len(vertices))]

    for tri in triangles:
        v1,v2,v3 = tri
        v1v = vertices[v1]
        v2v = vertices[v2]
        v3v = vertices[v3]
        vnormals[v1] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v2] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v3] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
    vnormals = [v.normalised() for v in vnormals]

    fout = open(off_filename,'w')
    print("nOFF", file=fout)
    print("3", file=fout)
    # last number is number of edges; unused by OFF format, but must be present
    print("%9d%9d%9d" %(len(vertices),len(triangles),0), file=fout)
    for v,vn in zip(vertices, vnormals):
        print("%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f" %(v.x, v.y, v.z, vn.x, vn.y, vn.z), file=fout)
    for t in triangles:
        print("    3 %8d%8d%8d" %(t[0],t[1],t[2]), file=fout)
    fout.close()

    return

def off2gts(off_filename, gts_filename):

    off_file = open(off_filename,'r')
    first_line = off_file.readline()
    num_verts, num_tris = [float(xx) for xx in first_line.split()]
    vertex_lines = [off_file.readline() for ii in range(num_verts)]
    tri_lines = [off_file.readline() for ii in range(num_tris)]

    vertices = [[float(xx) for xx in line.split()]
                for line in vertex_lines]
    triangles = [[int(xx) for xx in line.split()[1:4]]
                 for line in tri_lines]
    vnormals = [Vector(0,0,0) for ii in range(num_verts)]
    for tri in triangles:
        v1,v2,v3 = tri
        v1v = vertices[v1]
        v2v = vertices[v2]
        v3v = vertices[v3]
        vnormals[v1] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v2] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v3] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
    vnormals = [v.normalised() for v in vnormals]

    write_gts(gts_filename, vertices, triangles)

    return

def gts2msms(gts_filename, msms_filename):

    vertices, triangles = get_vertices_triangles_from_gts(gts_filename)

    vnormals = [Vector(0,0,0) for ii in range(len(vertices))]

    for tri in triangles:
        v1,v2,v3 = tri
        v1v = vertices[v1]
        v2v = vertices[v2]
        v3v = vertices[v3]
        vnormals[v1] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v2] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)
        vnormals[v3] +=  tnormal(v1v, v2v, v3v) * tarea(v1v,v2v,v3v)

    vnormals = [v.normalised() for v in vnormals]

    fout = open(msms_filename,'w')
    print("%d    %d" %(len(vertices),len(triangles)), file=fout)
    for ii,(v,vn) in enumerate(zip(vertices, vnormals)):
        print("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d  2" %(v.x, v.y, v.z, vn.x, vn.y, vn.z, 0, 0), file=fout)
    for t in triangles:
        print("%8d%8d%8d" %(t[0]+1,t[1]+1,t[2]+1), file=fout)
    fout.close()

    return

def align_mesh_to_reference(reference_mesh, mesh_to_realign):
    """Forces all vertices in the mesh_to_realign to lie on the surface defined by reference_mesh."""

    ref_vertices, ref_triangles = get_vertices_triangles_from_gts(reference_mesh)
    vertices, triangles = get_vertices_triangles_from_gts(mesh_to_realign)

    # get centre of set of points
    centre = Vector(0,0,0)
    for v in ref_vertices:
        centre += v
    centre = centre / len(ref_vertices)

    # get maximum xyz extent of set of points
    maxdim = max(max([ (v - centre).length() for v in vertices ]),
                 max([ (v - centre).length() for v in ref_vertices ])
                 )

    # insert reference vertices into an Octree
    max_neighbours = 27
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)

    for i,v in enumerate(ref_vertices):
        adaptive_tree.insert(Vector(v.x, v.y, v.z))

    new_vertices = []
    for i,v in enumerate(vertices):
        closest_original_vertex = adaptive_tree.get_uid_closest_to(v)
        #print "dist to closest vertex: ", (ref_vertices[closest_original_vertex] - v).length()
        possible_triangles = [t for t in ref_triangles if closest_original_vertex in t]

        if len(possible_triangles) == 0:
            print(v, centre, maxdim*2)
            print(v, closest_original_vertex, ref_vertices[closest_original_vertex], (ref_vertices[closest_original_vertex] - v).length())
            raise Exception

        vert_checklist = []
        for t in possible_triangles:
            for v_idx in t:
                if v_idx == closest_original_vertex: continue
                if v_idx not in vert_checklist:
                    vert_checklist.append(v_idx)

        def get_closest_vertex_from_list(v, vlist):
            vert_dists = [ ( (v - ref_vertices[vert_chk]).length(), vert_chk)
                           for vert_chk in vlist]
            vert_dists.sort()
            #print "sorted list: ",  vert_dists
            return vert_dists[0][1]

        second_closest_vertex = get_closest_vertex_from_list(v, vert_checklist)

        possible_triangles = [t for t in possible_triangles
                              if second_closest_vertex in t]
        vert_checklist = []
        for t in possible_triangles:
            for v_idx in t:
                if v_idx == closest_original_vertex: continue
                if v_idx == second_closest_vertex: continue
                if v_idx not in vert_checklist:
                    vert_checklist.append(v_idx)

        third_closest_vertex = get_closest_vertex_from_list(v, vert_checklist)

        closest_triangles = [t for t in possible_triangles
                             if third_closest_vertex in t]
        ct = closest_triangles[0]
        #print ct


        # at this point we have a vertex v which needs projecting onto a point
        # within the triangle defined by "ct"
        v1 = ref_vertices[ct[0]]
        v2 = ref_vertices[ct[1]]
        v3 = ref_vertices[ct[2]]
        norm = tnormal(v1,v2,v3)
        new_z_axis = norm
        new_x_axis = (v2-v1).normalised()
        new_y_axis = new_z_axis.cross(new_x_axis).normalised()
        pt_from_v1 = v - v1
        projected_pt = v1 + (new_x_axis * pt_from_v1.dot(new_x_axis)) + \
                     (new_y_axis * pt_from_v1.dot(new_y_axis))

        assert( (projected_pt - v1).dot(norm) < 1e-3 )
        assert( (projected_pt - v2).dot(norm) < 1e-3 )
        assert( (projected_pt - v3).dot(norm) < 1e-3 )

        new_vertices.append(projected_pt)
        #vertices[i] = projected_pt
#        if ((projected_pt - v).length() > 1e-3):
#            print "moved vertex by: ", (projected_pt - v).length()

    tnormals = [tnormal(new_vertices[t[0]],new_vertices[t[1]], new_vertices[t[2]])
                for t in triangles]

    write_gts(mesh_to_realign, new_vertices, triangles)

    return

def msms2xyzn(msms_vertex_filename, xyzn_filename):
    """Convert msms files to xyz-with-normal format"""

    vfloat = [line.split() for line in open(msms_vertex_filename,'r') if line[0] != '#']

    # delete the first lines of each list, which are column totals
    vfloat.remove(vfloat[0])

    vertices = []
    normals = []
    for vertex_descriptor in vfloat:
        x,y,z,xn,yn,zn = [float(v) for v in vertex_descriptor[:6]]
        vertices.append( Vector(x,y,z) )
        normals.append(  Vector(xn,yn,zn) )


    fout = open(xyzn_filename,'w')
    for v,n in zip(vertices, normals):
        print(v.x, v.y, v.z, n.x, n.y, n.z, file=fout)
    fout.close()

    return

def msms2gts(msms_vertex_filename, msms_face_filename, gts_filename):
    """Convert msms files to gts format (Gnu Triangulated Surface)"""

    vfloat = [line.split() for line in open(msms_vertex_filename,'r') if line[0] != '#']

    # delete the first lines of each list, which are column totals
    vfloat.remove(vfloat[0])

    vertices = []
    normals = []
    for vertex_descriptor in vfloat:
        x,y,z,xn,yn,zn = [float(v) for v in vertex_descriptor[:6]]
        vertices.append( Vector(x,y,z) )
        normals.append(  Vector(xn,yn,zn) )

    # get centre of set of points
    centre = Vector(0,0,0)
    for v in vertices:
        centre += v
    centre /= len(vertices)

    # get maximum xyz extent of set of points
    maxdim = max([ (v - centre).length() for v in vertices ])

    # insert vertices into an Octree
    max_neighbours = 27
    adaptive_tree = Octree(max_neighbours, Vector(centre.x, centre.y, centre.z), maxdim*2)

    duplicates = []
    ctr = 0
    mappingA = {}
    mappingB = {}
    for i,(v,n) in enumerate(zip(vertices,normals)):

        # check for duplicate already in tree
        pre_existing = adaptive_tree.check_for_duplicate(v)
        if (pre_existing == -1):
            mappingA[i] = ctr
            mappingB[ctr] = i
            adaptive_tree.insert(v)
            ctr += 1
        else:
            pre_existing_idx = mappingB[pre_existing]
            print("Clash between %d %s %s and %d %s %s" %(i, v, n, pre_existing_idx, vertices[pre_existing_idx], normals[pre_existing_idx]))
            duplicates.append(i)
            mappingA[i] = pre_existing

    vertices = [v for i,v in enumerate(vertices) if i not in duplicates]
    normals = [n for i,n in enumerate(normals) if i not in duplicates]

    tri = [line.split() for line in open(msms_face_filename,'r') if line[0] != '#']
    tri = [[int(float(t))-1 for t in tri_def] for tri_def in tri]
    tri.remove(tri[0])
    print("Number of vertices: ", len(vertices))
    print("Number of normals: ", len(normals))
    #new_tri = tri[:]
    new_tri = []
    for i,t in enumerate(tri):
        #print t[0]
        new_v1 = mappingA[t[0]]
        new_v2 = mappingA[t[1]]
        new_v3 = mappingA[t[2]]
        if (new_v1 in (new_v2, new_v3) or
            new_v2 in (new_v1, new_v3) or
            new_v3 in (new_v1, new_v2)
            ):
            print("skipped a dud triangle (%d) " %(i), t)
            continue
        new_tri.append([new_v1, new_v2, new_v3])

    # write the gts file
    write_gts(gts_filename, vertices, new_tri)

    return

def gts_wireframe(gts_file, kin_output):
    """Write a kinemage wireframe."""

    # get vertices and number of faces
    nv, ne, nf = get_gts_info(gts_file)

    print(nv, ne, nf)

    # get vertices
    vertices = get_vertices_from_gts(gts_file)

    gts = open(gts_file, 'r')
    get_next_line_not_comment(gts) # first line
    for skip_verts in range(nv):
        get_next_line_not_comment(gts)

    kin_out = open(kin_output, 'w')
    print("@kinemage", file=kin_out)
    print("@vectorlist", file=kin_out)

    for edge in range(ne):
        edge_pair = get_next_line_not_comment(gts).split()
        p1 = int(edge_pair[0]) - 1 # not zero indexed in gts file!
        p2 = int(edge_pair[1]) - 1
        v1 = vertices[p1]
        v2 = vertices[p2]
        x1,y1,z1 = [v1.x,v1.y,v1.z]
        x2, y2, z2 = [v2.x, v2.y, v2.z]
        print("{} P %f %f %f\n{} %f %f %f" %(x1, y1, z1, x2, y2, z2), file=kin_out)

    return

def gts_solid(gts_filename, kin_output, idx=0, num_indices=1):

    from random import choice
    import kintools

    # read entire gts file
    vertices, triangles = get_vertices_triangles_from_gts(gts_filename)

    # write kinemage
    kin_out = open(kin_output, 'w')
    print("@kinemage", file=kin_out)

    colour_names, col_str = kintools.colour_ranges(num_indices)
    print(col_str, file= kin_out)

    print("@trianglelist {triangles}", file=kin_out)

    for t in triangles:

        # create a line for this triangular face
        [v1, v2, v3] = t

        #colour = kintools.randomColour()
        colour = choice(colour_names[idx])

        a = (ax, ay, az) = vertices[v1]
        b = (bx, by, bz) = vertices[v2]
        c = (cx, cy, cz) = vertices[v3]

        print("X %f %f %f %f %f %f %s %f %f %f" % \
              (ax, ay, az, bx, by, bz, colour, cx, cy, cz), file=kin_out)

    kin_out.close()

def gts_lower_precision(gts_in, gts_out):

    input_lines = [xx.strip() for xx in open(gts_in,'r').readlines()]
    verts,edges,faces = [int(xx) for xx in input_lines[0].split()[:3]]

    fout = open(gts_out,'w')
    print(input_lines[0], file=fout)
    for line in input_lines[1:verts+1]:
        vx, vy, vz = [float(xx) for xx in line.split()]
        print("%9.3f %9.3f %9.3f" %(vx, vy, vz), file=fout)
    for line in input_lines[verts+1:]:
        print(line, file=fout)
    fout.close()
    return

#def gts_cleanup(infile, outfile):
    #gts.coarsen(infile, outfile, 0)
    #return

def gts_noclash(pqrfile, gtsfile, output_filename):

    # get charges
    import pqrtools
    atomlist = pqrtools.getCharges(pqrfile)

    # get mesh
    from mesh import GTS_Mesh
    mesh = GTS_Mesh(gtsfile)

    # enforce no charge clashes
    while True:
        if mesh.enforce_no_charge_clashes(atomlist,
                                          keep_gts_output=output_filename):
            break
        gts_cleanup(output_filename, output_filename)

    # all done
    return

#def gts_smooth(args):
    #print args
    #exe, filename_in, filename_out, lam, mu, niter = args
    #gts.smooth(filename_in, filename_out, float(lam),float(mu), int(niter))
    #return

def xyzr_to_vertices(xyzr_filename, gts_filename, probe_radius=1.4):

    master_vlist, master_tlist = get_spherical_mesh(2)

    spheres = open(xyzr_filename,'r').readlines()
    vlist = []
    tlist = []
    nlist = []
    for sph in spheres:
        offset = len(vlist)
        x,y,z,r = [float(xx) for xx in sph.split()]
        r += probe_radius
        for vv in master_vlist:
            n = Vector(vv)
            nlist.append(n)
            xx = (vv.x * r) + x
            yy = (vv.y * r) + y
            zz = (vv.z * r) + z
            vlist.append(Vector(xx,yy,zz))

        for v1_idx, v2_idx, v3_idx in master_tlist:
            tlist.append( (v1_idx+offset, v2_idx+offset, v3_idx+offset) )

    write_gts(gts_filename, vlist, tlist)

    return

def xyzn2xyzq(xyzn_filename, xyzq_filename, net_charge):

    pts = open(xyzn_filename,'r').readlines()
    q = net_charge / len(pts)
    fout = open(xyzq_filename,'w')
    for line in pts:
        x,y,z = [float(xx) for xx in line.split()[:3]]
        print(x,y,z,q,1.0, file=fout)
    fout.close()
    return

# use cases of gts_utils
if __name__ == "__main__":

    import sys
    mode = sys.argv[1]
    if (mode in locals()):
        locals()[mode](*sys.argv[2:])
    else:
        print("Can't run that: try one of these... ", locals())

