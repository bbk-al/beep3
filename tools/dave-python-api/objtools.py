#!/usr/bin/python

import string

def obj2kin(filename):

    # kintools contains the code to output kinemage file format
    import kintools
    
    faces, vertices, normals = getFacesVerticesNormals(filename)
    kin = kintools.triangleList(faces, vertices)
    return kin
    
def decodeObjFaceVertex(face_string):
    """Convert .obj file format face vertex string (face_string) 
       to a 0-based vertex index (e.g. 1//1 --> 0)"""

    idx = face_string.find("//")
    before_the_slashes = face_string[:idx]
    return string.atoi(before_the_slashes) - 1

def getFacesVerticesNormals(filename):
    
    objfile = open(filename, 'r')
    lines = objfile.readlines()

    faces = []
    vertices = []
    normals = []
    
    for line in lines:
        components = line.split()
        if len(components) > 1:

            if components[0] == 'vn':
                # vertex normal
                normal_vector = [string.atof(x) for x in components[1:]]
                normals.append(normal_vector)
        
            elif components[0] == 'v':
                # vertex
                vector = [string.atof(x) for x in components[1:]]
                vertices.append(vector)
        
            elif components[0] == 'f':
                # face 
                face = [decodeObjFaceVertex(s) for s in components[1:]]
                faces.append(face)
            
            else:
                pass
    
    return faces, vertices, normals

if __name__ == "__main__":
    
    import sys

    # if running as a script, call obj2kin
    try:
        filename = sys.argv[1] 
    
    except IndexError:
        print "No input filename specified"
        exit(1)

    print obj2kin(filename)
