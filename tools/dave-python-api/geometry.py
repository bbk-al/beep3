# vim: set fileencoding=UTF8 :
#
# geometry.py
#
# Defines the various geometrical elements for Boundary Element Method
# electrostatics.  (Formerly called triangle.py)
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

#
# TODO:
# - Allow for non-zero ionic effects (inverse screening length (i.e. kappa) > 0)
# - Test cases + doctests
# - docutils?
# - un-hack the Vector class (probably implement my own better version :-) )
#

import math
from random import choice

# MMTK/Scientific provides my Vector class
# IMPORTANT NOTE: I've hacked the Vector class in the underlying library to
# allow mutable Vectors (i.e. to allow += syntax etc.)
from _VECTOR import Vector
#from vector import Vector

# the linalg routines do some useful least-squares type stuff for
# interpolating the electric field vectors at the vertices of the Mesh.
# Requires: GNU Scientific Library + python bindings (pygsl).
# http://www.gnu.org/software/gsl/
# http://pygsl.sourceforge.net/
from pygsl import linalg

# I might Weave some bits in C++ which would require the following Scipy
# modules.
from scipy import weave
#from scipy.weave import converters

def normalised(vector):
    """Return normalised (length=1.0) Vector version of the input.

    This function is here to remove some ambiguity in the term 'normal' where
    I sometimes mean perpendicular to a plane, and sometimes I mean,
    normalised to length=1.0 (i.e. a unit vector)."""

    return Vector(vector).normal()

def apply_quaternion_to_vector(quaternion, vector):
    """Return a vector rotated by a quaternion, as a Vector."""

    # I nicked this from wikipedia :-)
    # http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

    # weaved! (woven?)
    try:
        a,b,c,d = [float(xx) for xx in quaternion]
    except TypeError:
        return vector
    v1, v2, v3 = [float(xx) for xx in vector]
    vnew = [0.0, 0.0, 0.0]

    code = """
    double t2 =   a*b;
    double t3 =   a*c;
    double t4 =   a*d;
    double t5 =  -b*b;
    double t6 =   b*c;
    double t7 =   b*d;
    double t8 =  -c*c;
    double t9 =   c*d;
    double t10 = -d*d;
    vnew[0] = 2.0*( (t8 + t10)*v1 + (t6 -  t4)*v2 + (t3 + t7)*v3 ) + v1;
    vnew[1] = 2.0*( (t4 +  t6)*v1 + (t5 + t10)*v2 + (t9 - t2)*v3 ) + v2;
    vnew[2] = 2.0*( (t7 -  t3)*v1 + (t2 +  t9)*v2 + (t5 + t8)*v3 ) + v3;
    """

    weave.inline(code,
                 ['a', 'b', 'c', 'd',
                  'v1', 'v2', 'v3',
                  'vnew'])
    
    return Vector(vnew)

def local_to_universal(vector, 
                       ref_pt_local, 
                       ref_pt_universe,
                       rotation=None):
    """Convert vector in local co-ordinates to Universe co-ords.

    ref_pt_local and ref_pt_universe are the SAME point in
    each coordinate frame; furthermore the reference point is the centre of
    rotation of the local coordinate frame, so this is the single point in the
    local coordinate frame which is not changed by rotations."""

    # reference points are mandatory!
    assert(ref_pt_local is not None)
    assert(ref_pt_universe is not None)
    
    # get the vector w.r.t the stationary point (common point of reference
    # between local/universe coordinate frames, which is also centre of
    # rotation).
    new_vector = vector - ref_pt_local
    
    # apply rotation (if passed in) (rotation is the rotation of the local
    # coordinate frame w.r.t universe frame)
    if rotation is not None:
        rotated = apply_quaternion_to_vector(rotation, new_vector)
    else:
        rotated = new_vector # no rotation applied
        
    # we know the coordinates of the reference point in universe coords (it
    # gets passed in to this function as reference_point_universe); the
    # universe coordinates of the vector passed in is the universe-frame
    # offset from this point.
    result = rotated + ref_pt_universe
        
    return result


class NodePatch(object):

    """An umbrella-like area around a Mesh Vertex, in Universe co-ordinates.

    Notes:

    1.) NodePatch objects are constructed at each Vertex of a mesh and
    incorporate an umbrella-like set of triangles made up of thirds of the
    full-sized mesh triangles surrounding this vertex. All that is required to
    construct a NodePatch is the central_vertex for the patch, which should
    already encapsulate all the information on the triangular mesh at that
    point, in the local co-ordinate frame of the mesh.

    2.) Importantly, the NodePatch object is in the Universe co-ordinate frame,
    rather than the mesh local coordinate frame. This is because NodePatch
    objects are used in the BEM code, which requires the relative positions of
    all the diffusing entities (or rather, their mesh) in the Universe frame,
    not the local entity frame."""

    def __init__(self, central_vertex, ref_pt_local, ref_pt_universe, rotation):
        """Constructor for NodePatch.

        Required arguments:

        central_vertex -- a Vertex object around which this NodePatch is
        constructed.

        rotation -- a Quaternion describing the rotation to be applied to the
        source points to get them into Universe co-ordinate frame.

        ref_pt_local -- centre of rotation in local coordinates.
        ref_pt_universe -- same point in universe coordinates.

        Notes:

        Central_vertex is the central vertex for the NodePatch, whose
        co-ordinates are in the local mesh frame, which will be translated by
        this function into the Universe co-ordinate frame."""

        # store the actual Vertex object about which this NodePatch is based
        # so that when we have BEM results we can iterate over the NodePatch
        # object list and assign f and h values to the original Vertex objects
        # in the Mesh.
        self._real_vertex = central_vertex

        self._xyz = local_to_universal(central_vertex, 
                                       ref_pt_local,
                                       ref_pt_universe,
                                       rotation)
        self._normal = apply_quaternion_to_vector(rotation, central_vertex.normal)

        # assert that the rotation of the unit normal didn't do anything funny!
        assert(abs(self._normal.length() - 1.0) < 1e-12)

        # create a list of the triangle elements making up the umbrella-shaped
        # NodePatch region. At this point, the triangles will all be in the
        # local co-ordinate frame.
        local_coordinate_triangles = self._calculate_triangles(central_vertex)

        # now convert the local_coordinate_triangles into Universe co-ordinate
        # frame by rotating and translating appropriately.
        self._triangles = []
        for t in local_coordinate_triangles:

            # translate the local coordinates for each triangle into Universe
            # coordinates, and create new Triangles from the result.
            verts = [local_to_universal(v, 
                                        ref_pt_local, 
                                        ref_pt_universe,
                                        rotation)
                     for v in t.vertices]
            normal_vector = apply_quaternion_to_vector(rotation, t.normal)

            # assert that the normal vector is indeed unit length
            assert(abs(self._normal.length() - 1.0) < 1e-12)

            new_triangle = Triangle(verts[0], 
                                    verts[1], 
                                    verts[2], 
                                    normal_vector)
            
            self._triangles.append(new_triangle)

        # calculate the area of this NodePatch
        self._area = sum([t.area for t in self._triangles])

        # assert that the area of the NodePatch in Universe co-ordinates
        # matches what we expect from the local co-ordinate frame.
        # i.e. that I haven't totally f***ed it up.
        expected_area = sum([t.area for t in central_vertex.triangles]) / 3.0
        #print abs(self.area - expected_area), self.area, expected_area
        assert(abs(self.area - expected_area) < expected_area * 1e-12)
        #self._area = expected_area

    def set_BEM_results(self, f, h):
        """Set the BEM results (f and h) for the Vertex of this NodePatch.

        f is the potential, h is the derivative of the potential in the normal
        direction at the vertex. After these two values are set, the Mesh
        level method calculateElectricFields should be called which uses the
        results to interpolate a 'least squares' fit of the Electric field at
        each Vertex."""

        # set phi and it's derivative in the normal direction on the Vertex
        self._real_vertex.f = f 
        self._real_vertex.h = h 

    @property
    def cofg(self):
        """Return the centre of gravity of the NodePatch."""
        acc = Vector(0.0,0.0,0.0) 
        for t in self.triangles:
            #acc += t.centre
            acc += (t.centre - self.xyz) * (t.area / self.area)
        acc /= len(self.triangles)
        
        return self.xyz + acc
        #return acc
        
    def compare_BEM_results(self, f, h):
        """Compare the passed in f and h value to those already set.
        
        Used to compare the GMRES method to the direct matrix inversion."""
        
        fdev = (abs(self._real_vertex.f - f) / self._real_vertex.f) * 100
        hdev = (abs(self._real_vertex.h - h) / self._real_vertex.h) * 100
        print "f: %f h: %f" %(fdev, hdev)

    @property
    def area(self):
        """Return the area for this NodePatch."""
        return self._area

    @property
    def xyz(self):
        """Return the coordinates of the centre of this NodePatch as Vector.
        
        In universe coordinate frame."""
        return self._xyz

    @property
    def normal(self):
        """Return the normal vector at the vertex of this NodePatch as Vector."""
        return self._normal

    @property
    def solid_angle(self):
        """Return the solid angle at the vertex of this NodePatch."""
        return self._real_vertex.solid_angle

    @property
    def triangles(self):
        """Return a the triangles which make up this NodePatch as a list."""
        return self._triangles

    def kinemage(self,colour_scale=1.0, num_colours=100, f=True):
        """Return kinemage description of all triangles in patch."""
        
        #from kintools import randomColour
        #colour = randomColour()

        def round_to_int(num):
            hi = math.ceil(num)
            low = math.floor(num)
            if hi-num > num-low:
                return int(low)
            else:
                return int(hi)
            
        if f:
            val = self._real_vertex.f
        else:
            val = self._real_vertex.h
        c_idx = int(round(num_colours*val / colour_scale))
        if c_idx < 1:
            c_idx = abs(c_idx)
            c_name = "red"
        else:
            c_name = "blue"

        colour = "%s_%d" %(c_name, min([c_idx,num_colours]))
            
        return [t.kinemage(colour=colour) for t in self.triangles]

    def kinemage_triangle_norms(self):
        """Return kinemage description for all triangle normals."""
        
        return [t.kinemageNormal() for t in self.triangles]
    
    @staticmethod
    def _calculate_triangles(vertex):
        """Return a list of triangle patches constructed from a Vertex.

        Give a Vertex, iterate over the triangles which surround that vertex
        and construct an umbrella-like collection of smaller triangles to make
        up this NodePatch. The vertex passed in should be in local
        coordinates, and the resulting triangles should subsequently be
        converted to Universe co-ordinates."""

        tlist = []

        average_a = vertex.area / len(vertex.triangles)
        
        for t in vertex.triangles:
            
            #if t.area < (average_a / 5.0): continue

            verts = t.vertices
            verts.remove(vertex)

            a = vertex
            b = verts[0]
            c = verts[1]

            midab = (b + a) / 2.0
            midac = (c + a) / 2.0

            tlist.append(Triangle(vertex, midab, t.centre, t.normal))
            tlist.append(Triangle(vertex, t.centre, midac, t.normal))

        return tlist

    def subdivided_quadrature_points(self):
        """Return recursively subdivided quadrature points for this NodePatch.

        Note: this function isn't really used anymore as the quadrature stuff
        is done within the BEM C++ Boost.Python library instead. Also we don't
        currently bother with recursive subdivision of the integration area
        around the central vertex of a NodePatch because it probably doesn't
        seem to really help with the singularity problem."""

        triangles = self.triangles

        all_quads = []
        all_wts = []

        for t in triangles:

            original_area = t.area
            total_area_so_far = 0.0
            last_triangle = t
            max_subdivides = 3
            ctr = 0
            while (True):
                ctr += 1
                for nt in last_triangle.subdivide():

                    # if the central node isn't in this subdivided triangle,
                    # then it is subdivided enough, output the quads and
                    # weights
                    if self.node not in nt.vertices:

                        total_area_so_far += nt.area

                        wt_normalisation = nt.area / self.area
                        quads, wts = nt.getQuadraturePoints()
                        wts = [w * wt_normalisation for w in wts]
                        #new_wts = []
                        #for q,w in zip(quads, wts):
                            #l = (q - self.node).length()
                            #if l < 1:
                                #w /= l
                            #new_wts.append(w)
                        #wts = new_wts
                        #wts = [w * wt_normalisation *
                               #(1.0 / (q - self.node).length())
                               #for w,q in zip(wts,quads)]

                        all_quads.extend(quads)
                        all_wts.extend(wts)

                    else:

                        last_triangle = nt

                if ctr == max_subdivides:
                    break;

                #if total_area_so_far > original_area:
                #    # at this point we've subdivided enough
                #    break

        return all_quads, all_wts

    def quadrature_points(self):
        """Return two lists of quadrature points and weights for this NodePatch.

        Iterate over the triangles which make up this NodePatch, and calculate
        the quadrature points on each triangle, then return a list of the
        quadrature points, and a list of the corresponding weights (as two
        separate lists).

        Integrations carried out with these quadrature points must be
        multiplied by 2 times the area of the entire NodePatch to take into
        account the conversion from actual triangle shape to the (0,0), (1,0),
        (0,1) standard parametric triangle on which the quadrature points are
        based.

        Note: This function isn't really used any more as the equivalent is
        carried out within the BEM C++ Boost.Python library."""

        all_quads = []
        all_wts = []

        for t in self.triangles:

            wt_normalisation = t.area / self.area

            quads, wts = t.getQuadraturePoints()
            wts = [w * wt_normalisation for w in wts]

            all_quads.extend(quads)
            all_wts.extend(wts)

        return all_quads, all_wts

class Vertex(Vector):

    """Represents a (three dimensional) Vertex of a BEM Mesh.

    Can be used directly as a Vector object, with some extra bonus attributes.

    This class is closely related to the Triangle class, which is defined in
    terms of three of these Vertex objects. The electrostatic properties of
    the Vertex object should be set directly (i.e. just set Vertex.f and
    Vertex.h) then call the calculate_electric_field method.

    The triangles list attribute of the Vertex class has to be filled in
    manually when you create Triangles from a Vertex. You can create a
    NodePatch object from a Vertex object, as long as the triangles have been
    setup correctly."""

    from constants import epsilon0
    kcal_per_mole_conversion = 96485.3383 * 1e-3 * 0.2388458966275 \
                                          * 1e10 * 1.6e-19 / epsilon0
    
    def __init__(self, *args):
        """Constructor for Vertex, which should be used like a Vector."""

        # call parent class init function
        Vector.__init__(self, *args)

        # some special stuff for Vertex that isn't in Vector
        self.triangles = [] # triangles of which this is a Vertex
        self.f = 0.0        # we'll store the electric potential here.
        self.h = 0.0        # normal deriv of potential here.

    @property
    def f_kcal_per_mole(self):
        """Return phi at this Vertex in kCal/mol. 
        
        Electrostatic potential in normal direction."""
        return self.f * self.kcal_per_mole_conversion # kCal/Mol
    
    @property
    def h_kcal_per_mole_angstrom(self):
        """Return d(phi)/dn at this Vertex in kCal/(Mol.Angstrom).
        
        Derivative of electrostatic potential in normal direction."""
        
        # NB: internal lengthscales are already in Angstroms, so the
        # conversion to Angstroms is just a multiplication by 1.0.
        return self.h * self.kcal_per_mole_conversion # kCal Mol^-1 Angstrom^-1
        
    def translate(self, xyz):
        """Move this Vertex by xyz.

        Uses the underlying __setstate__ method. Which is a bit naughty as
        it's not really a public method of the Vector class."""
        self.__setstate__(self + xyz)

    @property
    def area(self):
        """Return area associated with this vertex."""
        
        area = 0.0
        for t in self.triangles:
            area += t.area / 3.0
        return area

    @property
    def solid_angle_new(self):
        
        ave_centre = Vector(0.0,0.0,0.0)
        
        ang = 0.0
        for t in self.triangles:

            ave_centre += t.centre
            
            vset = t.vertices[:]
            vset.remove(self)
            
            vec1 = (vset[0] - self).normal()
            vec2 = (vset[1] - self).normal()
            ang += math.acos(vec1 * vec2)

        ave_centre /= len(self.triangles)
            
        #return ang

        if (ave_centre - self).length() < 1e-6:
            return ang
        
        if (ave_centre - self) * self.normal < 0:
            return (4.0*math.pi) - ang
        else:
            return ang
            
    @property
    def solid_angle(self):
        """Solid angle for this Vertex."""
        
        Ap = 0.0
        np = len(self.triangles)

        # helper function
        def find_tri_with_these_verts(tri, v1, v2):
            for t in self.triangles:
                if (t is not tri and
                    v1 in t.vertices and
                    v2 in t.vertices):
                    return t
            return None

        # force contiguity of the triangles we walk over
        first_triangle = self.triangles[0]
        done_list = [self]
        this_tri = first_triangle
        
        ctr = 0
        while True:
            ctr += 1
            vertlist = this_tri.vertices[:]
            vertlist = [v for v in vertlist if v not in done_list]
            
            next_t = find_tri_with_these_verts(this_tri, self, vertlist[0])
            done_list.append(vertlist[0])

            t1 = this_tri
            t2 = next_t
            
            # double check that the triangles have 2 common vertices
            count_common_vertices=0
            for v in t1.vertices:
                if v in t2.vertices:
                    count_common_vertices += 1
            assert(count_common_vertices == 2)
            
            norms = t1.normal * t2.normal
            if abs(norms) > 1.0: 
                norms /= abs(norms)
            ang = math.acos(-norms)
            Ap += ang
            
            if next_t is first_triangle:
                # we've walked around the nodepatch to where we started,
                # time to quit
                assert (ctr == np)
                break
            else:
                this_tri = next_t
        
        Ap -=  (np - 2)*math.pi
        while(Ap < 0.0):
            Ap += math.pi*4.0
        return Ap

    @property
    def normal(self):
        """Calculates the normal vector at this vertex.

        Computed by averaging over all triangles of which this vertex is a
        part, scaled by Area of each component triangles."""

        # holder for cumulative sum
        acc = Vector([0.0, 0.0, 0.0])

        # iterate over all triangles for which this Vertex object is
        # a vertex.  NB: self.triangles must have been populated already,
        # it's not done automatically in any way.
        for t in self.triangles:
            acc += t.normal * t.area # normal as in perpendicular

        # return the normalised vector acc
        return normalised(acc)

    def normalised(self):
        """Returns the normalised vector (length=1.0) for this Vertex."""

        # annoyingly the Vector.normal() method (Scientific.Geometry.Vector)
        # means normalised Vector, rather than normal in the perpendicular
        # sense. Which leads to some confusion in this code methinks... When I
        # say normal() I mean normal (perpendicular) and when I mean "length
        # 1.0" I say normalised(). So there.
        return super(Vertex, self).normal()

    def kinemageNormal(self, 
                       centre_of_rotation=None, 
                       rotation=None, 
                       translation=None):
        """Produce a kinemage polyline to represent the normal vector for this
        vertex"""

        here = local_to_universal(self,
                                  ref_pt_local=centre_of_rotation, 
                                  ref_pt_universe=translation,
                                  rotation=rotation)
        new_normal = apply_quaternion_to_vector(rotation, self.normal)

        (x1, y1, z1) = here
        (x2, y2, z2) = (here) + new_normal

        return "{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2)

    def calculate_electric_field(self):
        """Calculate the 'least squares' E vector at this Vertex.

        Assumes that the f and h attributes of the Vertex have already been
        set (probably by a call to NodePatch.set_BEM_results).

        Note that this function does some GSL calls (linalg module)."""

        # create a matrix of normal vectors for the triangles
        # around this vertex
        normals = linalg.zeros((len(self.triangles),3), linalg.Float)
        rhs = []

        for i, t in enumerate(self.triangles):

            # get normal vector for this triangle
            n = t.normal

            normals[i,0] = n.x()
            normals[i,1] = n.y()
            normals[i,2] = n.z()

            rhs.append(self.h)

        assert(len(self.triangles) > 2)
            
        if len(self.triangles) > 2:
            (U,V,S) = linalg.SV_decomp(normals)
            H = Vector(linalg.SV_solve(U, V, S, linalg.array(rhs, linalg.Float)))

        # normalised gives me the unit vector of H
        Hdir = normalised(H)
        
        E = Vector(0.0, 0.0, 0.0)
        for t in self.triangles:

            # get vertices of the triangle, removing the one we're at
            verts = t.vertices
            verts.remove(self)
            assert(len(verts) == 2) # otherwise this is a v strange triangle!

            # create a linalg (GSL) matrix
            normals = linalg.zeros((3,3), linalg.Float)

            normals[0,0] = Hdir.x()
            normals[0,1] = Hdir.y()
            normals[0,2] = Hdir.z()

            # the .normal() methods here normalise the lengths to 1.0. i.e.
            # they're calls to the Scientific.Geometry.Vector.normal() method
            normals[1,0] = (verts[0] - self).normal().x()
            normals[1,1] = (verts[0] - self).normal().y()
            normals[1,2] = (verts[0] - self).normal().z()

            normals[2,0] = (verts[1] - self).normal().x()
            normals[2,1] = (verts[1] - self).normal().y()
            normals[2,2] = (verts[1] - self).normal().z()

            rhs = []
            rhs.append(H.length())
            rhs.append( (verts[0].f - self.f) / (verts[0] - self).length() )
            rhs.append( (verts[1].f - self.f) / (verts[1] - self).length() )

            (U,V,S) = linalg.SV_decomp(normals)
            E += Vector(linalg.SV_solve(U, V, S, linalg.array(rhs, linalg.Float)))
            print U, V ,S

        self.E = E / len(self.triangles)
        #self.h = H * self.normal # bit of a sneaky hack

        return
    
    def get_tri_sharing_vertex(self, other_vertex, existing_tri):
        
        common_triangles = []
        for t in self.triangles:
            if self in t.vertices and other_vertex in t.vertices:
                common_triangles.append(t)

        common_triangles.remove(existing_tri)

        # check that there is exactly one
        assert(len(common_triangles)==1)
        
        return common_triangles[0]
    
                
class BadVertexException(Exception):
    
    def __init__(self, *args):
        Exception.__init__(self, *args)
        return
        
class Triangle(object):

    """A (planar) Triangle object intended to be used as a component of a Mesh
    for the Boundary Element Method."""

    quadrature_points = {}
    # (N,D)
    quadrature_points[(6,4)] = [ ((0.8168476,  0.09157621), 0.05497587),
                                 ((0.09157621, 0.8168476 ), 0.05497587),
                                 ((0.09157621, 0.09157621), 0.05497587),
                                 ((0.1081030,  0.4459485 ), 0.1116908),
                                 ((0.4459485,  0.1081030 ), 0.1116908),
                                 ((0.4459485,  0.4459485 ), 0.1116908) ]

    def __init__(self, a, b, c, normal=None):
        """Constructor for Triangle, requiring the vertices and a normal Vector.

        The Vertices can be in any order, they will be stored internally such
        that the vertices property of the Triangle object returns the vertices
        in a right-handed sense with respect to the normal vector supplied by
        the user.

        The normal Vector doesn't *have* to really be truly a normal, it's
        just used to determine the directionality of the normal vector, which
        is recalculated from the vertices in any case. i.e. the normal vector
        passed in controls the internal ordering of the vertices such that the
        normal vector points in the correct direction."""

        # vertices
        vertex_permutations = [[a,b,c],[b,c,a],[c,a,b]]

        # randomize the vertex ordering to remove any systematic bias in
        # calculations of e.g. the centre of the triangle
        #self.a, self.b, self.c = choice(vertex_permutations)

        # uncomment this bit if we don't want random vertex ordering...
        self.a = a
        self.b = b
        self.c = c
        
        if (self._calc_area() == 0.0):
            raise BadVertexException

        # figure out the normal vector for the triangle.
        # check against the normal at a vertex to get directionality
        # and flip vertex ordering if necessary
        self._normal_vector = self._calc_normal()        

        # if normal vector indicator is passed in
        if normal is not None:
            if normal.dot(self.normal) < 0:
                self.a, self.b, self.c = self.c, self.b, self.a # flip vertex order
            self._normal_vector = self._calc_normal()
            
        self._centre = self._calc_centre()
        self._area = self._calc_area()
        
        return

    def _force_recalc(self):
        """Force recalculation of triangle area, normal, centre to reflect
        changes to vertices."""

        # centre first
        self._centre = self._calc_centre()
        
        old = self._normal_vector
        self._normal_vector = self._calc_normal()
        
        # check we don't need to flip the vertex ordering
        if (self._normal_vector.dot(old)) < 0:
            self.a, self.b, self.c = self.c, self.b, self.a
            self._normal_vector = self._calc_normal()
            
        
        self._area = self._calc_area()
        return

    def subdivide_in_situ(self):
        """Subdivide this triangle into 3 triangles.  
        
        Returns the two new triangles in a list."""
        
        new_vert = Vertex(self.centre)
        
        old_norm = self.normal
        # de-register a vertex and make it the new central one
        # be respectful of the internal triangle list maintained by the vertex
        a,b,c = self.vertices
        c.triangles.remove(self)
        self.c = new_vert
        self._force_recalc()
        assert((self.normal-old_norm).length() < 1e-6)
        new_vert.triangles.append(self)
        
        # make two new triangles; register their vertices
        new_one = Triangle(a, new_vert, c, self.normal)
        for v in new_one.vertices:
            v.triangles.append(new_one)
        new_two = Triangle(c, new_vert, b, self.normal)
        for v in new_two.vertices:
            v.triangles.append(new_two)
        
        return new_vert, [new_one, new_two]
    
    def subdivide(self):
        """Return a list of 4 Triangle objects derived from this Triangle.

        The Triangle objects are defined by creating new vertices at the
        midpoint of each side of the existing Triangle and creating 4 new
        Triangles from the new vertex set."""

        new_triangles = []
        mid_ab = Vector((self.a + self.b) / 2.0)
        mid_ac = Vector((self.a + self.c) / 2.0)
        mid_bc = Vector((self.b + self.c) / 2.0)

        new_triangles.append(Triangle(mid_ac, mid_bc, self.c, self.normal))
        new_triangles.append(Triangle(self.a, mid_ab, mid_ac, self.normal))
        new_triangles.append(Triangle(mid_ab, self.b, mid_bc, self.normal))
        new_triangles.append(Triangle(mid_bc, mid_ac, mid_ab, self.normal))
        
        return new_triangles
    
    def _calc_normal(self):
        """Return the normal (perpendicular Vector to the plane.
        
        Calculated here and now."""
        return normalised( (self.b - self.a).cross(self.c - self.a) )
    
    @property
    def normal(self):
        """Return the normal (perpendicular) Vector to the plane.
        
        Pre-calculated and cached by __init__ function."""
        return self._normal_vector

    @property
    def vertices(self):
        """Return a list of the vertices for this Triangle.
        
        Could be Vector objects or full Vertex objects depending on usage."""
        return [self.a, self.b, self.c]

    def replace(self, old, new):
        """Replace vertex old with new. Retain ordering of vertex for consistent normal."""

        if old is self.a:
            self.a = new
        elif old is self.b:
            self.b = new
        elif old is self.c:
            self.c = new
        else:
            raise ValueError

        # update triangle lists within the vertex
        old.triangles.remove(self)
        new.triangles.append(self)
        
        # check normal is consistent
        old_normal = self.normal
        self._force_recalc()
        if self.normal * old_normal < 1.0:
            self.a,self.b,self.c = self.c, self.b, self.a
            self._normal_vector = self._calc_normal()

        return
    
    @property
    def allVertexPermutations(self):
        """Generator for all permutations of the vertex list.

        I have this to remove possible bias in the Vertex ordering."""

        def all_perms(str):
            if len(str) <=1:
                yield str
            else:
                for perm in all_perms(str[1:]):
                    for i in range(len(perm)+1):
                        yield perm[:i] + str[0:1] + perm[i:]

        return all_perms(self.vertices)

    @property
    def area(self):
        """Return the area of the Triangle.

        Precalculated in __init__ and cached."""
        return self._area

    def _calc_area(self):
        """Calculate the area of this Triangle."""

        vertex_set = self.vertices
        ac = Vector(vertex_set[0] - vertex_set[1])
        ab = Vector(vertex_set[0] - vertex_set[2])
        product = ac.cross(ab)
        area = 0.5 * product.length()
        return area

    def _calc_centre(self):
        """Calculate centre for later usage."""

        a, b, c = self.vertices

        ac_midpoint = (a + c) / 2.0
        b_to_ac_midpoint = ac_midpoint - b

        return b + (b_to_ac_midpoint * 2.0 / 3.0)
    
    @property
    def centre(self):
        """Return a Vector to the centre (centroid) of this Triangle.
        
        Precalculated in __init__ and cached."""
        return self._centre
        
    def hzToActual(self, h, z):
        """Convert parametric h,z coords to real 3D co-ordinates.

        Note: Probably redundant now; this function was used to generate the
        quadrature points but now that's done in the BEM C++ Boost.Python
        library."""

        a, b, c = self.vertices
        return (1.0 - h - z)*a + h*b + z*c

    @staticmethod
    def getRawQuadraturePoints(N,D):
        """Return list of Gaussian Quadrature points for unit parametric triangle.

        Note: Probably redundant now; this function was used to generate the
        quadrature points but now that's done in the BEM C++ Boost.Python
        library."""

        quads = Triangle.quadrature_points[ (N,D) ]
        return quads

    def getQuadraturePoints(self):
        """Return list of Gaussian Quadrature points in actual 3D space.

        Note: Probably redundant now; this function was used to generate the
        quadrature points but now that's done in the BEM C++ Boost.Python
        library."""

        quads = Triangle.getRawQuadraturePoints(6,4)
        pts = [self.hzToActual(h,z) for (h,z),wt in quads]
        wts = [wt for (h,z), wt in quads]

        return pts, wts

    @property
    def E(self):
        """Return Electric field Vector at the centre of this triangle.

        This is calculated from the sum of the E field at each Vertex, divided
        by the number of vertices (which had better be 3 for a triangle...).
        Note that the attribute E (Vector representing Electric Field) at each
        Vertex must have already been calculated and set, probably by a call
        to the higher-level Mesh function calculate_electric_fields()."""

        EField = Vector(0.0, 0.0, 0.0)
        for v in self.vertices:
            EField += v.E
        EField /= len(self.vertices)

        return EField

    @property
    def f(self):
        """Return the value of the potential (phi, or f) at centre of triangle.

        This is calculate from average of the values at the vertices. The
        Vertex attribute f *must* have already been calculated and set by the
        BEM methods."""
        return sum([v.f for v in self.vertices]) / len(self.vertices)

    def force(self, kappa, Dext):
        """Return a Vector of the force (in units of kT per Angstrom) on this
        triangle calculated via the Maxwell Stress Tensor.

        This method makes use the E property (Electric Field) for this
        Triangle, so note the requirement there that the field at each Vertex
        must have already been set before you can use this function.
        
        Note on units: up to this point, charges are multiples of the
        electronic charge; lengths are all in Angstroms; electronic
        permittivity is in relative not absolute units (i.e. Dext is ~80 for
        water)."""

        from universe import Universe
        
        E = self.E
        E2 = E * E
        Ex = E.x()
        Ey = E.y()
        Ez = E.z()

        # this is the MST premultiplier, in internal units (Angstroms,
        # relative permittivities etc.)
        mult = self.normal * self.area * Dext 

        # x stress
        stress_xx = Ex * Ex - 0.5*E2 - 0.5*kappa*kappa*self.f*self.f
        stress_xy = Ex * Ey
        stress_xz = Ex * Ez

        fx = mult * Vector(stress_xx, stress_xy, stress_xz)

        # y stress
        stress_yx = Ey * Ex
        stress_yy = Ey * Ey - 0.5*E2 - 0.5*kappa*kappa*self.f*self.f
        stress_yz = Ey * Ez

        fy = mult * Vector(stress_yx, stress_yy, stress_yz)

        # z stress
        stress_zx = Ez * Ex
        stress_zy = Ez * Ey
        stress_zz = Ez * Ez - 0.5*E2 - 0.5*kappa*kappa*self.f*self.f

        fz = mult * Vector(stress_zx, stress_zy, stress_zz)

        return Vector(fx, fy, fz)

    def kinemage(self, 
                 centre_of_rotation=Vector(0,0,0), 
                 rotation=None, 
                 translation=Vector(0,0,0), 
                 colour=None, 
                 scale=1.0):
        """Produce three lines in kinemage to represent this triangle.
        
        Note that the internal coordinates of the mesh may be converted to the
        universe coordinates by supplying a centre of rotation, rotation (as
        quaternion) and translation vector.
        
        Scale is the scale factor to apply to the raw internal units to get
        sensible sized values for kinemage."""

        # if no colour was passed in then get one randomly from kintools
        if colour is None:

            from kintools import randomColour
            colour = randomColour()
            
        # use lambda function for brevity in syntax
        cvt = lambda x: local_to_universal(x, 
                                           ref_pt_local=centre_of_rotation, 
                                           ref_pt_universe=translation,
                                           rotation=rotation)
        a = (ax, ay, az) = cvt(self.a) * scale
        b = (bx, by, bz) = cvt(self.b) * scale
        c = (cx, cy, cz) = cvt(self.c) * scale

        return "X %f %f %f %f %f %f %s %f %f %f" %(ax, ay, az,
                                                   bx, by, bz, colour,
                                                   cx, cy, cz )

    def kinemage_right_angled_triangles(self):
        """Kinemage triangles of two Right-Angled Triangles made from this.
        
        (Two right angled triangles sometimes give better results for the
        numerical quadrature algorithm.)"""

        t1,t2 = self.getTwoRightAngledTriangles()

        from kintools import randomPairPastels
        colour1, colour2 = randomPairPastels()

        return t1.kinemage_int(colour1) + "\n" + t2.kinemage_int(colour2)

    def kinemageNormal(self, 
                       centre_of_rotation=None, 
                       rotation=None, 
                       translation=None, 
                       scale=1.0):
        """Return kinemage polyline representing the normal vector."""

        # use lambda function for brevity in syntax
        cvt = lambda x: local_to_universal(x, 
                                           ref_pt_local=centre_of_rotation, 
                                           ref_pt_universe=translation,
                                           rotation=rotation)
        
        # c is the centre of the triangle (in universe coordinates)
        # rotated_normal is the normal vector, rotated to universe coords
        c = cvt(self.centre) * scale
        rotated_normal = apply_quaternion_to_vector(rotation, self.normal)
        (x1, y1, z1) = c
        (x2, y2, z2) = c + rotated_normal

        return "{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2)

    def kinemageFieldDirection(self, 
                       centre_of_rotation=None, 
                       rotation=None, 
                       translation=None, 
                       scale=1.0):
        """Return kinemage polyline representing the E Vector direction."""
        
        # use lambda function for brevity in syntax
        cvt = lambda x: local_to_universal(x, 
                                           ref_pt_local=centre_of_rotation, 
                                           ref_pt_universe=translation,
                                           rotation=rotation)
        
        # c is the centre of the triangle (in universe coordinates)
        # rotated_normal is the normal vector, rotated to universe coords
        c = cvt(self.centre) * scale
        rotated_field = apply_quaternion_to_vector(rotation, self.E.normal())
        (x1, y1, z1) = c
        (x2, y2, z2) = c + rotated_field*0.2
        
        return "{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2)

    def getTwoRightAngledTriangles(self):
        """Return two right angled triangles made from this triangle.

        Figures out which angle is the largest, then divides the triangle in
        two by cutting through that angle such that the new line joins the
        opposite edge at a right angle. (The angle is not necessarily
        bisected)."""

        a,b,c = self.vertices

        # sort the vertices by opposite edge length
        sideA = (c-b).length()
        sideB = (a-c).length()
        sideC = (b-a).length()
        lengths = [ (sideA, a), (sideB, b), (sideC, c) ]
        lengths.sort(reverse=True)

        # make vertex A opposite the largest side
        a = lengths[0][1]
        b = lengths[1][1]
        c = lengths[2][1]

        # now find out where the perpendicular bisector of side bc (the longest) lies
        direction_bc = (c-b).normal()

        # |ba| * |bc| * cos(angleABC) = length from b to new vertex, if |bc| is unity
        len_b_new_vertex = (a-b) * direction_bc
        new_vertex = b + (len_b_new_vertex * direction_bc)

        # these will be a pair of right angled triangles
        t1 = Triangle(a, b, new_vertex, self.normal)
        t2 = Triangle(c, a, new_vertex, self.normal)

        return (t1, t2)

# These are some useful functions derived from Numerical Recipes §21.6
import random
def point_inside_unit_circle():
    """Generate a random point within unit circle.
    
    Pick two randomly uniformly distributed numbers in range [-1,1]; reject
    any where x**2 + y**2 > 1.0.
    
    NB: The random.uniform function returns values up to, but not including,
    the upper boundary, so this is strictly slightly imperfect. Good enough
    for my purposes though."""
    
    while True:
        x = random.uniform(-1.0,1.0)
        y = random.uniform(-1.0,1.0)
        if (x*x + y*y) <= 1.0: break
    
    return (x,y)

def random_point_on_sphere(dimensions=3):
    """Return random point on sphere of unit radius in given number of dimensions."""

    if dimensions==3:

        u0,u1 = point_inside_unit_circle()
        useful = math.sqrt(1.0 - u0*u0 - u1*u1)
        
        x = 2.0*u0*useful
        y = 2.0*u1*useful
        z = 1.0 - 2.0*(u0*u0 + u1*u1)
        
        return (x,y,z)
    
    elif dimensions==4:
        
        u0,u1 = point_inside_unit_circle()
        u2,u3 = point_inside_unit_circle()
        
        x0 = u0
        x1 = u1
        useful = math.sqrt((1.0 - u0*u0 - u1*u1)/(u2*u2 + u3*u3))
        x2 = u2*useful
        x3 = u3*useful        
        
        return (x0,x1,x2,x3)
    
    else:
        # I haven't implemented routines for other dimensions.
        raise ValueError

def random_rotation_matrix():
    """Return a random rotation matrix in 3 dimensions."""
    
    from Scientific.Geometry import Tensor
    x0,x1,x2,x3 = random_point_on_sphere(dimensions=4)
    
    rot = Tensor([[1.0-2.0*(x1*x1 + x2*x2), 2.0*(x0*x1 - x3*x2), 2.0*(x0*x2 + x3*x1)],
                 [2.0*(x0*x1 + x3*x2), 1.0-2.0*(x0*x0 + x2*x2), 2.0*(x1*x2 - x3*x0)],
                 [2.0*(x0*x2 - x3*x1), 2.0*(x1*x2 + x3*x0), 1.0-2.0*(x0*x0 + x1*x1)]])
    
    return rot

def cosine_rule(a,b,c):
    """Solves angle A given sides a,b,c (a is opposite side to angle A)."""
    return math.acos((-a*a + b*b + c*c) / (2.0*b*c))

import random
def _rand_xyz(edge_length=100):
    """Generate a random vector within box of given edge_length."""
    return Vector([(random.random() - 0.5)*edge_length for i in range(3)])

def _rand_rot():
    """Generate a random rotation"""
    from geometry import random_rotation_matrix as rand_rot_matrix
    from vector import Rotation
    return Rotation(rand_rot_matrix()).asQuaternion()

def _rand_one_minus_one():
    
    import math
    n = math.floor(random.random() * 2.0)
    if n == 0: return -1
    else: return +1

def test_rand_rot():
    
    xx = Vector(1.0,0.0,0.0)
    from geometry import apply_quaternion_to_vector
    
    f = open("random_rotation_test.kin", "w")
    print >>f, "@kinemage"
    print >>f, "@dotlist"
    
    for ii in range(10000):
        rand_pt_on_surface_of_sphere = apply_quaternion_to_vector(_rand_rot(), xx)
        x,y,z = rand_pt_on_surface_of_sphere
        print >>f, "%f %f %f" %(x,y,z)

    f.close()

if __name__ == "__main__":

    # TODO: put some test cases here.
    test_rand_rot()
    pass
