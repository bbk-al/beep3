#
# mesh_advanced.py
#
# Experimental adaptive coarse quadrature
#
# Author: David Fallaize, University College London, 2008
# E-mail: drf33@cantab.net
#

from mesh import *
import octree
import random 
from Scientific.Geometry import Vector
from pygsl import linalg
from kintools import randomColour


class GTS_Adaptive_Mesh(MeshedObject):

    """GTS format mesh."""

    def __init__(self, gts_filename):
        
        # call the MeshedObject __init__ to set various internal lists.
        super(GTS_Adaptive_Mesh, self).__init__([],[],[])
        
        self._vertices = self.load_gtsfile(gts_filename)
        self._triangles = None

    def kinemage(self, ref_pt_local, ref_pt_universe, rotation=None, f=None, normals=True, charges=True):
        """Produce Kinemage for this Mesh.

        Pass in a file object to write directly to a file, otherwise this will
        return a string."""

        output = []
        
        ## colour list for kinemage
        #for hue,prefix in ((0,"red"),(240,"blue")):
            #for i, saturation in enumerate(range(0,101,1)):
                #colour_name = "%s_%d" %(prefix, i)
                #output.append("@hsvcolor {%s} %d %d 100" %(colour_name, hue, saturation))
        
        output.append("@vectorlist {vertex_normals} color=cyan off")
        for s in self.solution_pts:
            (x1, y1, z1) = s
            (x2, y2, z2) = s + s.new_normal
            output.append("{}P %f %f %f %f %f %f" %(x1, y1, z1, x2, y2, z2))

        output.append("@vectorlist {quad_pts}")            
        for s in self.solution_pts:
            colour = randomColour()                        
            (x1, y1, z1) = s
            (x2, y2, z2) = s + s.new_normal * 2.0
            output.append("{}P %f %f %f %s %f %f %f" %(x1, y1, z1, colour, x2, y2, z2))
            for q in s.quadrature_pts:
                    (x1, y1, z1) = q
                    (x2, y2, z2) = q + q.new_normal
                    output.append("{}P %f %f %f %s %f %f %f" %(x1, y1, z1, colour, x2, y2, z2))
            
        if charges:
            output.append("@spherelist {charges} color=purple off")
            for c in self.charges:
                x,y,z = c.position
                output.append("r=%f {} %f %f %f" %(c.radius,x,y,z))

        # turn the whole thing into a string
        result = "\n".join(output) + "\n"

        # try to write it to a file, otherwise return it as a string
        try:
            f.write(result)
        except AttributeError:
            return result

        return

    def calculate_electric_fields(self):
        """Calculate the Electric Field (E) at each Vertex of the Mesh.

        This function assumes that the f and h attributes of each Vertex has
        already been set (probably by calling set_BEM_results on a NodePatch
        object)."""

        for v in self.solution_pts:
            # create a matrix of normal vectors for the triangles
            # around this vertex
            normals = linalg.zeros((len(v.quadrature_pts),3), linalg.Float)
            rhs = []
            
            for i, q in enumerate(v.quadrature_pts):
            
                # get normal vector for this triangle
                n = q.new_normal
            
                normals[i,0] = n.x()
                normals[i,1] = n.y()
                normals[i,2] = n.z()
            
                rhs.append(v.rep.h)
            if len(v.quadrature_pts) > 2:
                (U,V,S) = linalg.SV_decomp(normals)
                H = Vector(linalg.SV_solve(U, V, S, linalg.array(rhs, linalg.Float)))
            else:
                H = v.rep_normal * v.rep.h
            
            # normalised gives me the unit vector of H
            Hdir = H.normal() # normalised as in length 1.0
            
            #E = Vector(0.0, 0.0, 0.0)
            #for q in v.quadrature_pts:
            
                ## get vertices of the triangle, removing the one we're at
                #verts = t.vertices
                #verts.remove(self)
                #assert(len(verts) == 2) # otherwise this is a v strange triangle!
            
                ## create a linalg (GSL) matrix
                #normals = linalg.zeros((3,3), linalg.Float)
            
                #normals[0,0] = Hdir.x()
                #normals[0,1] = Hdir.y()
                #normals[0,2] = Hdir.z()
            
                ## the .normal() methods here normalise the lengths to 1.0. i.e.
                ## they're calls to the Scientific.Geometry.Vector.normal() method
                #normals[1,0] = (verts[0] - self).normal().x()
                #normals[1,1] = (verts[0] - self).normal().y()
                #normals[1,2] = (verts[0] - self).normal().z()
            
                #normals[2,0] = (verts[1] - self).normal().x()
                #normals[2,1] = (verts[1] - self).normal().y()
                #normals[2,2] = (verts[1] - self).normal().z()
            
                #rhs = []
                #rhs.append(H.length())
                #rhs.append( (verts[0].f - self.f) / (verts[0] - self).length() )
                #rhs.append( (verts[1].f - self.f) / (verts[1] - self).length() )
            
                #(U,V,S) = linalg.SV_decomp(normals)
                #E += Vector(linalg.SV_solve(U, V, S, linalg.array(rhs, linalg.Float)))
            
            #self.E = E / len(self.triangles)
            v.rep.h = H * v.rep_normal # bit of a sneaky hack            
    
    def set_approximate_quadrature_points(self, simpler_mesh=None):

        # open the simpler mesh to get approximate solution locations
        if simpler_mesh is not None:
            simpler_pts = self.load_gtsfile(simpler_mesh)
            for s in simpler_pts:
                s.quadrature_pts=[]
            num_solution_pts = len(simpler_pts)
            
        # using approximate locations group the real vertices
        #num_solution_pts = int(0.005 * len(self._vertices))
        
        #num_random_solution_pts = 200

        # choose a random selection of vertices to be solution points
        # the remainder of vertices in the mesh are quadrature points
        #solution_idxs = [random.randint(num_random_solution_pts,len(self.vertices)-1) 
        #                 for i in range(num_random_solution_pts)]
        #solution_idxs.extend([i for i in range(num_random_solution_pts)])
        #num_solution_pts = len(solution_idxs)

        #solutions = [self.vertices[idx] for idx in solution_idxs]
        #for s in solutions:
        #    s.quadrature_pts = [s] # init. holder for quadrature points

        all_pts = self.vertices[:]
        all_pts.extend(simpler_pts)
            
        mins = [None,None,None]
        maxs = [None,None,None]
        
        for v in all_pts: #self.vertices:
            for coord in (0,1,2): # loop over x,y,z
                val = v[coord]
                if mins[coord] is None or val < mins[coord]:
                    mins[coord] = val
                if maxs[coord] is None or val > maxs[coord]:
                    maxs[coord] = val
        size = (Vector(maxs) - Vector(mins)).length()

        # init an octree for reasonably fast lookup of nearest solution point
        print "building octree"
        tree = octree.Octree(worldSize=size, worldCentre=(Vector(maxs)+Vector(mins))/2.0, max_objects_per_node=1)
        for v in simpler_pts: #solutions:
            tree.insertObject(v, v)
        
        print "starting classification of quad points (%d sol pts vs %d quad pts)" %(num_solution_pts, len(self.vertices))

        # iterate over quadrature points, find octree leaf node in which point lies
        # then traverse upwards until we encompass at least one solution point.
        # Then we need to simply check whether the quadrature point is nearest to that pt,
        # compared to all other possibilities in the vicinity (list of colleagues)
        for quad_idx in range(len(self.vertices)):
            #if quad_idx in solution_idxs: continue
            q = self.vertices[quad_idx]
            leaf = tree.find_leaf_at_position(q)
            while (len(leaf._data) == 0):
                leaf = leaf.parent

            # get closest solution pt
            possibilities = []
            for chk in leaf.colleagues:
                possibilities.extend(chk._data)
            if len(possibilities)==1:
                closest_sol_pt = possibilities[0]
            else:
                dists = [((p.obj-q).length(),p) for p in possibilities]
                min_dist, closest_sol_pt = min(dists)
            
            closest_sol_pt.quadrature_pts.append(q)
        print "done"

        #for s in solutions:
        #    print "sol pt: ",s, sum([q.dynamic_area for q in s.quadrature_pts]) + s.dynamic_area
        #    for q in s.quadrature_pts:
        #        print q,q.dynamic_area
        print "total s.a: ", sum([v.dynamic_area for v in self.vertices]) 

        self.solution_pts = [s for s in simpler_pts if len(s.quadrature_pts)>0]#solutions
        
        for s in self.solution_pts:
            
            area = sum([q.dynamic_area 
                        for q in s.quadrature_pts])
            blob_normal = Vector(0.0,0.0,0.0)
            blob_centroid = Vector(0.0,0.0,0.0)
            total_dist = 0.0
            for q in s.quadrature_pts:
                total_dist += q.length()
                blob_centroid += q*q.dynamic_area/area
                blob_normal += q.new_normal*q.dynamic_area/area
            #blob_centroid /= len(s.quadrature_pts)
            #blob_normal = blob_normal.normal()
            
            dists = [((q-blob_centroid).length(),q) for q in s.quadrature_pts]
            s.rep = min(dists)[1]
            s.rep_normal = s.rep.new_normal()
            print s.rep_normal
        
        return 

    def enforce_no_charge_clashes(self, atomlist, keep_gts_output=None):
        something_moved = False
        
        print "fixing vertex/charge clashes"
        for v in self._vertices:
            
            while True:
                unclean = False
                for at in atomlist:
                    
                    while True:
                        
                        dist = (v - at.position).length()
                        
                        if dist >= at.radius: 
                            break
                        
                        unclean = True
                        something_moved = True
                        #print dist, at.radius
                        #print "vertex within atom radius!"
                        direction = ((v-at.position)/dist + v.new_normal)/2.0
                        offset = v + direction*1e-3
                        #print offset - v

                        v.__setstate__(offset)
                        for t in v.triangles:
                            t._force_recalc()

                # loop until all vertices missing all atoms
                if not unclean:
                    break
        if keep_gts_output is not None:
            self.write_gts(keep_gts_output)
            
        return something_moved
    
    @staticmethod
    def load_gtsfile(filename):
        
        from gts_utils import get_next_line_not_comment    
        f = open(filename,"r")    
        
        vertices = []
        try:
        
            # get the first line, ignoring comment lines (starting with # or !)
            first_line = get_next_line_not_comment(f)
                
            # number of vertices, edges and faces
            nv, ne, nf = [int(xx) for xx in first_line.split()[:3]]
            #print nv, ne, nf
            while (len(vertices) < nv):
                
                # convert to Vertex object
                vx, vy, vz = [float(xx) for xx in get_next_line_not_comment(f).split()]
                new_vertex = geometry.Vertex(vx, vy, vz)
                new_vertex.normal_accum = Vector(0.0,0.0,0.0)
                new_vertex.norm_ctr = 0
                new_vertex.dynamic_area = 0.0
                vertices.append(new_vertex)
    
            # edges are defined as an ordered pair of vertices
            edges = []
            for i in range(ne):
                edge_verts = get_next_line_not_comment(f).split()
                edges.append( set([int(edge_verts[0])-1, int(edge_verts[1])-1]) )
            
            # triangle defs come after the edge defs
            tcentres = []
            for i in range(nf):
                e1, e2, e3 = [edges[int(xx)-1] for xx in get_next_line_not_comment(f).split()]
                
                # get the 3 vertices, in order
                v1_idx = e1.intersection(e2).pop()
                v2_idx = e2.intersection(e3).pop()
                v3_idx = e3.intersection(e1).pop()
                
                v1 = vertices[v1_idx]
                v2 = vertices[v2_idx]
                v3 = vertices[v3_idx]
        
                new_t = geometry.Triangle(v1, v2, v3)
                quad = geometry.Vertex(new_t.centre)
                quad.dynamic_area = new_t.area
                quad.new_normal = new_t.normal
                tcentres.append(quad)
                
                #v1.normal_accum += new_t.normal * (new_t.area / (new_t.centre-v1).length())
                #v2.normal_accum += new_t.normal * (new_t.area / (new_t.centre-v2).length())
                #v3.normal_accum += new_t.normal * (new_t.area / (new_t.centre-v3).length())
                
                #v1.norm_ctr += 1
                #v2.norm_ctr += 1
                #v3.norm_ctr += 1
    
            #for v in vertices:
                #v.new_normal = v.normal_accum.normal()
            #popularity = [(v.norm_ctr,v) for v in vertices]
            #popularity.sort()
            #popularity.reverse()
            #vertices = [v for p,v in popularity]

            vertices = tcentres
            
        except IOError:
            print "Bad GTS file."
            vertices = []
    
        finally:
            f.close()
        
        return vertices    

class SIMS_Adaptive_Mesh(GTS_Adaptive_Mesh):
    
    def __init__(self, filename):
        
        # call the MeshedObject __init__ to set various internal lists.
        #super(SIMS_Adaptive_Mesh, self).__init__([],[],[])
        
        self._vertices = self.load_sims_file(filename)
        self._triangles = None
        self._charges = None
        self._patches = None

    @staticmethod
    def load_sims_file(filename):

        f = open(filename, "r")
        
        vertices = []
        for line in f:
            if line.strip()[0] == "#": continue
            area, x,y,z, nx,ny,nz = [float(val) for val in line.split()[-7:]]
            v = geometry.Vertex(Vector(x,y,z))
            v.new_normal = Vector(nx, ny, nz).normal()
            v.dynamic_area = area
            vertices.append(v)
            
        return vertices

    def set_approximate_quadrature_points(self, simpler_mesh=None):

        num_random_solution_pts = len(self.vertices)

        # choose a random selection of vertices to be solution points
        # the remainder of vertices in the mesh are quadrature points
        solution_idxs = [random.randint(0,len(self.vertices)-1)
                         for i in range(num_random_solution_pts)]
        #solution_idxs.extend([i for i in range(num_random_solution_pts)])
        num_solution_pts = len(solution_idxs)

        solutions = [self.vertices[idx] for idx in solution_idxs]
        for s in solutions:
            s.quadrature_pts = [s] # init. holder for quadrature points

        mins = [None,None,None]
        maxs = [None,None,None]
        
        for v in self.vertices:
            for coord in (0,1,2): # loop over x,y,z
                val = v[coord]
                if mins[coord] is None or val < mins[coord]:
                    mins[coord] = val
                if maxs[coord] is None or val > maxs[coord]:
                    maxs[coord] = val
        size = (Vector(maxs) - Vector(mins)).length()

        # init an octree for reasonably fast lookup of nearest solution point
        print "building octree"
        tree = octree.Octree(worldSize=size, worldCentre=(Vector(maxs)+Vector(mins))/2.0, max_objects_per_node=1)
        for v in solutions:
            tree.insertObject(v, v)
        
        print "starting classification of quad points (%d sol pts vs %d quad pts)" %(num_solution_pts, len(self.vertices))

        # iterate over quadrature points, find octree leaf node in which point lies
        # then traverse upwards until we encompass at least one solution point.
        # Then we need to simply check whether the quadrature point is nearest to that pt,
        # compared to all other possibilities in the vicinity (list of colleagues)
        for quad_idx in range(len(self.vertices)):
            if quad_idx in solution_idxs: continue
            q = self.vertices[quad_idx]
            leaf = tree.find_leaf_at_position(q)
            while (len(leaf._data) == 0):
                leaf = leaf.parent

            # get closest solution pt
            possibilities = []
            for chk in leaf.colleagues:
                possibilities.extend(chk._data)
            if len(possibilities)==1:
                closest_sol_pt = possibilities[0]
            else:
                dists = [((p.obj-q).length(),p) for p in possibilities]
                min_dist, closest_sol_pt = min(dists)
            
            closest_sol_pt.quadrature_pts.append(q)
        print "done"

        #for s in solutions:
        #    print "sol pt: ",s, sum([q.dynamic_area for q in s.quadrature_pts]) + s.dynamic_area
        #    for q in s.quadrature_pts:
        #        print q,q.dynamic_area
        print "total s.a: ", sum([v.dynamic_area for v in self.vertices]) 

        self.solution_pts = []
        for s in [s for s in solutions if len(s.quadrature_pts)>0]:
            print s, s.new_normal
            area = sum([q.dynamic_area 
                        for q in s.quadrature_pts])
            if area == 0.0: continue
            blob_normal = Vector(0.0,0.0,0.0)
            blob_centroid = Vector(0.0,0.0,0.0)
            for q in s.quadrature_pts:
                blob_centroid += q*(q.dynamic_area/area)
                blob_normal += q.new_normal*(q.dynamic_area/area)
            blob_normal = blob_normal.normal()
            
            #dists = [((q-blob_centroid).length(),q) for q in s.quadrature_pts]
            #s.rep = min(dists)[1]
            #s.rep = blob_centroid
            #s.rep_normal = blob_normal
            
            new_pt = geometry.Vertex(blob_centroid)
            new_pt.new_normal = blob_normal
            new_pt.dynamic_area = area
            new_pt.quadrature_pts = s.quadrature_pts[:]
            self.solution_pts.append(new_pt)

        return     
    
if __name__=="__main__":

    def born_ion():
        
        import entities
        from vector import Vector, Rotation, Quaternion
    
        rad = 1.0
        detail = 2
        Dext = 80.0
        Dint = 1.0
        charge = 1.0
        kappa = 0.0
        
        origin = Vector(0.0, 0.0, 0.0) # the origin
    
        from math import pi
    
        # create diffusing entity
        de = entities.DiffusingEntity()
        import universe
        de.set_xyz(Vector(0.0,0.0,0.0))
        de.atoms = [entities.Atom(origin, rad)] 
        de.mesh = GTS_Adaptive_Mesh("sphere.gts")
        for v in de.mesh._vertices:
            v.__setstate__(v*rad)
            v.dynamic_area *= rad*rad
        de.mesh.set_approximate_quadrature_points("sphere-42.gts")
        de.mesh._charges = [entities.Charge(origin, charge, radius=rad)]
        
        import constants
        premult = -constants.Avogadro * (constants.elementary_charge ** 2) \
                / (1000.0*8.0*math.pi*constants.epsilon0*constants.__Angstroms)
        #analytic = -700.7725 * (1.0 - 1.0/Dext) * (charge*charge/rad) #Nathan's
        analytic = premult * (1.0/Dint - 1.0/Dext) * (charge*charge/rad) 
        
        import BEM
        BEM.solveElectrostatics([de],kappa,Dext,Dint)
        BEM.solve_energy([de],kappa,Dext,Dint)

        print "BORN SOLVATION ENERGY (kJ/mol): ", \
              de.energy*constants.Avogadro , analytic, 100.0 * (de.energy*constants.Avogadro - analytic) / analytic
        print "BORN SOLVATION ENERGY (kCal/mol): ", \
              de.energy*constants.Avogadro*0.2388, analytic * 0.2388

        # write a kinemage
        f = open("born.kin", 'w')
        print >>f, "@kinemage"
        de.mesh.kinemage(ref_pt_local=origin, ref_pt_universe=origin,f=f)
        f.close()                
        
        #print de.mesh.calculate_force(kappa,Dext,si_units=True) #* constants.kT * 0.2388 * constants.Avogadro
        #print de2.mesh.calculate_force(kappa,Dext,si_units=True)
        
        #de.mesh.write_gts("born_ion.gts")
        #de.mesh.gts_coarsen(0.9)
        
        return    
    
    def some_molecule(filename):
        
        import entities
        from vector import Vector, Rotation, Quaternion
    
        Dext = 80.0
        Dint = 2.0
        kappa = 0.0
        
        origin = Vector(0.0, 0.0, 0.0) # the origin
        
        #from Scientific.IO.PDB import Structure
        #struct = Structure(filename + ".pdb")
        
        from geometry import _rand_rot
        
        # create diffusing entity
        de = entities.DiffusingEntity()
        de.pqr_filename = filename + ".pqr"
        de.set_rotation(Rotation(Vector(1.0,0.0,0.0),0.0).asQuaternion())
        de.set_charges_from_pqr(filename + ".pqr", set_atoms_from_pqr=True)

        # gts mesh must already exist
        de.create_mesh(filename)
        
        centre = Vector(0.0,0.0,0.0)
        for atom in de.atoms:
            centre += atom.position
        de.set_centre_of_diffusion( centre / len(de.atoms) )
        de.set_xyz(de.centre_of_diffusion)
        
        f = open("preview.kin",'w')
        print >>f, "@kinemage"
        de.mesh.kinemage(ref_pt_local=centre, ref_pt_universe=de.xyz, f=f)
        f.close()

        import BEM
        import tempfile
        from shutil import move
        while(len(de.mesh.vertices) > 20):
        
            v = len(de.mesh.vertices)
                    
            BEM.solveElectrostatics([de],kappa,Dext,Dint)
            BEM.solve_energy([de],kappa,Dext,Dint)
    
            import constants
            print "SOLVATION ENERGY (kJ/mol): ", de.energy*constants.Avogadro
            print "SOLVATION ENERGY (kCal/mol): " , de.energy*constants.Avogadro*0.2388

            fresults = open("%s.txt" %(filename), "a")
            print >>fresults, "%d %f" %(v, de.energy*constants.Avogadro)
            fresults.close()

            f = open("%s-%d.kin" %(filename,v),'w')
            print >>f, "@kinemage"
            de.mesh.kinemage(ref_pt_local=de.centre_of_diffusion, ref_pt_universe=de.xyz, rotation=de.rotation, f=f)
            f.close()            
            
            # write out the converged electrostatic values -- might be useful...
            fh_vals = open("%s-%d.fh"%(filename, len(de.mesh.vertices)), "w")
            for v in de.mesh.vertices:
                
                # local and universe coords should coincide in these cases
                uni_v = x,y,z = de.convert_local_coordinate(v)
                assert((uni_v-v).length() < 1e-12)
                print >>fh_vals, "%f %f %f %9.9e %9.9e" %(x, y, z, v.f, v.h)
            fh_vals.close()
            
            # coarsen mesh - keep the gts file for later analysis
            gts_filename = tempfile.mktemp(prefix="%s-gts_mesh" %(filename), dir=".")
            de.mesh.gts_coarsen(200)
            while de.mesh.enforce_no_charge_clashes(de.atoms, keep_gts_output=gts_filename):
                de.mesh.gts_coarsen(0)
            move(gts_filename, "%s-%d.gts" %(filename, len(de.mesh.vertices)))
            
        return
    
    #import cProfile
    #cProfile.run('born_ion()','born_ion_profile')
    
    #born_ion()
    some_molecule("ubiquitin/ubiquitin")
    #some_molecule("1MAH-ache")
    #some_molecule("2CGA-mono")
    #some_molecule("ecm_model")
    #some_molecule("2OZO")
