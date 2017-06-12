import entities
from mesh import MeshedSphere
from math import pi
from vector import Vector, Rotation, Quaternion
import BEM
from geometry import _rand_rot

origin = Vector(0.0, 0.0, 0.0) # the origin

def create_born_ion_mesh(detail, radius, charge, xyz=Vector(0.0,0.0,0.0)):

    # create diffusing entity
    de = entities.DiffusingEntity()
    de.set_xyz(xyz)
    de.set_rotation(_rand_rot())
    de.set_centre_of_diffusion(origin)
    de.atoms = [entities.Atom(origin, radius)] 
    de.mesh = MeshedSphere(radius=radius, num_subdivides=detail)
    de.mesh._charges = [entities.Charge(Vector(0.0,0.0,0.0), charge, radius=1.0)]

    return de
    
def solve_born_ion(born_ion, Dext=80, Dint=1.0, kappa=0.0):

    BEM.solveElectrostatics([born_ion],kappa,Dext,Dint)
    BEM.solve_energy([born_ion],kappa,Dext,Dint)

    import constants
    print "BORN SOLVATION ENERGY (kJ/mol): ", \
          born_ion.energy*constants.Avogadro
    print "BORN SOLVATION ENERGY (kCal/mol): ", \
          born_ion.energy*constants.Avogadro*0.2388
    
    return born_ion.energy*constants.Avogadro

def born_test():

    # setup
    charge = 1.0
    rad = 3.0
    Dext = 80.0
    Dint = 1.0
    kappa = 1.0 / 7.92729
    
    # analytic solution
    import constants
    premult = -constants.Avogadro * (constants.elementary_charge ** 2) \
            / (1000.0*8.0*pi*constants.epsilon0*constants.Angstroms)
    analytic = premult * (1.0/Dint - 1.0/Dext) * (charge*charge/rad) 
    
    # create a mesh
    born_ion = create_born_ion_mesh(5, rad, charge)

    # solve solvation energy, then coarsen the mesh for next loop iteration
    #while(len(born_ion.mesh.vertices) > 50):

    vertices = len(born_ion.mesh.vertices)
    print "Vertices: %d" %(vertices)

    # solve solvation energy
    print "Solving ... "
    energy = solve_born_ion(born_ion, Dext, Dint, kappa)
    print energy, analytic
        
        #fresults = open("born.txt", "a")
        #print >>fresults, "%d %f %f" %(vertices, energy, analytic)
        #fresults.close()

        ## write a kinemage
        #f = open("born-%d.kin" %(vertices), 'w')
        #print >>f, "@kinemage"
        #born_ion.mesh.kinemage(ref_pt_local=born_ion.centre_of_diffusion,
                               #ref_pt_universe=born_ion.xyz,
                               #rotation=born_ion.rotation,
                               #f=f)
        #f.close()
        
        ## coarsen the mesh
        #born_ion.mesh.gts_coarsen(100)

def born_pair():
    
    rad = 1.0
    Dext = 80.0
    Dint = 1.0
    kappa = 0.0

    d = 10.0 # distance between centres
    
    # create a pair of ions
    ion1 = create_born_ion_mesh(5, rad*3, -1.0, xyz=Vector(-d/2.0,0.0,0.0))
    ion2 = create_born_ion_mesh(4, rad, +1.0, xyz=Vector(d/2.0,0.0,0.0))

    f = open("born_pair-asym.txt","w")
    f.close()

    from geometry import apply_quaternion_to_vector as rotate
    #import tempfile
    #from shutil import move
    from constants import Avogadro
    
    while(len(ion1.mesh.vertices) > 50 or len(ion2.mesh.vertices) > 50):

        for d_int in range(45,120,+5):
            d = 0.1*d_int
            ion1.set_xyz(Vector(-d/2.0,0.0,0.0))
            ion2.set_xyz(Vector(d/2.0,0.0,0.0))
            ion1.set_rotation(_rand_rot())
            ion2.set_rotation(_rand_rot())

            # Born energies for isolated ions
            BEM.solveElectrostatics([ion1],kappa,Dext,Dint)
            BEM.solve_energy([ion1],kappa,Dext,Dint)
            BEM.solveElectrostatics([ion2],kappa,Dext,Dint)
            BEM.solve_energy([ion2],kappa,Dext,Dint)
            isolated_total_e = Avogadro*(ion1.energy + ion2.energy)
            
            # solve energy for system of two interacting ions
            BEM.solveElectrostatics([ion1, ion2],kappa,Dext,Dint)
            BEM.solve_energy([ion1, ion2],kappa,Dext,Dint)
    
            from constants import Avogadro
            #print Avogadro*(ion1.energy + ion2.energy), ion2.mesh.calculate_force(kappa, Dext, si_units=True)
            total_e = Avogadro*(ion1.energy + ion2.energy)
            print total_e, isolated_total_e, total_e - isolated_total_e

            force_1 = ion1.mesh.calculate_force(kappa, Dext, si_units=True)            
            force_2 = ion2.mesh.calculate_force(kappa, Dext, si_units=True)
            rot_force_1 = rotate(ion1.rotation, force_1)
            rot_force_2 = rotate(ion2.rotation, force_2)
            print rot_force_1, rot_force_2
            
            f = open("born-pair-asym.txt","a")
            print >>f, "%f %d %d %f %f %f %f %f %f %f" %(d, 
                                          len(ion1.mesh.vertices), 
                                          len(ion2.mesh.vertices), 
                                          ion1.energy * Avogadro,
                                          ion2.energy * Avogadro,
                                          total_e,
                                          isolated_total_e,
                                          total_e - isolated_total_e,
                                          rot_force_1[0]*1e12, #force in pico-newtons
                                          rot_force_2[0]*1e12) #force in pico-newtons
            f.close()

            total_verts = len(ion1.mesh.vertices) + len(ion2.mesh.vertices)
            
            # write out the converged electrostatic values -- might be useful...
            fh_vals = open("born-pair-asym-solutions-%2.1f-%d.fh" %(d, total_verts), "w")
            for v in ion1.mesh.vertices:
                
                # local and universe coords should coincide in these cases
                uni_v = x,y,z = ion1.convert_local_coordinate(v)
                print >>fh_vals, "%f %f %f %f %f" %(x, y, z, v.f, v.h)

            for v in ion2.mesh.vertices:
                
                # local and universe coords should coincide in these cases
                uni_v = x,y,z = ion2.convert_local_coordinate(v)
                print >>fh_vals, "%f %f %f %f %f" %(x, y, z, v.f, v.h)
                
            fh_vals.close()            

            # write a kinemage
            f = open("born-pair-%2.1f-%d.kin" %(d, total_verts), 'w')
            print >>f, "@kinemage"
            ion1.mesh.kinemage(ref_pt_local=ion1.centre_of_diffusion, ref_pt_universe=ion1.xyz, rotation=ion1.rotation, f=f)
            ion2.mesh.kinemage(ref_pt_local=ion2.centre_of_diffusion, ref_pt_universe=ion2.xyz, rotation=ion2.rotation, f=f)
            f.close()
        
        if len(ion1.mesh.vertices) > 100:
            ion1.mesh.gts_coarsen(100)
        if len(ion2.mesh.vertices) > 100:
            ion2.mesh.gts_coarsen(25)
        
        # coarsen mesh - keep the gts file for later analysis
        #gts_filename = tempfile.mktemp(prefix="ion1-gts_mesh", dir=".")
        #ion1.mesh.gts_coarsen(100, keep_gts_output=gts_filename)
        #move(gts_filename, "ion1-%d.gts" %(len(ion1.mesh.vertices)))
        #gts_filename = tempfile.mktemp(prefix="ion1-gts_mesh", dir=".")
        #ion2.mesh.gts_coarsen(100, keep_gts_output=gts_filename)
        #move(gts_filename, "ion2-%d.gts" %(len(ion2.mesh.vertices)))
        
def born_triple():

    rad = 1.0
    Dext = 80.0
    Dint = 1.0
    kappa = 0.0
    
    d = 8.0 # distance between centres
    
    # create a triple of ions
    ion1 = create_born_ion_mesh(3, rad, -1.0, xyz=Vector(-d/2.0,0.0,0.0))
    ion2 = create_born_ion_mesh(3, rad, +0.0, xyz=Vector(0.0,0.0,0.0))
    ion3 = create_born_ion_mesh(3, rad, +1.0, xyz=Vector(d/2.0,0.0,0.0))

    from geometry import apply_quaternion_to_vector as rotate
    import tempfile
    from shutil import move
    
    ion1.set_rotation(_rand_rot())
    ion2.set_rotation(_rand_rot())
    ion3.set_rotation(_rand_rot())

    BEM.solveElectrostatics([ion1, ion2, ion3],kappa,Dext,Dint)
    BEM.solve_energy([ion1, ion2, ion3],kappa,Dext,Dint)

    #BEM.solveElectrostatics([ion1, ion3],kappa,Dext,Dint)
    #BEM.solve_energy([ion1, ion3],kappa,Dext,Dint)
    
    from constants import Avogadro
    #print Avogadro*(ion1.energy + ion2.energy), ion2.mesh.calculate_force(kappa, Dext, si_units=True)
    total_e = Avogadro*(ion1.energy + ion2.energy + ion3.energy)
    print total_e

    force_1 = ion1.mesh.calculate_force(kappa, Dext, si_units=True)            
    force_2 = ion2.mesh.calculate_force(kappa, Dext, si_units=True)
    force_3 = ion3.mesh.calculate_force(kappa, Dext, si_units=True)
    rot_force_1 = rotate(ion1.rotation, force_1)
    rot_force_2 = rotate(ion2.rotation, force_2)
    rot_force_3 = rotate(ion3.rotation, force_3)
    print rot_force_1, rot_force_2, rot_force_3
    
    f = open("born-triple.txt","a")
    print >>f, "%f %d %d %d %f %f %f %f %f %f %f" %(d, 
                                  len(ion1.mesh.vertices), 
                                  len(ion2.mesh.vertices), 
                                  len(ion3.mesh.vertices), 
                                  ion1.energy * Avogadro,
                                  ion2.energy * Avogadro,
                                  ion3.energy * Avogadro,
                                  total_e,
                                  rot_force_1[0]*1e12, #force in pico-newtons
                                  rot_force_2[0]*1e12, #force in pico-newtons
                                  rot_force_3[0]*1e12) #force in pico-newtons
    f.close()

    total_verts = len(ion1.mesh.vertices) + len(ion2.mesh.vertices) + len(ion3.mesh.vertices)
    
    # write a kinemage
    f = open("born-triple-%2.1f-%d.kin" %(d, total_verts), 'w')
    print >>f, "@kinemage"
    ion1.mesh.kinemage(ref_pt_local=ion1.centre_of_diffusion, ref_pt_universe=ion1.xyz, rotation=ion1.rotation, f=f)
    ion2.mesh.kinemage(ref_pt_local=ion2.centre_of_diffusion, ref_pt_universe=ion2.xyz, rotation=ion2.rotation, f=f)
    ion3.mesh.kinemage(ref_pt_local=ion3.centre_of_diffusion, ref_pt_universe=ion3.xyz, rotation=ion3.rotation, f=f)
    f.close()

if __name__=="__main__":

    import sys
    born_test()
    sys.exit()
    
    
    born_pair()
    
    
    pass
    
