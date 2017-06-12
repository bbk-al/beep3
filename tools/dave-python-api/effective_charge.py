from universe import _rand_xyz
from numpy import matrix, mat, linalg, eye, zeros, array 
from math import exp, pi
from constants import epsilon0, elementary_charge, __Angstroms, kT
from Scientific.Geometry import Vector
from opendx_utils import read_opendx_file, get_potential, real_xyz_to_grid_idx

# the units used by APBS for potentials are dimensionless, multiples of kT/e
kT_on_e = kT / elementary_charge

# the fitting kernel is exp(-kr)/r -- pre_factor takes care of the consts.
# NB: pre_factor includes conversion to units of kT/e (like APBS)
pre_factor = elementary_charge / (4.0*pi*epsilon0*__Angstroms*kT_on_e)

class PotentialValue(object):
    
    """Simple holder for a potential value at an position location."""
    
    def __init__(self, position, pot):
        """Init potential (pot) at a location (position).
        
        These become read-only properties of the PotentialValue class."""
        
        self.__position = position
        self.__potential = pot
        return
    
    @property
    def position(self):
        return self.__position
    
    @property
    def potential(self):
        return self.__potential

def solve(effective_charges, 
          opendx_filename,
          kappa, 
          Dext):

    """Evaluate effective charges to represent a potential field.
    
    Produces a least-squares fit for a set of effective charges to represent
    the 'actual' potential field defined by the real_potential_function,
    assuming that the effective charges are placed in uniform dielectric with
    ionic screening effect controlled by kappa.
    
    Parameters passed in to this function: 
    
    effective_charges -- Charge objects where the location of the charge is
    the location of this effective charge; the magnitude of the charge is
    modified by this function.
    
    opendx_filename -- APBS OpenDX format potential grid to fit against
    
    kappa -- screening constant
    
    Dext -- dielectric constant to use for the effective charge calculations
    (i.e. the dielectric constant of the solvent)
    
    """

    # get apbs results
    print "reading APBS OpenDX file..."
    gridparms, opendx_potentials = read_opendx_file(opendx_filename)
    print "done"
    
    # get potential field to fit against
    print "getting sample points"
    ref_potentials = get_ecm_reference_potentials(gridparms, opendx_potentials)
    #ref_potentials = get_APBS_reference_potentials(gridparms, opendx_potentials)
    test_potentials = ref_potentials # test against fitted points

    # utility function for coulombic potentials (with screening)
    def F(ri, rj):
        r = (ri - rj).length()
        return pre_factor * exp(-kappa * r) / (r*Dext)
    
    # build rhs vector
    print "calculating rhs"
    rhs_vals = []
    for effi in effective_charges:
        val = 0.0
        for integration_point in ref_potentials:
            val += F(effi.position, integration_point.position) * \
                      integration_point.potential
        rhs_vals.append(val)
    rhs = array(rhs_vals)

    # build fitting matrix (see ECM paper)
    print "building matrix"
    matrix = zeros( (len(effective_charges),len(effective_charges)) )
    
    for ii, effi in enumerate(effective_charges):
        for jj, effj in enumerate(effective_charges):

            # integrate over all sample points
            accum = 0.0
            for integration_point in ref_potentials:
                
                accum +=  F(effi.position, integration_point.position) *\
                          F(effj.position, integration_point.position)
            matrix[ii,jj] = accum
    
    eff_charge_mags = linalg.solve(matrix, rhs)

    # diagonalisation
#    eig_vals, eig_vecs = linalg.eig(matrix)
#    bvec = array([0.0 for ii in range( len(effective_charges) )])
#    for i in range(len(effective_charges)):
#        print "eig val %d: %f" %(i, eig_vals[i])
#        for k in range(len(effective_charges)):
#            bvec[i] += eig_vecs[k,i]*rhs_vals[k]
#           
#        bvec[i] /= eig_vals[i]
#        if i>5:
#            bvec[i] = 0.0
#            for k,ch in enumerate(effective_charges):
#                bvec[i] += eig_vecs[i,k]*ch.charge
#    eff_charge_mags = linalg.solve(eig_vecs.T, bvec)

    print eff_charge_mags

    #p = linalg.Permutation(4)
    #LU, P, signum = linalg.LU_decomp(matrix)
    
    #gsl_linalg_LU_decomp (&m.matrix, p, &s)
    #eff_charge_mags = linalg.LU_solve(LU, P, rhs)

    #print "running SVD"
    #(U,V,S) = linalg.SV_decomp(matrix)
    
    #print "solving for least squares fit"
    #eff_charge_mags = linalg.SV_solve(U, V, S, rhs)

    for eff, result in zip(effective_charges, eff_charge_mags):
        eff.set_charge(result)
        
    print "evaluating fit quality"
    denominator = 0.0
    numerator = 0.0
    for integration_point in test_potentials:
        denominator += integration_point.potential * integration_point.potential
        int_pot_accum = 0.0
        for eff in effective_charges:
            int_pot_accum += eff.charge * F(eff.position, integration_point.position)
        numerator += abs(integration_point.potential - int_pot_accum) ** 2
    quality = 1.0 - (numerator / denominator)
    print "fit quality = %9.9f" %(quality)
    print "sum of effective charge magnitudes: %9.9f" \
          %(sum([ch.charge for ch in effective_charges]))

    return effective_charges

def get_random_reference_potentials(real_potential_function):
    """Return set of random reference potentials, computed using
    real_potential_function (supplied by the caller."""
    
    # I could collapse this into a completely unreadable list comprehension.
    # Today I choose not to.
    potentials = []
    ref_pts = spherical_shell_sample(7.0, 13.0)
    print "Number of sample points to fit against: ", len(ref_pts)
    for loc in ref_pts: # HACK!
        # this should return dimensionless value for potential 
        # in multiples of kT/e
        pot = real_potential_function(loc) 

        # put this reference potential in the list
        potentials.append( PotentialValue(loc, pot) )

    return potentials
    
def get_APBS_reference_potentials(grid, opendx_matrix):
    """Returns a reference potentials from APBS grid results."""
    
    potentials = []
    rmin = 2.0+1.9+7.0
    rmax = rmin + 6.0
    r2max = rmax*rmax
    r2min = rmin*rmin
    print grid.x_pts, grid.y_pts, grid.z_pts
    print grid.deltax(), grid.deltay(), grid.deltaz()
    print grid.origin_x, grid.origin_y, grid.origin_z
    xx = grid.origin_x
    for xctr in range(grid.x_pts):
        xx2 = xx*xx
        yy = grid.origin_y
        for yctr in range(grid.y_pts):
            
            yy2 = yy*yy
            zz = grid.origin_z
            for zctr in range(grid.z_pts):
                
                r2 = xx2 + yy2 + (zz*zz)
                if ((r2min < r2) and ( r2 < r2max)):
                    # NB: APBS potential values are in units of kT/e
                    pv = PotentialValue(Vector(xx,yy,zz), 
                                        opendx_matrix[xctr,yctr,zctr])
                    potentials.append( pv )
                
                zz += grid.deltaz()
            yy += grid.deltay()
        xx += grid.deltax()
    print "Number of sample points to fit against: ", len(potentials)
        
    return potentials

def get_ecm_reference_potentials(grid, opendx_matrix):
    """Returns a reference potentials from APBS grid results."""
    
    potentials = []
    print grid.x_pts, grid.y_pts, grid.z_pts
    print grid.deltax(), grid.deltay(), grid.deltaz()
    print grid.origin_x, grid.origin_y, grid.origin_z
    
    ecm_pts = open("ecm_model.pts",'r')
    for line in ecm_pts:
        x,y,z,pot = [float(val) for val in line.split()]
        xidx, yidx, zidx = real_xyz_to_grid_idx(grid, x, y, z)
        dx_val = get_potential(opendx_matrix, grid, x,y,z)        
        pv = PotentialValue(Vector(x,y,z), pot/0.592)
        if dx_val != 0.0:
            print 100.0*(pot - dx_val)/pot
        potentials.append(pv)

    f = open("ecm_pts.kin","w")
    print >>f, "@kinemage"
    print >>f, "@dotlist {sample points from ECM}"
    for pv in potentials:
        x,y,z = pv.position
        print >>f, "%f %f %f" %(x,y,z)
    f.close()        
        
    print "Number of sample points to fit against: ", len(potentials)
    return potentials

def spherical_shell_sample(start_radius, end_radius, num_subdivides=5):
    """Returns a set of points within a spherical shell.
    
    Thickness controlled by the start_radius and end_radius parameters. Tries
    to guess the radial point density based on the level of circumferential
    subdivision."""
    
    point_cloud = []
    
    import _GTSMESH as gts
    import tempfile
    from mesh import get_vertices_from_gts
    import math
    from geometry import _rand_rot, apply_quaternion_to_vector
    
    spherical_file = tempfile.mktemp(suffix=".gts", 
                                     prefix="gts_sphere_", 
                                     dir=".")
    gts.write_spherical_mesh(num_subdivides, spherical_file)
    vertices, num_faces = get_vertices_from_gts(spherical_file)
    
    # cleanup
    import os
    os.remove(spherical_file)
    
    # guess radial spacing from discretization level of the vertices on
    # surface of sphere: area = spacing**2 * sqrt(3)/4.0
    mean_rad = (start_radius+end_radius)/2.0
    mean_area = 4.0*math.pi*(mean_rad**2) / num_faces
    spacing = math.sqrt(4.0*mean_area/math.sqrt(3.0))
    
    shell_thickness = end_radius - start_radius
    num_layers = int(math.ceil(shell_thickness / spacing)) + 1
    spacing_actual = shell_thickness / num_layers

    radius = start_radius
    for ii in range(num_layers):
        #print radius
        #rot = _rand_rot()
        for v in vertices:
            #pt = apply_quaternion_to_vector(rot, v)
            pt = v
            point_cloud.append(pt.normal() * radius)
        
        radius += spacing_actual
    
    return point_cloud

# I used this little utility function to create a PQR file corresponding to
# the test model molecule described in Effective Charges for Macromolecules in
# Solvent (Gabdoulline and Wade 1996)
#def conv_raw():

    #import pdbtools
    #lines = []
    #f_in = open("ecm_model.raw","r")
    #f_out = open("ecm_model.pqr", "w")
    #for raw_line in f_in.readlines():
        #pdb_line = ["ATOM",len(lines)+1,"X","XXX",1]    
        #pdb_line.extend(raw_line.split())
        #print >>f_out, pdbtools.to_pdb_format(pdb_line)
    #f_in.close()
    #f_out.close()
    
def test_spherical_shell():
    
    point_cloud = spherical_shell_sample(7.0,8.0,5)
    
    f = open("spherical_shell_test.kin","w")
    print >>f, "@kinemage"
    print >>f, "@dotlist {spherical shell}"
    for pt in point_cloud:
        x,y,z = pt
        print >>f, "%f %f %f" %(x,y,z)
    f.close()

def test_ecm_model():
    
    filename = "ecm_model"
    
    import pqrtools

    Dext = 78.0
    Dint = 2.0
    kappa = 0.0
    
    from geometry import _rand_rot
    
    getCharges
    # create diffusing entity
    de = entities.DiffusingEntity()
    de.pqr_filename = filename + ".pqr"
    de.set_xyz(origin)
    #de.set_rotation(_rand_rot())
    de.set_charges_from_pqr(filename + ".pqr", set_atoms_from_pqr=True)

    #de.create_mesh(filename)
    de.create_effective_charges(kappa, Dext)

    print de.charges
    print de.effective_charges
        
    return

def pqr2ecm(filename, Dext, kappa):
    
    from pqrtools import getCharges
    from entities import Charge
    
    # filenames
    opendx_filename = "%s.apbs.dx" %(filename)
    pqr_filename = "%s.pqr" %(filename)
    out_filename = "%s.ecm" %(filename)
    
    print "APBS OpenDX file: %s" %(opendx_filename)
    print "PQR file: %s" %(pqr_filename)
    print "Output ECM file: %s" %(out_filename)
    
    # get effective charge locations
    charges = [c for c in getCharges(pqr_filename) 
               if c.charge != 0.0]
    
    # this will assign effective charges to the Charge objects in eff_ch
    solve(charges, 
          opendx_filename,
          kappa, 
          Dext)

    print "completed effective charge calcs"
    
    print "writing ECM (pseudo-PQR format) file"
    pqr_str = Charge.writeChargesAsPQR(charges)
    f = open(out_filename, "w")
    print >>f, "#REMARK   1 Effective Charges for kappa=%f, Dext=%f" \
                      %(kappa, Dext) 
    print >>f, pqr_str
    f.close()

    # print out the charges before/after
    print [c for c in getCharges(pqr_filename)]
    print [c for c in getCharges(out_filename)]
    
    return

if __name__ == "__main__":
    
    #test_spherical_shell()
    #test_ecm_model()
    import sys
    pqr2ecm(sys.argv[1], 78.0, 0.0)
