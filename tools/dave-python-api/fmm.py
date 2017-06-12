from Scientific.Geometry import Vector
import octree
#from pygsl import sf as pygsl_sf
from collections import deque
import math
import numpy
import yukawa

class FMM(object):

    def __init__(self, octree):
        """Constructor for FMM."""

        # number of terms might be 9 or 18 depending on the level of precision
        # compiled into the FMM library -- the function query_nterms should
        # return the appropriate value
        nterms = yukawa.qnterms()
        
        # I compiled the library with nterms=9 so it'd better be that :-)
        assert(nterms == 9 or nterms == 18)

        self.nterms = nterms # terms in multipole/exponential expansions

        # the octree holds the charges and their positions
        self.octree = octree

        # assert unit cubic octree
        assert(self.octree.edge_length == 1.0)

    def calc_beta_scaling(self, beta, max_depth):
        """Calculates beta at each level, plus over/underflow scale factors."""

        self.overflow_scale_factor = {}

        # this is what they do in fortran code
        if beta > 1.0:
            self.overflow_scale_factor[1] = 1.0
        else:
            self.overflow_scale_factor[1] = beta

        for level in range(2,max_depth+1):

            # over-under flow scaling
            self.overflow_scale_factor[level] = \
                self.overflow_scale_factor[level-1] / 2.0

    def calc_mp_shift_coefficients(self, beta, max_depth):
        """Calculates the MP z-direction shift coefficients for given level."""

        nterms = self.nterms
        self.mp_shifts = {}

        for level in range(2,max_depth+1):

            r0 = math.sqrt(3.0) / (2.0 ** (level))
            c0 = numpy.zeros((nterms+1,nterms+1,nterms+1),'d')
            info = 0

            c0,info = yukawa.ympshftcoef(self.overflow_scale_factor[level],
                                         beta,
                                         r0,
                                         c0,
                                         info)

            # assert that there weren't any errors
            assert(info == 0)

            self.mp_shifts[level] = c0

        return

    @staticmethod
    def calc_rotation_matrices(nterms):
        """Fast rotation matrices for FMM."""

        carray = numpy.zeros((4*nterms+1,4*nterms+1),'d')
        rdpi2  = numpy.zeros((nterms+1,nterms+1,2*nterms+1),'d')
        rdmpi2 = numpy.zeros((nterms+1,nterms+1,2*nterms+1),'d')
        rdsq3  = numpy.zeros((nterms+1,nterms+1,2*nterms+1),'d')
        rdmsq3 = numpy.zeros((nterms+1,nterms+1,2*nterms+1),'d')

        # call to Huang code
        rdpi2,rdmpi2,rdsq3,rdmsq3 = yukawa.yhrotgen(carray,
                                                    rdpi2,
                                                    rdmpi2,
                                                    rdsq3,
                                                    rdmsq3)

        return [rdpi2,rdmpi2,rdsq3,rdmsq3]

    @staticmethod
    def calc_scale_factors(nterms):
        """Calculate factorial scale factors (used in YHFORMMP)."""

        C = numpy.zeros((nterms+1,nterms+1),'d')
        C = yukawa.yhfrmini(C)
        #print "binomial scale factors (c)", C
        return C

    def get_huang_locals(self, beta):
        """Return local expansions calculated using Huang Fortran FMM code."""

        num_unknowns = len(self.octree.master_list)

        # lists to hold the effective charges
        charges = numpy.zeros(num_unknowns,'d')
        locations = numpy.zeros((3,num_unknowns),'d')

        for idx, wrapped_data in enumerate(self.octree.master_list):

            locations[:,idx] = wrapped_data.position
            charges[idx] = wrapped_data.obj.charge

        nlev = self.octree.max_depth

        # scaled kappa because this is all in a unit cube, rather than the
        # real units of the universe.
        #print "starting FMM calculations (%d unknowns)" %(num_unknowns)
        import bem_fmm
        loc, pot = bem_fmm.run_yukawa_fmm(charges, locations, nlev, beta)

        return loc, pot

    def get_locals_pot_field(self, beta, charge_magnitudes, results):
        
        self.get_local_expansions(beta, 
                                  charge_magnitudes,
                                  results)
        self.evaluate_FMM_potential_at_particles(beta)
        return

    def huang_reference_test(self, beta):
        
        nlev = self.octree.max_depth
        natoms = len(self.octree.master_list)
        pot = numpy.zeros(natoms,'d').T
        field = numpy.zeros((natoms,3), 'd').T
        charge = numpy.array([wrapped_obj.obj.charge 
                              for wrapped_obj in self.octree.master_list]).T
        zat = numpy.array([wrapped_obj.position 
                           for wrapped_obj in self.octree.master_list]).T
    
        import time
        start_huang = time.clock()
    
        # call reference implementation in huang fortran library
        pot, field, ier = yukawa.fmmyuk_uni(beta,
                                            zat,
                                            charge,
                                            pot,
                                            field,
                                            nlev)    
    
        print "Huang FMM: ", time.clock() - start_huang
    
        return pot, field

    def solve(self, beta, calc_local_interactions=True):
        """Solve the electrostatics by fast multipole method.
        
        set calc_local_interactions to False to not calculate the """

        import time
        start_python_fmm = time.clock()

        # scale factors and rotation matrices for multipole expansions
        self.calc_beta_scaling(beta, self.octree.max_depth)        
        self.factorial_scale_factors = self.calc_scale_factors(self.nterms)
        self.rotation_matrices = self.calc_rotation_matrices(self.nterms)

        # pre-calculate shift coefficients for multipoles
        self.calc_mp_shift_coefficients(beta, self.octree.max_depth)
 
        # init holders for expansions
        self.init_expansion_holders()

        # put the charges in the right places
        #for wrapped_obj in self.octree.master_list:
        #    wrapped_obj.charge = wrapped_obj.obj.charge
        
        # upward pass
        # TODO: logging module for messages
        #print "Starting upward pass of tree."
        self.calculate_multipole_expansions(beta)
        self.translate_multipole_expansions()

        # downward pass
        #print "Starting downward pass."
        self.calculate_local_expansions(beta)
        
        node = self.octree.find_leaf_at_position(Vector(0,0,0))
        print "Python Local expansions at ", node.position, "\n", node.local_expansions

        # calculate FMM potential
        #print "Starting FMM potential evaluation."
        self.evaluate_FMM_potential_at_particles(beta)

        # only include the local interactions if they're wanted (under normal
        # circumstances, of course the local interactions are wanted; but in
        # context of BEM actually we just want the potential and field
        # contributions from the far field, near field are integrated
        # separately and more accurately elsewhere)
        if calc_local_interactions:
            # direct (neighbour interactions)
            #print "Calculating nearest neighbour contributions."
            self.add_local_interactions(beta)

        #print "Python FMM: ", time.clock() - start_python_fmm

        # get huang fortran fmm expansions
        #start_huang = time.clock()
        #loc, pot = self.get_huang_locals()
        #print "Huang FMM: ", time.clock() - start_huang

        #print "Done!"
        
        # clear up allocated arrays in the octree
        self.delete_expansion_holders()

        return
    
    def init_expansion_holders(self):
        """Init. numpy arrays for each cube to hold mp and local exps."""

        # process from top of tree downwards
        cube_list = [self.octree]
        nterms = self.nterms

        while(len(cube_list) > 0):

            cube = cube_list.pop()
            cube.mpole = numpy.zeros((nterms+1, nterms+1), 'complex')
            cube.local_expansions = numpy.zeros((nterms+1, nterms+1),
                                                'complex')
            cube.done_upward_pass = False

            # append children to list of cubes to process
            cube_list.extend(cube.children)

    def delete_expansion_holders(self):
        """Destroy numpy arrays in each cube which hold mp and local exps."""

        # process from top of tree downwards
        cube_list = [self.octree]
        nterms = self.nterms

        while(len(cube_list) > 0):

            cube = cube_list.pop()
            del(cube.mpole)
            del(cube.local_expansions)
            del(cube.done_upward_pass)

            # append children to list of cubes to process
            cube_list.extend(cube.children)
            
    def calculate_multipole_expansions(self, beta):
        """For each fine level box, calculate the multipole expansions."""

        # alias
        nterms = self.nterms

        for cube in self.octree.leaves:

            # multipole coefficients matrix
            x0y0z0 = numpy.array(cube.position)
            charge_locations = numpy.array([wrp.position
                                            for wrp in cube._data]).T
            charge_magnitudes = numpy.array([wrp.charge
                                            for wrp in cube._data]).T
            mpole = numpy.zeros((nterms+1, nterms+1), 'complex')
            p_workspace = numpy.zeros((nterms+1,nterms+1),'d')
            overflow_scale_factor = self.overflow_scale_factor[cube.level]
            # empty cubes will cause the fortran code to bail
            if len(charge_magnitudes) > 0:

                cube.mpole = yukawa.yformmp(beta,
                                            x0y0z0,
                                            charge_locations,
                                            charge_magnitudes,
                                            mpole,
                                            overflow_scale_factor,
                                            p_workspace,
                                            self.factorial_scale_factors)
            else:
                cube.mpole = mpole
                
        return

    def translate_multipole_expansions(self):
        """Translate the multipole expansions up through the octree.

        NB: Huang code does this efficiently by use of rotation matrices."""

        nterms = self.nterms

        # get hold of the rotation matrices
        [rdpi2, rdmpi2, rdsq3, rdmsq3] = self.rotation_matrices

        # init queue for upward pass
        queue = deque(self.octree.leaves)
        ctr = 0

        while len(queue) > 0:

            # pop the leftmost cube in the queue
            cube = queue.popleft()

            # skip it if we've already done it
            if cube.done_upward_pass:
                continue
            else:
                cube.done_upward_pass = True

            parent = cube.parent

            # if parent is None then this is the top level cube. Don't bother.
            if parent is None: continue

            ctr += 1
            #print "%d processing a cube..." %(ctr)

            # we shouldn't be dealing with empty cubes!
            #assert(len(cube._data) > 0)

            # after we've done all low level cubes, we'll process the parents,
            # so append them to the queue as we go along (if haven't already
            # done so). Initialise the multipole_expansions dictionary
            # attribute on the parent cubes at the same time.

            mpole_out = numpy.zeros((nterms+1, nterms+1), 'complex')
            shift_coeffs = self.mp_shifts[cube.level]

            # get octant id according to Huang's numbering
            child_octant_id = cube.Huang_octant_id

            if child_octant_id > 4:

                # is in the upper plane: use +125degree y rotation
                rotation_matrix = rdmsq3
                child_octant_id -= 4
            else:

                # lower half plane: use +54degree y rotation
                rotation_matrix = rdsq3

            # allocate workspace for the Huang fortran stuff
            workspace = numpy.zeros((nterms+1, nterms+1), 'complex')

            mpole_out = yukawa.ympshift(child_octant_id,
                                        cube.mpole,
                                        mpole_out,
                                        workspace,
                                        shift_coeffs,
                                        rotation_matrix)

            parent.mpole += mpole_out

            # could check for whether already in queue, but faster to set
            # a 'done flag' and skip it at the start of the loop if
            # already processed
            queue.append(parent)

            # converted all multipole expansions for this child to parent
            # level, pop the next cube off the queue. Note that because we
            # started with a list of all leaf nodes, and we're appending the
            # parents to a deque as we go, we should complete each layer of
            # the tree in turn.

        return

    @staticmethod
    def get_pre_calcs(nterms, scall, betascal):

        nlambs = nterms

        quad_nodes = numpy.zeros((1,nlambs),'d')   # rlams in fortran
        quad_weights = numpy.zeros((1,nlambs),'d') # whts in fortran
        quad_nodes, quad_weights = yukawa.vwts(quad_nodes,
                               quad_weights,
                               nlambs)

        numfour = numpy.zeros((1,nlambs),'d')
        numphys = numpy.zeros((1,nlambs),'d')
        numfour = yukawa.numthetahalf(numfour, nlambs)
        numphys = yukawa.numthetafour(numphys, nlambs)

        for i in range(nlambs):
            test1 = quad_nodes[0,i]
            test2 = math.sqrt(test1*test1 + 2.0*test1*betascal)
            indd = i
            mmax = numphys[0,i]
            for jj in range(i,nlambs):
                if (test2 <= quad_nodes[0,jj]):
                    indd = jj
                    break
                else:
                    mmax = max(numphys[0,jj], mmax)

            numphys[0,i] = max(mmax, numphys[0,indd])

        nexptot = numfour.sum()
        nexptotp = numphys.sum() / 2

        # slightly hacky -- this is size of exponentials array needed for any
        # given number of lambda terms (i.e. p)
        fsize = math.ceil(10**(math.log(nterms+2,2)) / 5.2465)

        fexpe = numpy.zeros((1,fsize),'complex')
        fexpo = numpy.zeros((1,fsize),'complex')
        fexpback = numpy.zeros((1,fsize),'complex')
        fexpe, fexpo, fexpback = yukawa.ymkfexp(numfour,
                                                numphys,
                                                fexpe,
                                                fexpo,
                                                fexpback)

        rlsc = numpy.zeros((nterms+1,nterms+1,nlambs),'d')

        # get scaled legendre functions
        rlsc = yukawa.yrlscini(scall,
                               betascal,
                               rlsc,
                               quad_nodes)

        # get translation operators
        xs = numpy.zeros((3,nexptotp),'complex')
        ys = numpy.zeros((3,nexptotp),'complex')
        zs = numpy.zeros((3,nexptotp),'double')
        xs, ys, zs = yukawa.ymkexps(betascal,
                        quad_nodes,
                        numphys,
                        xs,
                        ys,
                        zs)
        #
        # Multipole to local expansions
        #

        pre_calcs = [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,
                 quad_weights,betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs]

        return pre_calcs

    @staticmethod
    def convert_mp_to_exp(pre_calcs, mpole):
        """Return the multipole expansion for this cube as plane wave."""

        [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,quad_weights,
         betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs] = pre_calcs

        # convert multipole to exponential
        mexpup   = numpy.zeros((1,nexptot),'complex')
        mexpdown = numpy.zeros((1,nexptot),'complex')
        mexpup, mexpdown = yukawa.ympoletoexp(mpole,
                                              numfour,
                                              mexpup,
                                              mexpdown,
                                              rlsc)

        # up
        mexpuphys = numpy.zeros((1,nexptotp),'complex')
        mexpuphys = yukawa.yftophys(mexpup,
                                    quad_nodes,
                                    numfour,
                                    numphys,
                                    mexpuphys,
                                    fexpe,
                                    fexpo)

        # down
        mexpdphys = numpy.zeros((1,nexptotp),'complex')
        mexpdphys = yukawa.yftophys(mexpdown,
                                    quad_nodes,
                                    numfour,
                                    numphys,
                                    mexpdphys,
                                    fexpe,
                                    fexpo)

        return mexpuphys, mexpdphys

    @staticmethod
    def pw_to_local(pre_calcs, wave_up, wave_down, nterms, ytop):

        [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,
         quad_weights,betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs] = pre_calcs

        fmode_up   = numpy.zeros((1,nexptot), 'complex')
        fmode_down = numpy.zeros((1,nexptot), 'complex')
        fmode_up = yukawa.yphystof(fmode_up,
                                   quad_nodes,
                                   numfour,
                                   numphys,
                                   wave_up,
                                   fexpback)
        fmode_down = yukawa.yphystof(fmode_down,
                                     quad_nodes,
                                     numfour,
                                     numphys,
                                     wave_down,
                                     fexpback)

        local_exp = numpy.ones((nterms+1, nterms+1), 'complex')
        local_exp = yukawa.yexptolocal(betascal,
                                       rlsc,
                                       local_exp,
                                       quad_nodes,
                                       quad_weights,
                                       ytop,
                                       numfour,
                                       fmode_up,
                                       fmode_down)

        return local_exp

    def calculate_local_expansions(self, beta):
        """Calculate the local expansions for each level of the tree."""

        [rdpi2, rdmpi2, rdsq3, rdmsq3] = self.rotation_matrices

        # alias: consistency with Huang code
        nterms = nlambs = self.nterms

        total_translation_time = 0.0
        
        # start at the top
        current_level = [self.octree]
        for level in range(1,self.octree.max_depth):

            #print "starting level %d" %(level)
            import time
            s1 = time.clock()

            # pre-calcs for this level
            scall = self.overflow_scale_factor[level+1]
            betascal = beta * (2.0**(-level))
            pre_calcs = self.get_pre_calcs(nterms, scall, betascal)
            nexptotp = pre_calcs[3]
            xs = pre_calcs[12]
            ys = pre_calcs[13]
            zs = pre_calcs[14]

            # shift coeffs for this level
            r0 = math.sqrt(3.0) / (2.0 ** (level+1))
            dc = numpy.zeros((nterms+1,nterms+1,nterms+1),'d')
            info = 0

            dc,info = yukawa.ylcshftcoef(self.overflow_scale_factor[level],
                                         beta,
                                         r0,
                                         dc,
                                         info)

            # assert that there weren't any errors
            assert(info == 0)
            
            #print "done shift coeffs %f" %(time.clock() - s1)
            s1 = time.clock()

            this_levels_children = []
            for cube in current_level:

                # translate local expansion here to children
                for child in cube.children:

                    #
                    # local translation to next level
                    #

                    # get octant id according to Huang's numbering
                    child_octant_id = child.Huang_octant_id

                    # IMPORTANT: This is the opposite rotations from MP
                    # translations because we're pointing from parent to child
                    # rather than child to parent
                    if child_octant_id > 4:
                        # is in the upper plane: use +54degree y rotation
                        rotation_matrix = rdsq3
                        child_octant_id -= 4
                    else:
                        # is in the lower plane: use +125degree y rotation
                        rotation_matrix = rdmsq3

                    # allocate workspace for the Huang fortran stuff
                    local_out = numpy.zeros((nterms+1, nterms+1), 'complex')
                    workspace = numpy.zeros((nterms+1, nterms+1), 'complex')

                    local_out = yukawa.ylcshift(child_octant_id,
                                                cube.local_expansions,
                                                local_out,
                                                workspace,
                                                dc,
                                                rotation_matrix)

                    # store the parentally-translated local expansion
                    child.local_expansions = local_out

                    # init holders for the plane wave expansions
                    child.pw_up = numpy.zeros((1,nexptotp),'complex')
                    child.pw_down = numpy.zeros((1,nexptotp),'complex')
                    child.pw_north = numpy.zeros((1,nexptotp),'complex')
                    child.pw_south = numpy.zeros((1,nexptotp),'complex')
                    child.pw_east = numpy.zeros((1,nexptotp),'complex')
                    child.pw_west = numpy.zeros((1,nexptotp),'complex')
                    
                this_levels_children.extend(cube.children)

            #print "done downward translations  %f" %(time.clock() - s1)
            s1 = time.clock()
                    
### FAST VERSION ###
### DOESN'T ACTUALLY WORK YET THOUGH...###

            ### MP-PW conversions and translations
            #print "MP-PW conversions and translations"
            #for cube in current_level:

                #self.process_up_down(pre_calcs,cube)
                #self.process_north_south(pre_calcs,cube,rdmpi2,nterms)
                #self.process_east_west(pre_calcs,cube,rdpi2,nterms)

            #print "time: %f" %(time.clock() - s1)
                
### END FAST VERSION ###

### SLOW VERSION ###


            #print "starting UP-DOWN" 

            # UP-DOWN
            for cube in this_levels_children:

                up,down = self.convert_mp_to_exp(pre_calcs,cube.mpole)

                if False and level==2 and cube.Huang_octant_id==3 and cube.parent.Huang_octant_id==3:
                    print "level: ", level, " octant ", cube.Huang_octant_id, " cube ", cube.position, " upward:\n", up
                
                for other in cube.uplist:
                    dx,dy,dz = FMM.offset(other, cube)
                    other.pw_up += self.translate_wave_up(up,xs,ys,zs,dx,dy,dz)

                for other in cube.downlist:
                    dx,dy,dz = FMM.offset(other, cube)
                    other.pw_down += self.translate_wave_down(down,xs,ys,zs,dx,dy,dz)

            total_translation_time += (time.clock() - s1)
            #print "up-down %f" %(time.clock() - s1)
            s1 = time.clock()
            #print "starting NORTH-SOUTH"
                     
            ## NORTH-SOUTH
            #for cube in this_levels_children:

                #rot_mpole = numpy.zeros((nterms+1,nterms+1),'complex')
                #workspace = numpy.zeros((nterms+1,nterms+1),'complex')
                #rot_mpole = yukawa.rotztoy(cube.mpole,
                                           #workspace,
                                           #rot_mpole,
                                           #rdmpi2)
                #north,south = self.convert_mp_to_exp(pre_calcs,rot_mpole)
                #for other in cube.northlist:
                    #dy,dz,dx = FMM.offset(other, cube)
                    #other.pw_north += self.translate_wave_up(north,xs,ys,zs,dx,dy,dz)

                #for other in cube.southlist:
                    #dy,dz,dx = FMM.offset(other, cube)
                    #other.pw_south += self.translate_wave_down(south,xs,ys,zs,dx,dy,dz)
            
            #total_translation_time += (time.clock() - s1)                    
            ##print "north-south %f" %(time.clock() - s1)
            #s1 = time.clock()

            ##print "starting EAST-WEST"
                    
            ## EAST-WEST
            #for cube in this_levels_children:
                

                #rot_mpole = numpy.zeros((nterms+1,nterms+1),'complex')
                #rot_mpole = yukawa.rotztox(cube.mpole,
                                           #rot_mpole,
                                           #rdpi2)
                #east,west = self.convert_mp_to_exp(pre_calcs,rot_mpole)
                #for other in cube.eastlist:

                    #dz,dy,dx = FMM.offset(other, cube)
                    #dx=-dx
                    #other.pw_east += self.translate_wave_up(east,xs,ys,zs,dx,dy,dz)

                #for other in cube.westlist:
                    #dz,dy,dx = FMM.offset(other, cube)
                    #dx=-dx
                    #other.pw_west += self.translate_wave_down(west,xs,ys,zs,dx,dy,dz)

            #total_translation_time += (time.clock() - s1)
            ##print "east-west %f" %(time.clock() - s1)
            
                                                              
#### END SLOW VERSION ###

            s1 = time.clock()
            #print "converting plane waves"
            
            # convert exponential waves to local expansions
            for cube in this_levels_children:

                # alias
                ytop = self.factorial_scale_factors

                # UP/DOWN
                cube.local_expansions += self.pw_to_local(pre_calcs,
                                                          cube.pw_up,
                                                          cube.pw_down,
                                                          nterms,
                                                          ytop)

                if level==2 and cube.parent.Huang_octant_id==3 and cube.Huang_octant_id==3:
                    print "level: ", level, " octant ", cube.Huang_octant_id, " cube ", cube.position, " local expansions:\n", cube.local_expansions
                
                # NORTH/SOUTH
                unrotated_local_exp = self.pw_to_local(pre_calcs,
                                                       cube.pw_north,
                                                       cube.pw_south,
                                                       nterms,
                                                       ytop)
                rotated_local_exp = numpy.zeros((nterms+1,nterms+1),'complex')
                workspace = numpy.zeros((nterms+1,nterms+1),'complex')
                rotated_local_exp = yukawa.rotytoz(unrotated_local_exp,
                                                   workspace,
                                                   rotated_local_exp,
                                                   rdpi2)
                cube.local_expansions += rotated_local_exp
                if level==2 and cube.parent.Huang_octant_id==3 and cube.Huang_octant_id==3:
                    print "level: ", level, " octant ", cube.Huang_octant_id, " cube ", cube.position, " local expansions:\n", cube.local_expansions

                
                # EAST/WEST
                unrotated_local_exp = self.pw_to_local(pre_calcs,
                                                       cube.pw_east,
                                                       cube.pw_west,
                                                       nterms,
                                                       ytop)
                rotated_local_exp = numpy.zeros((nterms+1,nterms+1),'complex')
                rotated_local_exp = yukawa.rotztox(unrotated_local_exp,
                                                   rotated_local_exp,
                                                   rdmpi2)
                cube.local_expansions += rotated_local_exp
                if level==2 and cube.parent.Huang_octant_id==3 and cube.Huang_octant_id==3:
                    print "level: ", level, " octant ", cube.Huang_octant_id, " cube ", cube.position, " local expansions:\n", cube.local_expansions

            #
            # End of MP-Local conversion/translation
            #

            # done this discretization level of cubes

            # travel down to the next level of discretization (i.e. process
            # the children of this layer)
            current_level = this_levels_children

            #print "plane wave conversions: %f" %(time.clock() - s1)
            
            #print "done level"
            
        #print "total time spent moving expansions: %f" %(total_translation_time)

        # now we should have done all the local expansions all the way down to
        # the finest level.  Miller time.
        return
    
    def evaluate_FMM_potential_at_point(self, beta, position):
        """Evaluate local expansion at particle positions at finest level of
        the tree."""

        level = self.octree.max_depth
        scale_factor = self.overflow_scale_factor[level]
        
        cube = self.octree.find_leaf_at_position(position)
        
        x0y0z0 = numpy.array(cube.position).T

        point = numpy.array(position).T
        workspace = numpy.zeros((self.nterms+2,self.nterms+2),'d')
        pot = 0.0
        field = numpy.zeros((1,3),'d') # first deriv of pot (field)
        field2 = numpy.zeros((3,3),'d') # derivative of field

        pot, field, field2 = yukawa.bem_calc(beta,
                                             cube.local_expansions,
                                             x0y0z0,
                                             point,
                                             pot,
                                             field,
                                             field2,
                                             scale_factor,
                                             workspace)

        return potential

    def evaluate_FMM_potential_at_particles(self, beta):
        """Evaluate local expansion at particle positions at finest level of
        the tree."""

        level = self.octree.max_depth
        scale_factor = self.overflow_scale_factor[level]

        for cube in self.octree.leaves:
            
            matrix_data = cube.local_expansions
            def Local_nm(n,m):
                if abs(m) > abs(n):
                    return complex(0.0,0.0)
                if m < 0:
                    return ((-1)**m)*Local_nm(n,-m)
                else:
                    return matrix_data[n,m]
            
            x0y0z0 = numpy.array(cube.position).T
            for obj in cube._data:

                point = numpy.array(obj.position).T
                workspace = numpy.zeros((self.nterms+2,self.nterms+2),'d')
                pot = 0.0
                field = numpy.zeros((1,3),'d') # first deriv of pot (field)
                field2 = numpy.zeros((3,3),'d') # derivative of field

                pot, field, field2 = yukawa.bem_calc(beta,
                                                     cube.local_expansions,
                                                     x0y0z0,
                                                     point,
                                                     pot,
                                                     field,
                                                     field2,
                                                     scale_factor,
                                                     workspace)

                obj.potential = pot
                obj.field = field
                obj.field2 = field2

        # that does it for the fast multipole bit of the potential.
        return

    @staticmethod
    def calculate_potential(obj1, obj2, kappa):
        """Calculate potential contribution at obj1 from obj2."""
        
        r = (obj2.position - obj1.position).length()
        potential = obj2.obj.charge * math.exp(-r*kappa) / (r)
        return potential

    def add_local_interactions(self, beta):
        """Add on the screened coulomb potential from local neighbours."""

        # iterate over all the bottom level cubes
        for cube in self.octree.leaves:
            # for each particle position in each cube...
            for wrapped_obj in cube._data:

                accum = 0.0

                # add on the potential from all other points in all other
                # near-neighbour cubes (colleagues of this cube).
                for other_cube in cube.colleagues:
                    for other_obj in other_cube._data:
                        if wrapped_obj is not other_obj:
                            accum += FMM.calculate_potential(wrapped_obj,
                                                             other_obj,
                                                             beta)

                wrapped_obj.potential += accum

        return

    @staticmethod
    def offset(target_cube, ref_cube):
        diff = (target_cube.position - ref_cube.position) / ref_cube.edge_length
        return [int(d) for d in diff]

    @staticmethod
    def get_mpoles(cube):
        """Multipoles expansions for the 8 children of this cube. In correct
        order for use in Jingfang's Fortran code."""
        return [cube.wsd.mpole, cube.esd.mpole,
                cube.wnd.mpole, cube.end.mpole,
                cube.wsu.mpole, cube.esu.mpole,
                cube.wnu.mpole, cube.enu.mpole]
    
    @staticmethod
    def process_up_down(pre_calcs, cube):

        [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,quad_weights,
         betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs] = pre_calcs

        mpole1,mpole2,mpole3,mpole4,mpole5,mpole6,mpole7,mpole8 = \
              FMM.get_mpoles(cube)

        up_total = numpy.zeros((1,nexptotp),'complex')
        up_1234 = numpy.zeros((1,nexptotp),'complex')
        down_total = numpy.zeros((1,nexptotp),'complex')
        down_5678 = numpy.zeros((1,nexptotp),'complex')
        mexpup   = numpy.zeros((1,nexptot),'complex')
        mexpdown = numpy.zeros((1,nexptot),'complex')
        mexpuphys = numpy.zeros((1,nexptotp),'complex')
        mexpdphys = numpy.zeros((1,nexptotp),'complex')

        up_total, up_1234, down_total, down_5678 = \
            yukawa.ymkud_df(mpole1,mpole2,mpole3,mpole4,
                    mpole5,mpole6,mpole7,mpole8,
                    quad_nodes,numfour,numphys,
                    mexpup,mexpdown,mexpuphys,mexpdphys,
                    up_total,up_1234,down_total,down_5678,
                    xs,ys,zs,fexpe,fexpo,rlsc)

        # translate exponential waves to uplist
        for other in cube.wsd.uplist:
            dx,dy,dz = FMM.offset(other, cube.wsd)
            other.pw_up += FMM.translate_wave_up(up_total,xs,ys,zs,dx,dy,dz)

        # translate exponential waves to downlist
        for other in cube.wsu.downlist:
            dx,dy,dz = FMM.offset(other, cube.wsd)
            other.pw_down += FMM.translate_wave_down(down_total,xs,ys,zs,dx,dy,dz)

        # subtract overcounted cubes -- this could be optimized a bit
        for child in cube.children:
            up,down = FMM.convert_mp_to_exp(pre_calcs,child.mpole)
            for other in cube.wsd.uplist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dx,dy,dz = FMM.offset(other, child)
                    other.pw_up -= FMM.translate_wave_up(up,xs,ys,zs,dx,dy,dz)
            for other in cube.wsu.downlist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dx,dy,dz = FMM.offset(other, child)
                    other.pw_down -= FMM.translate_wave_down(down,xs,ys,zs,dx,dy,dz)

        return

    @staticmethod
    def process_north_south(pre_calcs, cube, rdminus, nterms):

        [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,quad_weights,
         betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs] = pre_calcs

        mpole1,mpole2,mpole3,mpole4,mpole5,mpole6,mpole7,mpole8 = \
              FMM.get_mpoles(cube)

        mrotate = numpy.zeros((nterms+1,nterms+1),'complex')
        mwork = numpy.zeros((nterms+1,nterms+1),'complex')

        mexnall = numpy.zeros((1,nexptotp),'complex')
        mexn1256 = numpy.zeros((1,nexptotp),'complex')
        mexn12 = numpy.zeros((1,nexptotp),'complex')
        mexn56 = numpy.zeros((1,nexptotp),'complex')
        mexsall = numpy.zeros((1,nexptotp),'complex')
        mexs3478 = numpy.zeros((1,nexptotp),'complex')
        mexs34 = numpy.zeros((1,nexptotp),'complex')
        mexs78 = numpy.zeros((1,nexptotp),'complex')
        mexpnof = numpy.zeros((1,nexptot),'complex')
        mexpsof = numpy.zeros((1,nexptot),'complex')
        mexpnphys = numpy.zeros((1,nexptotp),'complex')
        mexpsphys = numpy.zeros((1,nexptotp),'complex')

        [mexnall,mexn1256,mexn12,mexn56,
        mexsall,mexs3478,mexs34,mexs78] = \
            yukawa.ymkns_df(mpole1,mpole2,mpole3,mpole4,
                    mpole5,mpole6,mpole7,mpole8,
                    quad_nodes,numfour,numphys,
                    mexpnof,mexpsof,mexpnphys,mexpsphys,
                    mrotate,mwork,rdminus,
                    mexnall,mexn1256,mexn12,mexn56,
                    mexsall,mexs3478,mexs34,mexs78,
                    xs,ys,zs,fexpe,fexpo,rlsc)

        # translate exponential waves to northlist
        for other in cube.wsd.northlist:
            dy,dz,dx = FMM.offset(other, cube.wsd)
            other.pw_north += FMM.translate_wave_up(mexnall,xs,ys,zs,dx,dy,dz)

        # translate exponential waves to southlist
        for other in cube.wnd.southlist:
            dy,dz,dx = FMM.offset(other, cube.wsd)
            other.pw_south += FMM.translate_wave_down(mexsall,xs,ys,zs,dx,dy,dz)

        # subtract overcounted cubes -- this could be optimized a bit
        for child in cube.children:
            rot_mpole = numpy.zeros((nterms+1,nterms+1),'complex')
            workspace = numpy.zeros((nterms+1,nterms+1),'complex')
            rot_mpole = yukawa.rotztoy(child.mpole,
                                       workspace,
                                       rot_mpole,
                                       rdminus)
            north,south = FMM.convert_mp_to_exp(pre_calcs,rot_mpole)
            for other in cube.wsd.northlist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dy,dz,dx = FMM.offset(other, child)
                    other.pw_north -= FMM.translate_wave_up(north,xs,ys,zs,dx,dy,dz)
            for other in cube.wnd.southlist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dy,dz,dx = FMM.offset(other, child)
                    other.pw_south -= FMM.translate_wave_down(south,xs,ys,zs,dx,dy,dz)

        return

    @staticmethod
    def process_east_west(pre_calcs, cube, rdplus, nterms):

        [numfour,numphys,nexptot,nexptotp,rlsc,quad_nodes,quad_weights,
         betascal,scall,fexpe,fexpo,fexpback,xs,ys,zs] = pre_calcs

        mpole1,mpole2,mpole3,mpole4,mpole5,mpole6,mpole7,mpole8 = \
              FMM.get_mpoles(cube)

        mrotate = numpy.zeros((nterms+1,nterms+1),'complex')

        mexeall = numpy.zeros((1,nexptotp),'complex')
        mexe1357 = numpy.zeros((1,nexptotp),'complex')
        mexe13 = numpy.zeros((1,nexptotp),'complex')
        mexe57 = numpy.zeros((1,nexptotp),'complex')
        mexe1 = numpy.zeros((1,nexptotp),'complex')
        mexe3 = numpy.zeros((1,nexptotp),'complex')
        mexe5 = numpy.zeros((1,nexptotp),'complex')
        mexe7 = numpy.zeros((1,nexptotp),'complex')
        mexwall = numpy.zeros((1,nexptotp),'complex')
        mexw2468 = numpy.zeros((1,nexptotp),'complex')
        mexw24 = numpy.zeros((1,nexptotp),'complex')
        mexw68 = numpy.zeros((1,nexptotp),'complex')
        mexw2 = numpy.zeros((1,nexptotp),'complex')
        mexw4 = numpy.zeros((1,nexptotp),'complex')
        mexw6 = numpy.zeros((1,nexptotp),'complex')
        mexw8 = numpy.zeros((1,nexptotp),'complex')
        mexpeof = numpy.zeros((1,nexptot),'complex')
        mexpwof = numpy.zeros((1,nexptot),'complex')
        mexpephys = numpy.zeros((1,nexptotp),'complex')
        mexpwphys = numpy.zeros((1,nexptotp),'complex')

        [mexeall,mexe1357,mexe13,mexe57,mexe1,mexe3,mexe5,mexe7,
        mexwall,mexw2468,mexw24,mexw68,mexw2,mexw4,mexw6,mexw8] = \
                yukawa.ymkew_df(mpole1,mpole2,mpole3,mpole4,
                                mpole5,mpole6,mpole7,mpole8,
                                quad_nodes,numfour,numphys,
                                mexpeof,mexpwof,mexpephys,mexpwphys,
                                mrotate,rdplus,
                                mexeall,mexe1357,mexe13,mexe57,
                                mexe1,mexe3,mexe5,mexe7,
                                mexwall,mexw2468,mexw24,mexw68,
                                mexw2,mexw4,mexw6,mexw8,
                                xs,ys,zs,fexpe,fexpo,rlsc)

        # translate exponential waves to eastlist
        for other in cube.wsd.eastlist:
            dz,dy,dx = FMM.offset(other, cube.wsd)
            other.pw_east += FMM.translate_wave_up(mexeall,xs,ys,zs,-dx,dy,dz)

        # translate exponential waves to westlist
        for other in cube.esd.westlist:
            dz,dy,dx = FMM.offset(other, cube.wsd)
            other.pw_west += FMM.translate_wave_down(mexwall,xs,ys,zs,-dx,dy,dz)
            
        # subtract overcounted cubes -- this could be optimized a bit
        for child in cube.children:
            rot_mpole = numpy.zeros((nterms+1,nterms+1),'complex')
            rot_mpole = yukawa.rotztox(child.mpole,
                                       rot_mpole,
                                       rdplus)
            east,west = FMM.convert_mp_to_exp(pre_calcs,rot_mpole)
            for other in cube.wsd.eastlist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dz,dy,dx = FMM.offset(other, child)
                    other.pw_east -= FMM.translate_wave_up(east,xs,ys,zs,-dx,dy,dz)
            for other in cube.esd.westlist:
                if other in child.colleagues:
                    # subtract the contribution from this child
                    dz,dy,dx = FMM.offset(other, child)
                    other.pw_west -= FMM.translate_wave_down(west,xs,ys,zs,-dx,dy,dz)

        return

    @staticmethod
    def translate_wave_up(wave, xs, ys, zs, dx, dy, dz):

        # assert z-offset is non-zero
        assert(dz > 0)

        mul_array = numpy.ones((1,zs.shape[1]),'complex')
        mul_array *= zs[dz-1,:]

        if dy > 0:
            mul_array *= ys[dy-1,:]
        elif dy < 0:
            mul_array *= ys[(-dy)-1,:].conjugate()

        if dx > 0:
            mul_array *= xs[dx-1,:]
        elif dx < 0:
            mul_array *= xs[(-dx)-1,:].conjugate()

        # element-wise multiplication (not a dot or cross product)
        return wave * mul_array

    @staticmethod
    def translate_wave_down(wave, xs, ys, zs, dx, dy, dz):

        # assert z-offset is non-zero
        assert(dz < 0)

        mul_array = numpy.ones((1,zs.shape[1]),'complex')
        mul_array *= zs[(-dz)-1,:]

        if dy > 0:
            mul_array *= ys[dy-1,:].conjugate()
        elif dy < 0:
            mul_array *= ys[(-dy)-1,:]

        if dx > 0:
            mul_array *= xs[dx-1,:].conjugate()
        elif dx < 0:
            mul_array *= xs[(-dx)-1,:]

        # element-wise multiplication (not a dot or cross product)
        return wave * mul_array

class testObj(object):

    def __init__(self, position, charge):

        self.position = position
        self.charge = charge
        self.potential = complex(0.0,0.0)

# put some good tests here
if __name__ == "__main__":

    import random, time
    import _FMM
    from _BEM import Vector as newVector
    from geometry import random_point_on_sphere

    # consts
    num_test_objs = 1000
    num_direct_tests = 20
    world_size = 1.0
    kappa = 1.0

    nterms = 9 # 6 digit accuracy
    new_fmm = _FMM.FMM_Octree(newVector(0,0,0),world_size,50,nterms)
    
    start = time.clock()

    real_pot = 0.0
    charge = 1.0
    num_charges = 1000
    print "Adding %d charges on surface of sphere..." %(num_charges)
    for xx in xrange(num_charges):
        
        x,y,z = random_point_on_sphere(3)
        v = newVector(x,y,z) / 2.1
        new_fmm.add_charge(v, charge)
        #dist = v.length()
        #real_pot += charge*math.exp(-kappa*dist) / dist
    
    # tree to store the test in
    #tree = octree.Octree(world_size,max_objects_per_node=100)
    #fmm = FMM(tree)

    ## helper function
    #def r(num=3):
        #return [(random.random()*world_size) - (world_size/2.0)
                #for i in range(num)]

    ## insert some test charges into the world
    #real_pot = 0.0
    #for i in range(num_test_objs):
        #charge = random.choice([-1,+1])
        ##charge = 1.0
        #loc = r()
        ##loc = [0.0,-0.499999999,0.0]
        #loc_vector = Vector(loc)
        ##tree.insertObject( testObj(loc_vector,charge) )
        #new_fmm.add_charge(newVector(loc), charge)
        #dist = loc_vector.length()        
        #real_pot += charge*math.exp(-kappa*dist) / dist

    #tree._max_depth = 5
    #tree.optimize()

    #for leaf in tree.leaves:
        #for up in leaf.uplist:
            #assert(leaf in up.downlist)
        #for down in leaf.downlist:
            #assert(leaf in down.uplist)
        #for far_up in leaf.parent.far_uplist:
            #assert(leaf.parent.wsd in far_up.parent.far_downlist)
            ##assert(leaf.parent.wsu in far_up.parent.near_downlist)
        #for far_up in leaf.parent.near_uplist:
            #assert(leaf.parent.wsd in far_up.parent.far_downlist)
            
    # now solve the electrostatics
    #fmm.solve(kappa)
    #print fmm.evaluate_FMM_potential_at_point(kappa, Vector(0,0,0))
    print "Solving..."
    new_fmm.solve(kappa)
    ch = new_fmm.get_charge(0)
    p = ch.get_position()

    real_pot = 0.0
    for i in range(1,new_fmm.get_num_charges()):
        other_ch = new_fmm.get_charge(i)
        dist = (other_ch.get_position() - p).length()
        real_pot += charge*math.exp(-kappa*dist) / dist
        
    
    print "C++ FMM result: ", new_fmm.calculate_potential_at_xyz(ch.get_position(),kappa)
    print "answer should be: ", real_pot
    
    #print "Python FMM time: ", time.clock() - start
    
    ## huang reference implementation (pure fortran)
    #ref_pot, ref_field = fmm.huang_reference_test(kappa)
    
    #for i, obj in enumerate(fmm.octree.master_list[0:num_direct_tests]):
        #pot = 0.0
        #for other_obj in fmm.octree.master_list:
            #if obj is not other_obj:
                #pot += FMM.calculate_potential(obj, other_obj, kappa)
        #error = 100.0 * (abs(pot - obj.potential) / pot)
        #x,y,z = obj.position
        #leaf = fmm.octree.find_leaf_at_position(obj.position)
        #print "FMM potential: %f actual: %f, huang %f, error: %f %% at %f %f %f (%f) (%d, %d)" \
              #%(obj.potential,pot,ref_pot[i],error,x,y,z, obj.position.length(), len(leaf.interaction_list),len(leaf.northlist))

    #print "finished."

