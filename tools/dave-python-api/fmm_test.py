from _OCTREE import FMM_Octree
from _BEM import Vector
import random
import math

fmm = FMM_Octree(Vector(0,0,0), 1.0, 4)

beta=1e-10 # effectively zero screening
real_pot = 0.0
for i in range(100):

    v = [random.random() - 0.5,random.random() - 0.5,random.random() - 0.5]
    vec = Vector(v)
    r = vec.length()
    #v = [0.1,0,0]
    ch = 1.0
    real_pot += math.exp(-beta*r)/r
    fmm.add_charge(ch, vec)

print "solving..."
fmm.solve(beta)
print "evaluating..."
print fmm.calculate_potential_at_xyz(0.0,0.0,0.0,beta), real_pot
