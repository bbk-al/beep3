#!/usr/bin/python
from random import uniform
from pybeep import Vector
from geometry import _rand_rot, _rand_xyz
import sys
from string import atoi, atof

shell_radius = atof(sys.argv[1])
n = atoi(sys.argv[2])
radius = 1.0
diameter = radius * 2
origin = Vector(0,0,0)
centres = []
while (len(centres) < n):
    new_pt = _rand_xyz(shell_radius*2)
    if (new_pt.length() > (shell_radius-radius)): continue
    ok = True
    for xx,d in centres:
        if ((new_pt - xx).length2() < 1):
            ok = False
            break
    if ok:
        centres.append((new_pt,diameter))
    if (len(centres) % 1000 == 0): print len(centres)

output = open("packed_sphere_%d_%d.xyzqr" %(shell_radius, n),'w')
kin = open("packed_sphere_%d_%d.kin" %(shell_radius, n),'w')
print >>kin, "@kinemage"
print >>kin, "@spherelist"
for ctr,(pt,d) in enumerate(centres):
    print >>output, "%f %f %f %f %f" %(pt.x, pt.y, pt.z, uniform(-1,1), 1.)
    print >>kin, "r=1.0 {} %f %f %f" %(pt.x, pt.y, pt.z)

