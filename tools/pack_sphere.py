#!/usr/bin/env python3
from random import uniform
from pybeep import Vector
from geometry import _rand_rot, _rand_xyz
import sys

shell_radius = float(sys.argv[1])
n = int(sys.argv[2])
radius = 1.0
diameter = radius * 2
origin = Vector(0,0,0)
centres = []
while (len(centres) < n):
	new_pt = _rand_xyz(shell_radius*2)
	if (new_pt.length() > (shell_radius-radius)):
		continue
	ok = True
	for xx,d in centres:
		if ((new_pt - xx).length2() < radius):
			ok = False
			break
	if ok:
		centres.append((new_pt,diameter))
	if (len(centres) % 1000 == 0):
		print(len(centres),end='\r')

print(f"\r{len(centres)}\n")
output = open(f"packed_sphere_{shell_radius}_{n}.xyzqr",'w')
kin = open(f"packed_sphere_{shell_radius}_{n}.kin",'w')
print("@kinemage", file=kin)
print("@spherelist", file=kin)
for ctr,(pt,d) in enumerate(centres):
	print(f"{pt.x} {pt.y} {pt.z} {uniform(-1,1)} {radius}", file=output)
	print(f"r={radius} {} {pt.x} {pt.y} {pt.z}", file=kin)

