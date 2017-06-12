#!/usr/bin/python

base_filename="1D3Z"

counter = 1
fout = open("%s-%d.pdb" %(base_filename, counter), 'w')
for line in open(base_filename+".pdb", 'r').readlines():
    if (line[:6] == "ENDMDL"):
        fout.close()
        counter += 1
        fout = open("%s-%d.pdb" %(base_filename, counter), 'w')
    if (line[:4] == "ATOM"):
        print >>fout, line.strip()
fout.close()
