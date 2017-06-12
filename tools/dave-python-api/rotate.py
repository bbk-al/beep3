#!/usr/bin/python
import sys
import os
import random
import tempfile
import stat

def rotate(input_filename, output_filename):
    
    pymol_exec = "pymol -cpq"
    
    pymol_string = """load %s.pdb
rotate [%f,%f,%f],%f
save %s.pdb
""" %(input_filename, random.random(), random.random(), random.random(), random.random()*360, output_filename)

    infile=tempfile.mktemp()
    open(infile,"w").write(pymol_string)
    command="cat " + infile + " | " + pymol_exec

    if os.path.exists(output_filename+".pdb"):
        os.remove(output_filename+".pdb")

    os.system(command)

    while not os.path.exists(output_filename+".pdb"):
	pass
    
    while os.stat(output_filename+".pdb")[stat.ST_SIZE] == 0:
	pass

    print "done"
if __name__ == "__main__":

    filename = sys.argv[1]
    try:
        output_filename = sys.argv[2]
    except IndexError:
        output_filename = filename + "_rotated"
    rotate(filename, output_filename)
   
   
