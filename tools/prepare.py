#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pybeep import *

import tarfile
from tarfile import TarFile
import os
import shutil
import tempfile
from os.path import basename, splitext
from pathlib import Path
import time
from centre import calculate_mass


def precalculate_mesh(mtz_filename, mesh_filename, pqr_filename, centre_filename=""):

    # create a tar file to contain mesh
    mtz = tarfile.open(mtz_filename, "w:gz")

    temp_dir = tempfile.mkdtemp(suffix=".beep-prepare")
    
    try:
        
        # figure out if the pqr_filename really is pqr or if have been passed
        # in an xyzqr instead
        base_filename, ext = splitext(basename(pqr_filename)) # remove the extension
        stem_filename = os.path.join(temp_dir, base_filename)
        if ext == ".xyzqr":
            xyzqr_filename = pqr_filename
            pqr_filename = ""
            mtz.add(xyzqr_filename, basename(xyzqr_filename))
        elif ext == ".pqr":
            xyzqr_filename = "%s.xyzqr" %(stem_filename)
            import pqr2xyzqr
            pqr2xyzqr.pqr2xyzqr(pqr_filename, xyzqr_filename)
            mtz.add(pqr_filename, basename(pqr_filename))
            mtz.add(xyzqr_filename, basename(xyzqr_filename))
        else:
            print("Expected either a .xyzqr or .pqr file as input")
            raise BaseException
        
        mtz.add(mesh_filename, basename(mesh_filename))

        # from here the stem filename corresponds to whatever the gts file was called
        stem_filename, ext = splitext(basename(mesh_filename))
        assert(ext == ".gts" or ext == ".off")
        stem_filename = os.path.join(temp_dir, stem_filename)
        
        mesh = Mesh(mesh_filename, xyzqr_filename, True, False, False)

        energies_filename = "%s.energies" %(stem_filename)
        mesh.write_energies(energies_filename)
        mtz.add(energies_filename, basename(energies_filename))
        
        if (centre_filename == ""):
            centre_filename = "%s.centre" %(stem_filename)
            pdbhfile = Path(f"{base_filename}H.pdb")
            if pdbhfile.is_file():
                (m,com) = calculate_mass(str(pdbhfile))
                print(com, file=open(centre_filename,'w'))
                print(f"Using centre of mass from {str(pdbhfile)}")
            else:
                print(mesh.get_centre(), file=open(centre_filename,'w'))
                print(f"Using charge centre")
        mtz.add(centre_filename, basename(centre_filename))
        
        # Check for the mesh2 file
        meshplusfile = Path(f"{base_filename}+{ext}")
        mesh2file = Path(f"{stem_filename}2{ext}")
        e2file = Path(f"{stem_filename}2.energies")
        print(f"meshplusfile {str(meshplusfile)}: {meshplusfile.is_file()}")
        if meshplusfile.is_file():
            meshplus = Mesh(str(meshplusfile), xyzqr_filename, True, False,True)
            st = time.perf_counter()
            mesh2 = Mesh(mesh.create_mesh2(meshplus, True))
            print(f"Meshing time {time.perf_counter()-st}")
            mesh2.write_mesh(str(mesh2file))
            mtz.add(str(mesh2file), basename(str(mesh2file)))
            mesh2.write_energies(str(e2file))
            print(f"{mesh2.num_node_patches} mesh2 patches")
            kin = mesh.kinemage_meshing()
            out = open(f"{base_filename}2.kin",'w')
            print(kin,file=out)
            print(f"output {base_filename}2.kin, {len(kin)}")
            out.close()
            mtz.add(str(e2file), basename(str(e2file)))

        # create xml file
        from xml.dom.minidom import Document
        doc = Document()
        
        # top level element
        contents = doc.createElement("contents")
        doc.appendChild(contents)
        
        def create_xml_node(tag, filename):
            t = doc.createElement(tag)
            txt = doc.createTextNode(basename(filename))
            t.appendChild(txt)
            contents.appendChild(t)

        if (pqr_filename != ""): create_xml_node("pqr",pqr_filename)
        create_xml_node("mesh",mesh_filename)
        if mesh2file.is_file():
            create_xml_node("mesh2",str(mesh2file))
        create_xml_node("xyzqr",xyzqr_filename)
        create_xml_node("centre",centre_filename)
        create_xml_node("energies",energies_filename)
        if e2file.is_file():
            create_xml_node("energies2",str(e2file))
        
        # Print our newly created XML
        xml_filename = os.path.join(temp_dir, "definition.xml")
        definition_xml = open(xml_filename,'w')
        print(doc.toxml(), file=definition_xml)
        definition_xml.close();
        
        # add definition.xml to the .mtz tarfile
        mtz.add(xml_filename, basename(xml_filename))

        energy_f = 0.0
        energy_h = 0.0
        epsilon_ratio = 40.0
        for np in mesh.node_patches:
            np.f = 1.0
            np.h = 1.0
            energy_f += - np.f*np.energy_coefficient_f;
            energy_h += np.h*epsilon_ratio*np.energy_coefficient_h 
        print("energy_f: ", energy_f)
        print("energy_h: ", energy_h)
        print("total energy: ", energy_f+energy_h)
        
    finally:
        
        shutil.rmtree(temp_dir)
    
    return

if __name__=="__main__":
    import sys
    if (len(sys.argv) >= 5): centre_filename = sys.argv[3]
    else: centre_filename = ""
    precalculate_mesh(sys.argv[-1], sys.argv[1], sys.argv[2], centre_filename)
