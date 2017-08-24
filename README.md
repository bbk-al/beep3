# beep3
libBEEP source at v0.3

v0.3 includes short range effects;  v0.2 is electrostatics only but updated from the original 2010 baseline.
The main sources created or modified for v0.3 or 0.2 are:

bem/meshing.h and .cpp - created to allow use of meshing and collision detection functionality from within both Mesh and MeshInstance


Major updates for short-range effects and to support simulator:

beep/beep.h and .cpp - implements the Beep class as a top-level container

bem/mesh.h and .cpp - implements the reference Mesh class

bem/mesh_instance.h and .cpp - implements the MeshInstance class

bem/node_patch.h and .cpp - implements the BasicNodePatch and NodePatch classes


Extensions and corrections:

bem/mesh_tarball.h - MeshTarball class, updated to use Boost::filesystem version 3

bem/spharm.cpp - updated to use latest GSL version (replacing deprecated functions)

common/charge.h - Charge template class

common/math_vector.h - Vector and Quaternion template class

fmm/fmm-octree.h - generalise octree distance calculation

opencl/opencl_handler.h and .cpp - fix major memory leak

pybeep/pybeep.cpp - improve modularity and add/change functions made available to python


In tools, all python scripts were upgraded to python3, and pybeep.py and prepare.py were adapted to support additional meshing.
