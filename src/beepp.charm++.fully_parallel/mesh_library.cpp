#include "prerequisites.h"

#include "mesh_library.decl.h"
#include "mesh_library.h"

MeshLibrary::MeshLibrary(std::vector<std::string> filenames) {

	std::cout << "MeshLibrary on processor " << CkMyPe() << " address: " << this << std::endl;
    lock = CmiCreateLock(); 
    CmiLock(lock);
	// for each filename passed in (which corresponds to a gzipped-tar'ed mesh in .mtz file,
	// create an entry in this 'mesh library' so everyone knows what is what.
	for (std::vector<std::string>::const_iterator it=filenames.begin(), end=filenames.end(); it != end; ++it)
	{
		insert(*it);
	}
    CmiUnlock(lock);
	//std::cout << "Added " << mesh_library.size() << " meshes to MeshLibrary on Processor " << CkMyPe() << std::endl;
}

void MeshLibrary::insert(const std::string mesh_tar_filename)
{
	// unzip the mesh tarball and instantiate a mesh into the mesh library
	boost::shared_ptr<Mesh> mesh_ptr(new Mesh(mesh_tar_filename));
	mesh_library.push_back(mesh_ptr);
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
MeshLibrary::MeshLibrary(CkMigrateMessage *msg) { }

#include "mesh_library.def.h"
