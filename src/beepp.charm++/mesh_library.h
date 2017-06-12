#ifndef __MESH_LIBRARY_H__
#define __MESH_LIBRARY_H__

#include "prerequisites.h"
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

class MeshLibrary : public NodeGroup
{

private:

    std::vector< boost::shared_ptr<Mesh> > mesh_library;

    /// entry methods ///
    void insert(const std::string mesh_tar_filename);

    /// lock for avoiding concurrency nastiness
    CmiNodeLock lock;
    
public:

    /// Constructors ///
    MeshLibrary(std::vector<std::string> filenames);
    MeshLibrary(CkMigrateMessage *msg);
    ~MeshLibrary() {
        CmiDestroyLock(lock);
    }

    // Not an entry method-- clients should get the address of the 'local' mesh library
    // instance and call this directly, as you would a normal C++ object.
    inline const Mesh& get_mesh(int idx) const {
    	//std::cout << "MeshLibrary " << CkMyPe() << " : get_mesh (idx=" << idx << " of " << mesh_library.size() << ")" << std::endl;
        CmiLock(lock);
    	assert(idx < mesh_library.size());
        const Mesh& retval = *(mesh_library[idx]);
        CmiUnlock(lock);
        
    	return retval;
    }
    
    
};

#endif // __MESH_LIBRARY_H__
