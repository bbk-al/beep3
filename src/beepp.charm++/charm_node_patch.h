/*
* charm_node_patch.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef CHARM_NODE_PATCH_H_
#define CHARM_NODE_PATCH_H_

#include "prerequisites.h"
#include "../bem/node_patch.h"
#include "mesh_library.decl.h"

//class CProxy_MeshLibrary;

class CharmNodePatch : public BasicNodePatch
{

public:

    CharmNodePatch();
    CharmNodePatch(const CProxy_MeshLibrary& mesh_library_proxy,
                   const BasicNodePatch& other, 
                   unsigned int _lib_idx,
                   const Vector& _centre_of_rotation, 
                   const Quaternion& _rotation, 
                   const Vector& _xyz_offset,
                   unsigned int num_quad_points,
                   unsigned int num_qual_points);
    virtual ~CharmNodePatch() {}

    virtual const Mesh& get_ref_mesh() const;
    virtual void pup(PUP::er &p);

    virtual void change_coordinate_frame(QuadList& qps) const
    {
        for (QuadList::iterator it=qps.begin(), end=qps.end(); it != end; ++it)
        {
            it->change_coordinate_frame(centre_of_rotation, rotation, xyz_offset);
        }
    }

private:

    CProxy_MeshLibrary mesh_lib_proxy;
    unsigned int lib_idx;
    Vector centre_of_rotation;
    Quaternion rotation;
    Vector xyz_offset;

};

#endif /* NODE_PATCH_H_ */
