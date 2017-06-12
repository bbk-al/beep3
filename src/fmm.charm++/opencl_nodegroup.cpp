/*
* opencl_nodegroup.cpp
*
*  Created on: 6 Sep 2010
*      Author: david
*/

#include "vanilla_fmm_prerequisites.h"
#include <charm++.h>
#include "opencl_nodegroup.decl.h"
#include "opencl_nodegroup.h"

OpenCL_NodeGroup::OpenCL_NodeGroup() {
#ifdef OPENCL
    lock = CmiCreateLock();
    ocl_handler = new OpenCL_Handler();
#endif
}

OpenCL_NodeGroup::~OpenCL_NodeGroup() {
#ifdef OPENCL
    CmiDestroyLock(lock);
    delete ocl_handler;
#endif
}

long OpenCL_NodeGroup::add_to_queue(OpenCL_WorkBlob* ptr, bool track)
{
    long retval = 0;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->add_work_to_queue(ptr, track);
    CmiUnlock(lock);
#endif
    return retval;
}

long OpenCL_NodeGroup::add_to_queue_blocking(OpenCL_WorkBlob* ptr, bool track)
{
    long retval = 0;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->add_work_to_queue_blocking(ptr, track);
    CmiUnlock(lock);
#endif
    return retval;
}

void OpenCL_NodeGroup::wait_for_completion()
{
#ifdef OPENCL
    // this will block until the OpenCL handler is empty
    CmiLock(lock);
    ocl_handler->wait_until_idle();
    CmiUnlock(lock);
#endif
    return;
}

int OpenCL_NodeGroup::pending() {

#ifdef OPENCL
    CmiLock(lock);
    int pending = static_cast<int>(ocl_handler->pending());
    CmiUnlock(lock);
    return pending;
#else
    return 0;
#endif
}

bool OpenCL_NodeGroup::item_has_finished(long id, bool remove_if_done) 
{
    bool retval = true;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->item_has_finished(id, remove_if_done);
    CmiUnlock(lock);
#endif
    return retval;

}


#include "opencl_nodegroup.def.h"
