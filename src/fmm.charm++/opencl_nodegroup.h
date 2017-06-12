/*
 * opencl_nodegroup.h
 *
 *  Created on: 6 Sep 2010
 *      Author: david
 */

#ifndef OPENCL_NODEGROUP_H_
#define OPENCL_NODEGROUP_H_


#include "../common/charm_prereq.h"
#include <charm++.h>
#include <vector>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/thread/mutex.hpp>

class OpenCL_WorkBlob; // fwd decl.
class OpenCL_Handler; // fwd decl.

#ifdef OPENCL
#include "../opencl/opencl_handler.h"
#include "../opencl/opencl_workblob.h"
#endif

class OpenCL_NodeGroup : public CBase_OpenCL_NodeGroup {

public:

	 OpenCL_NodeGroup();
	 ~OpenCL_NodeGroup();
     long add_to_queue(OpenCL_WorkBlob* res_ptr, bool track=false);
     long add_to_queue_blocking(OpenCL_WorkBlob* res_ptr, bool track=false);

	 void wait_for_completion();
	 int pending();
     bool item_has_finished(long id, bool remove_if_done);
     
// 	 inline OpenCL_Handler& get_ocl_handler() { 
//          std::cerr << "Warning: accessing the ocl_handler directly is not thread-safe!" << std::endl;
//          return *ocl_handler; 
//      }

     inline OpenCL_Handler& get_ocl_handler() const { 
         return *ocl_handler; 
     }

private:

	 CmiNodeLock lock;
	 OpenCL_Handler* ocl_handler;

};

#endif /* OPENCL_NODEGROUP_H_ */
