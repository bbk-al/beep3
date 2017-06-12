/*
 * opencl_nodegroup.h
 *
 *  Created on: 6 Sep 2010
 *      Author: david
 */

#ifndef OPENCL_NODEGROUP_H_
#define OPENCL_NODEGROUP_H_


#include "prerequisites.h"
#include <vector>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/thread/mutex.hpp>
#include "../bem/local_integrations.h"

class OpenCL_WorkBlob; // fwd decl.
class OpenCL_Handler; // fwd decl.

//class OpenCL_CkCB : public OpenCL_WorkBlob
//{
//
//public:
//	OpenCL_CkCB(CkCallbackResumeThread& cb) : callback(cb) {
//		OpenCL_WorkBlob::auto_delete = true;
//	}
//	~OpenCL_CkCB() { callback.send(); }
//	std::string opencl_source_filename() { return "opencl_fmm_kernels.cl"; }
//	void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program) { return; }
//	void copy_out_results() { return; }
//	void enqueue(cl_command_queue& ocl_queue) { return; }
//	bool finished() const { return true; }
//private:
//	CkCallbackResumeThread callback;
//};

class OpenCL_NodeGroup : public CBase_OpenCL_NodeGroup {

public:

	 OpenCL_NodeGroup(unsigned int num_patches);
	 ~OpenCL_NodeGroup();
	 void run_bem(CkArrayIndexOctreeIndexer idxer, double kappa);
	 void precalc_bem(CkArrayIndexOctreeIndexer idxer, double kappa);
     long add_to_queue(OpenCL_WorkBlob* res_ptr, bool track=false);
     long add_to_queue_blocking(OpenCL_WorkBlob* res_ptr, bool track=false);

	 void collate_bem_results(CkCallback cb);
	 void wait_for_completion();
	 int pending();
	 bool item_has_finished(long id, bool remove_if_done); 
     inline const OpenCL_Handler& get_ocl_handler() const { return *ocl_handler; }

private:

     bool m_ready;
     
	 CmiNodeLock lock;
	 OpenCL_Handler* ocl_handler;

	 unsigned int total_num_patches;
	 boost::scoped_array<double> results;
	 std::vector<LintArray_Size> explicit_integrations;
	 boost::mutex results_mutex;

};

#endif /* OPENCL_NODEGROUP_H_ */
