/*
* opencl_workblob.h
*
*  Created on: 22 Sep 2010
*      Author: david
*/

#ifndef OPENCL_WORKBLOB_H_
#define OPENCL_WORKBLOB_H_

#include <vector>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/thread/mutex.hpp>

#include <CL/opencl.h>

class OpenCL_WorkBlob : public boost::noncopyable
{

public:

    class AllocationException : public std::exception {};
    class EnqueueException : public std::exception {};

    friend class OpenCL_Handler;
    virtual ~OpenCL_WorkBlob() {}
    virtual const unsigned char* get_opencl_source() { return NULL; }
    virtual void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program)=0;
    virtual void copy_out_results()=0;
    virtual void enqueue(cl_command_queue& ocl_queue)=0;
    virtual bool finished() const=0;

    virtual size_t device_mem_reqd() const=0;

    inline void set_lock(boost::mutex& mx)
    {
        lock_ptr.reset(new boost::mutex::scoped_lock(mx));
    }

protected:
    bool auto_delete;
    long qid;
    
private:
    boost::scoped_ptr<boost::mutex::scoped_lock> lock_ptr;


};

#endif /* OPENCL_WORKBLOB_H_ */

