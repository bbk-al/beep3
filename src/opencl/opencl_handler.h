// *********************************************************************
// opencl_handler.h
//
// OpenCL handler -- manages a threaded opencl wrapper
//
// *********************************************************************

#ifndef OPENCL_HANDLER_H_
#define OPENCL_HANDLER_H_
#ifdef OPENCL

#include <iosfwd>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <list>
#include <exception>
#include <queue>
#include <CL/opencl.h>
#include "opencl_workblob.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/date_time/posix_time/posix_time_duration.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>


class OpenCL_Handler : boost::noncopyable
{

public:

    //static const size_t RUN_QUEUE_SIZE = GPU_QUEUE_SIZE;
    //static const size_t ALLOC_QUEUE_SIZE = GPU_QUEUE_SIZE;
    
    OpenCL_Handler();
    ~OpenCL_Handler();
    static void context_error(const char *errinfo,
                                const void *private_info,
                                size_t cb,
                                void *user_data);

    inline size_t get_max_malloc_size() const { return max_ocl_malloc_size; }

    long add_work_to_queue(OpenCL_WorkBlob* res_ptr, bool tracked=false);
    long add_work_to_queue_blocking(OpenCL_WorkBlob* res_ptr, bool tracked=false);
    void wait_until_idle();

    bool item_has_finished(long id, bool remove_if_done=true) 
    {
        if (id == 0) {
            std::cerr << "User is asking about OpenCL completion of a work item with qid == 0, but that's not a valid qid: " <<
                         " you probably did not set the tracked=true flag when calling add_work_to_queue(), or discarded/overwrote" <<
                         " the return value. Either way it's wrong so I'm going to throw an exception." << std::endl;
            throw std::exception();
        }
        
        // get a lock on the queue
        boost::mutex::scoped_lock l(fifo_mutex);
        
        std::list<long>::iterator find_it = std::find(completed_tracked.begin(), completed_tracked.end(), id);
        if (find_it == completed_tracked.end())
        {
            return false;
        }
        
        // if the item has completed, remove it from the completion list (to avoid wasting
        // a few bytes!!)
        if (remove_if_done)
        {
            completed_tracked.erase(find_it);
        }
        
        return true;
        
    }
    
    inline size_t pending() {
        size_t retval=0;
        {
            boost::mutex::scoped_lock l(fifo_mutex);
            retval = fifo_queued.size() + fifo_running.size() + fifo_allocated.size(); // + fifo_completed.size();
            //std::cout << "pending: " << retval << std::endl;
        }
        return retval;
    }

private:

    // seems pretty unlikely that we'll encounter a box with more than 16 GPUs...
    static const unsigned int MAX_DEVICES = 16;
    
    // the add_work_to_queue functions return a uid of the work item
    // this can be queried to check for completion of that work item
    long qid;

   // couter to keep track of total number of items run, which triggers
   // a refresh of the OpenCL context after a certain number of items
   // (otherwise we get an OUT_OF_RESOURCES error- not sure why, possibly
   // the OpenCL context doesn't free the memory when I think it does?)
   size_t total_run_count; 

    // OpenCL vars
    cl_context ocl_gpu_context;
    cl_platform_id* ocl_platform_ids;  // Platform ID's
    cl_device_id ocl_devices[MAX_DEVICES];           // OpenCL device
    unsigned int num_devices;
    std::vector<cl_command_queue> ocl_queues;       // OpenCL command queue
    std::vector<cl_command_queue>::iterator ocl_qit; // Command queue iterator

    // OpenCL Device limits
    size_t max_ocl_malloc_size;

    // Thready vars
    // flag to shout "stop!!" at the thread
    volatile bool m_stoprequested;
    volatile bool m_ready;

    boost::mutex fifo_mutex; // controls access to the fifo's
    std::queue< OpenCL_WorkBlob* > fifo_queued;
    std::queue< OpenCL_WorkBlob* > fifo_allocated;
    std::queue< OpenCL_WorkBlob* > fifo_running;
    
    // once upon a time, when I had a separate copy-delete thread I needed
    // the fifo_completed queue to hold completed items.  But I simplified 
    // things so there's only one do_work thread, so this fifo is no longer
    // needed.
    //std::queue< OpenCL_WorkBlob* > fifo_completed; 
    
    std::list<long> completed_tracked; // items which are tracked get their qid put in here on completion

    // the boost thread which spins around processing the fifo queue of work
    boost::thread m_worker_thread;
    boost::thread m_deleter_thread;

    boost::mutex altering_context_mutex;
    
    
    // OpenCL Programs, keyed by the char* of their source code 
    //(which is presumably a static char array somewhere, so unique 
    // to that block of opencl code)
    typedef std::map<const unsigned char*, cl_program> PMap;
    PMap prog_map; 

    // init functions to get a list of all OpenCL devices on the system and associate them
    // with a context, create queues etc.
    void init_opencl();
    void init_context_and_queues();
    void init_queues();
    
    // utility functions to deal with compiling and hash-tabling opencl programs
    cl_program get_program(const unsigned char* source_key);
    static cl_program compile_program(std::string prog_source, cl_context& context, unsigned int num_devices, cl_device_id device_ids[]);
    
    // probably obsolete for our purposes as it is inconvenient to read
    // the opencl source code from file- better to have it as a char array
    static std::string load_program(const std::string& source_filename);
    
    // the function in which the thread spins around
    void do_work();

    // a function to refresh the opencl Context - forces memory deallocation
    void refresh_context();

    // I once had a separate thread running to handle the copy/delete
    // work, but it didn't really help and complicates the threadyness
    // needlessly
    //void do_copy_and_delete();

};
#endif
#endif

