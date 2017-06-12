// *********************************************************************
// opencl_handler.cpp
//
// OpenCL handler -- manages a threaded opencl wrapper
//
// *********************************************************************
#ifdef OPENCL
#include "opencl_handler.h"

#include <iostream>
#include <string>
#include <cstring> // i need strstr
#include <fstream>
#include <boost/bind.hpp>
#include <boost/date_time/posix_time/posix_time_duration.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

OpenCL_Handler::OpenCL_Handler() : qid(1), total_run_count(0), m_stoprequested(false),
                m_ready(false)
                //m_worker_thread(boost::bind(&OpenCL_Handler::do_work, this))
{
    // populate the context list
    init_opencl();

    // allow the do_work thread to proceed
    {
        boost::mutex::scoped_lock lock(fifo_mutex);
        m_ready = true;
    }

    // create the worker thread now
    m_worker_thread = boost::thread(boost::bind(&OpenCL_Handler::do_work, this));
    //m_deleter_thread = boost::thread(boost::bind(&OpenCL_Handler::do_copy_and_delete, this));

}

OpenCL_Handler::~OpenCL_Handler()
{
    m_stoprequested = true;
    m_worker_thread.join(); // kill off the worker thread
    //m_deleter_thread.join();

    boost::mutex::scoped_lock fifo_lock(fifo_mutex);
    // now delete anything in the queues
    while (fifo_queued.size())
    {
        delete fifo_queued.front();
        fifo_queued.pop();
    }
    while (fifo_allocated.size())
    {
        delete fifo_allocated.front();
        fifo_allocated.pop();
    }
    while (fifo_running.size())
    {
        delete fifo_running.front();
        fifo_running.pop();
    }

//     while (fifo_completed.size())
//     {
//         delete fifo_completed.front();
//         fifo_completed.pop();
//     }


    // now delete the OpenCL resources
    // programs
    for (PMap::iterator it=prog_map.begin(), end=prog_map.end(); it != end; ++it)
    {
        clReleaseProgram(it->second);
    }
    prog_map.clear();

    // command queues
    for (std::vector<cl_command_queue>::iterator it=ocl_queues.begin(), end=ocl_queues.end(); it != end; ++it)
    {
        clReleaseCommandQueue(*it);
    }

    // now delete allocated resources
    delete[] ocl_platform_ids;

    if(ocl_gpu_context) clReleaseContext(ocl_gpu_context);
}

void OpenCL_Handler::context_error(const char *errinfo,
                            const void *private_info,
                            size_t cb,
                            void *user_data)
{
    std::cerr << "OpenCL Context Error: " << errinfo << std::endl;
}

long OpenCL_Handler::add_work_to_queue(OpenCL_WorkBlob* res_ptr, bool tracked)
{
    // get a lock on the fifo queue
    boost::mutex::scoped_lock l(fifo_mutex);

    long tracking_qid = (tracked) ? ++qid : 0;
    res_ptr->qid = tracking_qid;

    //std::cout << "Adding " << res_ptr << " (" << tracking_qid << ")" << std::endl;
    
    fifo_queued.push(res_ptr);

    return tracking_qid;
}

long OpenCL_Handler::add_work_to_queue_blocking(OpenCL_WorkBlob* res_ptr, bool tracked)
{
    assert(res_ptr->auto_delete == true);
    boost::mutex blocking_mutex;
    res_ptr->set_lock(blocking_mutex);

    long tracking_qid = (tracked) ? ++qid : 0;
    res_ptr->qid = tracking_qid;

    // now get a lock on the fifo queue, and submit the work
    {
        boost::mutex::scoped_lock l(fifo_mutex);
        fifo_queued.push(res_ptr);
    }

    // now try to get a lock on the blocking_mutex- this
    // will not be possible until the res_ptr has been deleted
    boost::mutex::scoped_try_lock tlock(blocking_mutex);
    //std::cout << "Try lock loop" << std::endl;
    while (tlock.try_lock() == false)
    {
        boost::this_thread::sleep(boost::posix_time::microseconds(10.0));
    }
    //std::cout << "Try lock loop succeeded" << std::endl;

    return tracking_qid;
}

void OpenCL_Handler::wait_until_idle()
{
    while(pending() > 0)
    {
        // nothing to do- suspend the thread for a little while
        boost::this_thread::sleep(boost::posix_time::microseconds(10.0));
    }
}

// get a list of all OpenCL devices on the system and create a context handler for each
void OpenCL_Handler::init_opencl()
{
    //std::cerr << "Init'ing OpenCL\n";
    
    // Get OpenCL platform count - allow up to 1024 to be returned as a 
    // first guess (the spec should allow us to just pass a NULL pointer, 
    // but that causes a memory leak within the NVIDIA OpenCL library,
    // hence this unsightly hack.
    cl_uint num_platforms;
    cl_platform_id* dummy = new cl_platform_id[1024];
    clGetPlatformIDs (1024, dummy, &num_platforms);
    delete[] dummy;

    // Get all OpenCL platform IDs
    ocl_platform_ids = new cl_platform_id[num_platforms];
    clGetPlatformIDs (num_platforms, ocl_platform_ids, NULL);

    // Select NVIDIA platform (this example assumes it IS present)
    // TODO:: what to do if more than one or none at all??
    char cBuffer[1024];
    std::vector<cl_uint> platform_ids;
    for(cl_uint i = 0; i < num_platforms; ++i)
    {
        clGetPlatformInfo (ocl_platform_ids[i], CL_PLATFORM_NAME, 1024, cBuffer, NULL);
        //std::cerr << "OpenCL platform info: " << cBuffer << std::endl;
        if(strstr(cBuffer, "NVIDIA") != NULL || 
           strstr(cBuffer, "ATI Stream") != NULL || 
           strstr(cBuffer, "AMD Accelerated Parallel Processing") != NULL)
        {
            platform_ids.push_back(i);
        }
    }

    // sanity check on nvidia platforms
    if (platform_ids.size() != 1) {
        std::cerr << "Error init_opencl(), Line " << static_cast<int>(__LINE__) << " in file " << __FILE__ << " !!!\n";
        std::cerr << "Got " << static_cast<int>(platform_ids.size()) << " OpenCL-compatible (NVIDIA / ATI Stream) platforms.\n";
        throw std::exception();
    }

    cl_uint cl_num_devices;

    // clear the device list array
    for (int ii=0; ii < MAX_DEVICES; ++ii) {
        ocl_devices[ii] = NULL;
    }

    //Get a GPU device on Platform (this example assumes one IS present)
    clGetDeviceIDs(ocl_platform_ids[platform_ids[0]], CL_DEVICE_TYPE_GPU, MAX_DEVICES, ocl_devices, &cl_num_devices);
    num_devices = static_cast<unsigned int>(cl_num_devices);
    //std::cout << "OpenCL found " << num_devices << " eligible devices" << std::endl;
    
    if (num_devices == 0) {
        std::cerr << "Error init_opencl(), Line " << static_cast<int>(__LINE__) << " in file " << __FILE__ << " !!!\n";
        std::cerr << "Failed to create any OpenCL contexts.  No devices??\n";
        throw std::exception();
    }
#ifdef SINGLE_GPU_ONLY
    if (num_devices > 1) {
        num_devices = 1;
    }
#endif
    // get the device constraints
    cl_ulong max_size;
    for (unsigned int ii=0; ii < cl_num_devices; ++ii)
    {
        cl_int err = clGetDeviceInfo(ocl_devices[ii], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), (void *)&max_size, NULL);
        assert(err == CL_SUCCESS);
        //std::cerr << "MAX_MEM_ALLOC_SIZE on device " << ii << " is " << max_size << std::endl;
        if (ii==0 || static_cast<size_t>(max_size) < max_ocl_malloc_size) {
            max_ocl_malloc_size = static_cast<size_t>(max_size);
        }
    }

    init_context_and_queues();
}

void OpenCL_Handler::init_context_and_queues()
{
    //Create a context
    int err;
    ocl_gpu_context = clCreateContext(0, num_devices, ocl_devices, &OpenCL_Handler::context_error, NULL, &err);
    if (err != CL_SUCCESS) { std::cout << "Error in init_context_and_queues: " << err << std::endl; }
    assert(err == CL_SUCCESS);

    // init OpenCL queue
    init_queues();

}

void OpenCL_Handler::init_queues()
{
    // since this can modify the context, need to ensure only one
    // thread at a time is using this function
    boost::mutex::scoped_lock state_lock(altering_context_mutex);
    for (unsigned int ctr=0; ctr < num_devices; ++ctr)
    {
        int err=0;

        // Create a command-queue
        //std::cerr << "clCreateCommandQueue...";
        cl_command_queue ocl_q = clCreateCommandQueue(ocl_gpu_context, ocl_devices[ctr], 0, &err);
        if (err != CL_SUCCESS)
        {
            std::cerr << "Error " << err << " in clCreateCommandQueue, Line " << __LINE__ << " in file " << __FILE__ << " !!!\n";
            throw std::exception();
        }
        ocl_queues.push_back(ocl_q);
        //std::cerr << "done\n";
    }

    // init the queue iterator (round-robins the work if more than one device found)
    ocl_qit = ocl_queues.begin();
    assert(ocl_qit != ocl_queues.end());

    return;
}

cl_program OpenCL_Handler::get_program(const unsigned char* source)
{
    PMap::const_iterator pit = prog_map.find(source);
    if (pit == prog_map.end())
    {
        // since this will modify the object, need to get a mutex lock
        boost::mutex::scoped_lock l(altering_context_mutex);

        std::string opencl_source((char*)source);

        // it is technically possible that another thread will have populated
        // the program map with the program in time it took us to get
        // the mutex lock.  So check again
        pit = prog_map.find(source);
        if (pit == prog_map.end())
        {
            cl_program new_program = compile_program(opencl_source, ocl_gpu_context, num_devices, ocl_devices);
            prog_map[source] = new_program;
            return new_program;
        }
    }

    return pit->second;
}

std::string OpenCL_Handler::load_program(const std::string& source_filename)
{

    // Read the OpenCL kernel in from source file
    //std::cerr << "Loading Source from (" << source_filename.c_str() << ")...";
    std::stringstream buffer;
    std::ifstream file( source_filename.c_str() );
    if ( file )
    {
        buffer << file.rdbuf();
        file.close();
    }
    else
    {
        std::cerr << "Unable to read " << source_filename << std::endl;
        throw std::exception();
    }
    //std::cerr << "done\n";

    return buffer.str();
}

cl_program OpenCL_Handler::compile_program(std::string prog_source, cl_context& context, unsigned int num_devices, cl_device_id device_ids[])
{
    // going to create this and then return it by value (I suspect the cl_program
    // is implemented as a pointer to something within the opencl context anyway, 
    // so this isn't too wasteful. Besides this function shouldn't be invoked very
    // often)
    
    const char* prog_sources_ptrs[1] = {(const char*)prog_source.c_str()}; // ouch
    const size_t size_array[1] = {prog_source.size()};
    //std::cerr << "Program source size: " << size_array[0] << "\n";

    // Create the program
    //std::cerr << "clCreateProgramWithSource...";
    int err;
    cl_program ocl_program = clCreateProgramWithSource(context, 1, prog_sources_ptrs, size_array, &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error " << err << " in clCreateProgramWithSource, Line " << __LINE__ << " in file " << __FILE__ << " !!!\n";
        throw std::exception();
    }
    //std::cerr << "done\n";

    // Build the program with 'mad' Optimization option
    #ifdef MAC
        std::string flags = "-cl-fast-relaxed-math -DMAC -cl-single-precision-constant";
    #else
        std::string flags = "-cl-fast-relaxed-math -cl-single-precision-constant";
    #endif
    //std::cerr << "clBuildProgram...";
    err = clBuildProgram(ocl_program, num_devices, device_ids, flags.c_str(), NULL, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error " << err << " in clBuildProgram, Line " << __LINE__ << " in file " << __FILE__ << " !!!\n";
        
        char build_output[2048];
        size_t size_out;
        clGetProgramBuildInfo(ocl_program, device_ids[0], CL_PROGRAM_BUILD_LOG, 2048, build_output, &size_out);
        std::cerr << std::string(build_output) << std::endl;

        throw std::exception();
    }
    //std::cerr << "done\n";

    return ocl_program;
}

void OpenCL_Handler::refresh_context()
{
    std::cout << "Refreshing context..." << std::endl;
    // programs
    for (PMap::iterator it=prog_map.begin(), end=prog_map.end(); it != end; ++it)
    {
        clReleaseProgram(it->second);
    }
    prog_map.clear();

    // command queues
    for (std::vector<cl_command_queue>::iterator it=ocl_queues.begin(), end=ocl_queues.end(); it != end; ++it)
    {
        clReleaseCommandQueue(*it);
    }
    ocl_queues.clear();

    // context
    if (ocl_gpu_context) clReleaseContext(ocl_gpu_context);
    ocl_gpu_context = NULL;

    // Reinit context and queues
    init_context_and_queues();
}

// the function in which the thread spins around
void OpenCL_Handler::do_work()
{
    long memory_count = 0l;
    size_t total_memory_count = 0l;
    const size_t MAX_DEVICE_MEMORY = num_devices*512*1024*1024l;
    
    // sleep until the object has initialised the queue, otherwise we'll get into a pickle
    while(true)
    {
        boost::this_thread::sleep(boost::posix_time::milliseconds(10.0));
        {
            boost::mutex::scoped_lock l(fifo_mutex);
            if (m_ready || m_stoprequested) { break; }
        }
    }

    try
    {

        while (!m_stoprequested)
        {
        
            // wrap access to the fifo_queued with some mutex locking
            OpenCL_WorkBlob* front_runner=NULL;
            size_t queued=0, allocated=0, running=0;
            {
                boost::mutex::scoped_lock l(fifo_mutex);
                queued = fifo_queued.size();
                allocated = fifo_allocated.size();
                running = fifo_running.size();

                // something is running; check to see if it's finished
                if (running > 0) {
                    front_runner = fifo_running.front();
                }
            }

            if (queued == 0 && running == 0 && allocated == 0)
            {
                // nothing to do- suspend the thread for a little while
                boost::this_thread::sleep(boost::posix_time::milliseconds(10.0));
            }
            else
            {
                //std::cout << "Start- running: " << running << " allocated: " << allocated << " queued: " << queued << std::endl;
                bool full_memory = false; // keep track of how much memory has been allocated 

                try {
                    while (full_memory == false && queued > 0) 
                    {
                        OpenCL_WorkBlob* res_ptr;
                        {
                            boost::mutex::scoped_lock l(fifo_mutex);
                            res_ptr = fifo_queued.front();
                        }

                        size_t mem_req = res_ptr->device_mem_reqd();

                        if (mem_req >= MAX_DEVICE_MEMORY) {
                            std::cout << "Item is too large to allocate memory: " << mem_req << " bytes" << std::endl;
                            throw std::exception();
                        }


                        if ((memory_count+mem_req) >= MAX_DEVICE_MEMORY) {
                            full_memory = true;
                            //std::cout << "Full memory: " << memory_count << " mem_req: " << mem_req << " (total so far: " << total_memory_count << " bytes) items in queues:" << queued << "/" << allocated << "/" << running << std::endl;
                            continue;
                        }

                        // The Resource block need's init'ing prior to execution on the GPU
                        // This is where the memory and stuff gets allocated/copied.
                        // Nb: can throw exceptions, so don't alter the fifos/counters yet...
                        const unsigned char* id = res_ptr->get_opencl_source();
                        if (id)
                        {
                            cl_program program = get_program(id);
                            res_ptr->allocate(ocl_gpu_context, program);
                        }

                        memory_count += mem_req;
                        total_memory_count += mem_req;

                        // ok now that the thing has allocated memory successfully, can
                        // change the fifo's / counters
                        {
                            boost::mutex::scoped_lock l(fifo_mutex);
                            fifo_allocated.push(res_ptr);
                            fifo_queued.pop();
                        }

                        ++allocated;
                        --queued;
                    }
                }
                catch (OpenCL_WorkBlob::AllocationException) {
                    boost::this_thread::sleep(boost::posix_time::milliseconds(1.0));
                }

                //std::cout << "After alloc loop- running: " << running << " allocated: " << allocated << " queued: " << queued << std::endl;
                try {
                    // otherwise there's something useful we can be doing...
                    while (allocated > 0)
                    {
                        // pop something off the 'allocated' fifo_queue
                        // add it to the running queue
                        OpenCL_WorkBlob* res_ptr;
                        {
                            boost::mutex::scoped_lock l(fifo_mutex);
                            assert(fifo_allocated.size() == allocated);
                            res_ptr = fifo_allocated.front();
                        }
                        
                        // NB: enqueue might throw an exception, so don't change the fifo's yet...
                        res_ptr->enqueue(*ocl_qit);
                        clFlush(*ocl_qit);
                        if (++ocl_qit == ocl_queues.end()) { ocl_qit = ocl_queues.begin(); }

                        // ok it's in the queue successfully, so update the fifos and counters
                        {
                            boost::mutex::scoped_lock l(fifo_mutex);
                            fifo_running.push(res_ptr);
                            fifo_allocated.pop();
                        }
                        ++running;
                        --allocated;

                        ++total_run_count;

                        // if the one we just added was the first in the queue, set
                        // the front runner pointer
                        if (running == 1) { front_runner = res_ptr; }
                    }
                }
                catch (OpenCL_WorkBlob::EnqueueException) {
                    boost::this_thread::sleep(boost::posix_time::milliseconds(1.0));
                }

                //std::cout << "After run loop- running: " << running << " allocated: " << allocated << " queued: " << queued << std::endl;
                try {
                    while (running > 0 && front_runner->finished())
                    {
                        //boost::mutex::scoped_lock fifo_lock(fifo_mutex);
                        //fifo_completed.push(res_ptr);

                        front_runner->copy_out_results();
                        long tracking_qid = front_runner->qid; // this is a uid tracking items in the queue

                        // decrement the amount of memory
                        memory_count -= front_runner->device_mem_reqd();
                        assert(memory_count >= 0);

                        if (front_runner->auto_delete) {
                            delete front_runner;
                        }

                        boost::mutex::scoped_lock fifo_lock(fifo_mutex);
                        fifo_running.pop();
                        --running;
                        
                        // the 'zero' tracking_qid means, don't track the
                        // item, so don't bother entering in completion list
                        if (tracking_qid != 0) {
                            completed_tracked.push_back(tracking_qid);
                        }
                        
                        if (running > 0) {
                            front_runner = fifo_running.front();
                        }
                        else
                        {
                            front_runner = NULL;
                        }
                    }
                
                }
                catch (...) 
                {
                    std::cout << "Memory: " << memory_count << " (total so far: " << total_memory_count << " bytes) items in queues:" << queued << "/" << allocated << "/" << running << std::endl;
                    assert(false);
                    throw;
                }
            }
        }
    }
    catch (std::exception)
    {
        std::cerr << "do_work caught an exception!" << std::endl;
        throw;
    }

    return;
}

// could run this code in a separate thread to handle the copy/deleting bit of do_work
// but it didn't really help much, so is commented out.  Haven't totally deleted it though
// just in case it's useful in future.
#if 0
void OpenCL_Handler::do_copy_and_delete()
{

    // sleep until the object has initialised the queue, otherwise we'll get into a pickle
    while(!m_ready && !m_stoprequested)
    {
        boost::this_thread::sleep(boost::posix_time::microseconds(10.0));
    }

    try
    {

        while (!m_stoprequested)
        {

            size_t completed=0;
            {
                boost::mutex::scoped_lock l(fifo_mutex);
                completed = fifo_completed.size();
            }

            if (completed == 0)
            {
                // nothing to do- suspend the thread for a little while
                boost::this_thread::sleep(boost::posix_time::microseconds(10.0));
            }
            else
            {
                // otherwise there's something useful we can be doing...
                while (completed > 0)
                {
                    // pop something off the 'allocated' fifo_queue
                    // add it to the running queue
                    OpenCL_WorkBlob* res_ptr;
                    {
                        boost::mutex::scoped_lock l(fifo_mutex);
                        res_ptr = fifo_completed.front();
                    }

                    res_ptr->copy_out_results();
                    long tracking_qid = res_ptr->qid; // this is a uid tracking items in the queue
                    
                    if (res_ptr->auto_delete) {
                        delete res_ptr;
                    }

                    {
                        boost::mutex::scoped_lock l(fifo_mutex);
                        fifo_completed.pop();
                        
                        // the 'zero' tracking_qid means, don't track the
                        // item, so don't bother entering in completion list
                        if (tracking_qid != 0) {
                            completed_tracked.push_back(tracking_qid);
                        }
                    }

                    --completed;
                }
            }
        }
    }
    catch (std::exception)
    {
        std::cerr << "do_copy_and_delete caught an exception!" << std::endl;
    }

    return;
}
#endif

#endif // ifdef OPENCL

