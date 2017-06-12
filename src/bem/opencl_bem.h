/*
* opencl_bem.h
*
*  Created on: 22 Sep 2010
*      Author: david
*/

#ifndef OPENCL_BEM_H_
#define OPENCL_BEM_H_

#ifdef OPENCL

#include <vector>
#include <string>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/thread/mutex.hpp>

// fwd decls.
class BasicNodePatch;
class LocalIntegrations;

#include <CL/opencl.h>
#include "opencl_bem_structs.h"
#include "../opencl/opencl_workblob.h"

typedef std::vector<const BasicNodePatch*> PPList;
typedef boost::shared_ptr< PPList > PatchPtrList;

// This is a holder which keeps QuadPoints in scope -- will
// get deleted by the OpenCL handler like a usual WorkBlob, 
// and when it does ... **ping** the references to the 
// many many many quad points in existence magically vanish 
// and they disappear rather than clogging up all the memory.
class QuadPointCache : public OpenCL_WorkBlob
{
public:
        
    QuadPointCache(std::vector< boost::shared_ptr<QuadList> >& cached)
    {
        stored.insert(stored.begin(), cached.begin(), cached.end());
    }

    virtual void allocate(const cl_context&, const cl_program&) {}
    virtual void copy_out_results() {}
    virtual void enqueue(cl_command_queue&) {}
    virtual bool finished() const { return true; }
    virtual size_t device_mem_reqd() const { return 0; }

    virtual ~QuadPointCache() {}

private:
    
    std::vector< boost::shared_ptr<QuadList> > stored;
        
};

class BEM_Resources : public OpenCL_WorkBlob
{

public:

    static const unsigned int xdim = OPENCL_LOCAL_GROUP_DIM_SOURCES;
    static const unsigned int ydim = OPENCL_LOCAL_GROUP_DIM_TARGETS;

    BEM_Resources(PatchPtrList incoming_srcs,
            PatchPtrList incoming_targs,
            double incoming_kappa,
            LocalIntegrations* incoming_final_destination);

    virtual ~BEM_Resources();
    
    // this function returns the source code (as array of chars, read in from .cl.h file, generated
    // by xxd (see above). NB: it is a pure virtual function in the base class
    const unsigned char* get_opencl_source();

    void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program);
    bool finished() const;
    inline void copy_out_results();
    void enqueue(cl_command_queue& ocl_queue);
    size_t device_mem_reqd() const;

private:

    void init_OpenCL_structs();
    void allocate_device_mem(const cl_context ocl_context);
    void init_kernel(const cl_context& ocl_context,
                    const cl_program& ocl_program);

private:

    // BEM parameters
    cl_float kappa;
    cl_uint num_sources;
    cl_uint num_targs;
    cl_uint num_quads;
    size_t num_results;

    // Where to put the final results
    LocalIntegrations* final_destination;
    size_t results_size;

    PatchPtrList stashed_source_patches;
    PatchPtrList stashed_target_patches;

    // host memory
    boost::shared_array<OpenCL_NodePatch> host_srcs;
    boost::shared_array<OpenCL_NodePatch> host_targs;
    boost::shared_array<OpenCL_QuadPoint> host_quads;
    //boost::shared_array<OpenCL_LocalIntegrationResult> host_results;

    // device memory
    cl_mem cmDevResults;
    cl_mem cmDevSrcs;
    cl_mem cmDevTargs;
    cl_mem cmDevQuads;

    // OpenCL control stuff
    cl_kernel m_kernel;
    
    // Events -- track kernels executed
    static const unsigned int NUM_EVENTS = 2;
    cl_event m_events[NUM_EVENTS];
    
    // Quadrature point cacheing -- this list of shared_ptrs
    // prevents deletion of the underlying QuadList until all
    // references in all BEM_Resource objects have finished
    // with them
    //std::vector<boost::shared_ptr<QuadList> > qp_cache;
};

class BEM_OnDemand_Resources : public OpenCL_WorkBlob
{

public:

    static const unsigned int xdim = OPENCL_LOCAL_GROUP_DIM_SOURCES;
    static const unsigned int ydim = OPENCL_LOCAL_GROUP_DIM_TARGETS;

    BEM_OnDemand_Resources(PatchPtrList incoming_srcs,
                            PatchPtrList incoming_targs,
                            double incoming_kappa,
                            const double* incoming_f,
                            const double* incoming_h,
                            double* output_f,
                            double* output_h,
                            bool auto_del=true,
                            boost::mutex *results_mutex=NULL);

    virtual ~BEM_OnDemand_Resources();
    
    const unsigned char* get_opencl_source();

    void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program);
    bool finished() const;
    void copy_out_results();
    void enqueue(cl_command_queue& ocl_queue);
    size_t device_mem_reqd() const;

private:

    void init_OpenCL_structs();
    void allocate_device_mem(const cl_context ocl_context);
    void init_kernel(const cl_context& ocl_context,
                    const cl_program& ocl_program);

    // BEM parameters
    cl_float kappa;
    cl_uint num_sources;
    cl_uint num_targs;
    cl_uint num_quads;
    size_t block_results_size;

    // Where to put the final results
    double* final_destinations_f;
    double* final_destinations_h;
    boost::mutex *results_mutex_ptr; // (optional) a mutex to lock before copying the results

    // arrays of src/target patches
    PatchPtrList stashed_source_patches;
    PatchPtrList stashed_target_patches;

    const double* stashed_incoming_f;
    const double* stashed_incoming_h;

    // host memory
    boost::shared_array<OpenCL_NodePatch> host_srcs;
    boost::shared_array<OpenCL_NodePatch> host_targs;
    boost::shared_array<OpenCL_QuadPoint> host_quads;
    boost::shared_array<cl_float2> host_results;
    boost::shared_array<cl_float> host_fh;

    // device memory
    cl_mem cmDevBlockResults;
    cl_mem cmDevResults;
    cl_mem cmDevFH;
    cl_mem cmDevSrcs;
    cl_mem cmDevTargs;
    cl_mem cmDevQuads;

    // OpenCL control stuff
    cl_kernel m_kernel;
    cl_kernel m_blk_results_kernel;

    // Events -- track kernels executed
    static const unsigned int NUM_EVENTS = 3;
    cl_event m_events[NUM_EVENTS];

    size_t szGlobalWorkSize[2];
    size_t szLocalWorkSize[2];

    cl_uint blk_xdim;
    cl_uint blk_ydim;

    // Quadrature point cacheing -- this list of shared_ptrs
    // prevents deletion of the underlying QuadList until all
    // references in all BEM_Resource objects have finished
    // with them
    //std::vector<boost::shared_ptr<QuadList> > qp_cache;
    
};

class SingularBEM : public OpenCL_WorkBlob
{

public:

    static const unsigned int xdim = 16;
    
    SingularBEM(PatchPtrList incoming_srcs,
                double incoming_kappa,
                LocalIntegrations* incoming_final_destination);

    virtual ~SingularBEM();
    
    const unsigned char* get_opencl_source();

    void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program);
    bool finished() const;
    inline void copy_out_results();
    void enqueue(cl_command_queue& ocl_queue);
    size_t device_mem_reqd() const;

private:

    void init_OpenCL_structs();
    void allocate_device_mem(const cl_context ocl_context);
    void init_kernel(const cl_context& ocl_context,
                    const cl_program& ocl_program);

private:

    // BEM parameters
    cl_float kappa;
    cl_uint num_patches;
    cl_uint num_quads;
    size_t num_results;

    // Where to put the final results
    LocalIntegrations* final_destination;
    size_t results_size;

    PatchPtrList stashed_source_patches;
    
    // host memory
    boost::shared_array<OpenCL_NodePatch> host_srcs;
    boost::shared_array<OpenCL_QuadPoint> host_quads;
    //boost::shared_array<OpenCL_LocalIntegrationResult> host_results;

    // device memory
    cl_mem cmDevResults;
    cl_mem cmDevSrcs;
    cl_mem cmDevQuads;

    // OpenCL control stuff
    cl_kernel m_kernel;
    
    // Events -- track kernels executed
    static const unsigned int NUM_EVENTS = 2;
    cl_event m_events[NUM_EVENTS];

    
    // Quadrature point cacheing -- this list of shared_ptrs
    // prevents deletion of the underlying QuadList until all
    // references in all BEM_Resource objects have finished
    // with them
    std::vector<boost::shared_ptr<QuadList> > qp_cache;
    
};

#endif // ifdef OPENCL

#endif /* OPENCL_RESOURCES_H_ */

