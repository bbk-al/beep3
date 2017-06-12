/*
* opencl_resources.h
*
*  Created on: 22 Sep 2010
*      Author: david
*/

#ifndef OPENCL_FMM_H_
#define OPENCL_FMM_H_
#ifdef OPENCL

// these will get pulled in via the opencl_fmm.cpp which #includes the openecl_fmm_kernels.cl.h
extern unsigned char opencl_fmm_kernels_cl[];
extern unsigned int opencl_fmm_kernels_cl_len;

// these control the blocksizes on the GPU
// VERY IMPORTANT: don't change these without also changing them in opencl_fmm_kernels.cl
// and then updating opencl_fmm_kernels.cl.h by running "xxd -i" on the cl file to create
// the header file.
#define OPENCL_LOCAL_GROUP_DIM_EVAL_PTS 8
#define OPENCL_LOCAL_GROUP_DIM_CHARGES  8
#define BLOCK_RESULTS_COLLECTION_WIDTH (OPENCL_LOCAL_GROUP_DIM_EVAL_PTS*OPENCL_LOCAL_GROUP_DIM_CHARGES)

#include <vector>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/thread/mutex.hpp>

#include "eval_pt.h"
#include "../common/charge.h"

#include <CL/opencl.h>
#include "../opencl/opencl_handler.h"
#include "../opencl/opencl_workblob.h"

// Probable limit on maximum threads per block-
// could query the GPU to get this directly, but
// probably not a massively critical param here
#define MAX_BLOCKS_PER_GRID 64*1024

using fmm::EvalPt;
using fmm::EvalPt_2ndDerivs;

template <typename EvalPtType>
class FMM_Derivatives_Level
{};

template <>
class FMM_Derivatives_Level<EvalPt_2ndDerivs>
{
public:
    typedef cl_float16 CL_VECTORIZED_FLOAT_TYPE;
    static const std::string fmm_kernel_name;
    static const std::string block_results_kernel_name;
};

template <>
class FMM_Derivatives_Level<EvalPt>
{
public:
    typedef cl_float4 CL_VECTORIZED_FLOAT_TYPE;
    static const std::string fmm_kernel_name;
    static const std::string block_results_kernel_name;
};

template<typename EvalPtType>
class FMM_Resources : public OpenCL_WorkBlob
{

public:

    static const unsigned int xdim = OPENCL_LOCAL_GROUP_DIM_EVAL_PTS;
    static const unsigned int ydim = OPENCL_LOCAL_GROUP_DIM_CHARGES;
    
    // this defines the appropriate float type
    typedef typename FMM_Derivatives_Level<EvalPtType>::CL_VECTORIZED_FLOAT_TYPE CL_VECTORIZED_FLOAT_TYPE;
    
    FMM_Resources(const std::vector<EvalPtType*>& eval_pts,
                    const std::vector<Charge*>& charges,
                    double kappa_scaled);
    virtual ~FMM_Resources();

    // this function returns the source code (as array of chars, read in from .cl.h file, generated
    // by xxd. NB: it is a pure virtual function in the base class
    const unsigned char* get_opencl_source();
    
    void allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program);
    bool finished() const;
    void copy_out_results();
    void enqueue(cl_command_queue& ocl_queue);

    // this calculates the number of chunks that a lump of FMM work needs to be split into
    // in order to fit within the limits of the GPU
    static unsigned int calc_x_chunksize(size_t num_eval_pts, size_t num_charges, size_t max_malloc_size)
    {
    
        size_t block_results_width = sizeof(CL_VECTORIZED_FLOAT_TYPE)*BLOCK_RESULTS_COLLECTION_WIDTH*(static_cast<unsigned int>(ceil(static_cast<double>(num_eval_pts)/(BLOCK_RESULTS_COLLECTION_WIDTH))));
        
        // arbitrary cutoff here- if the block width is greater than 1MB per row, then split it
        unsigned int x_chunks = static_cast<unsigned int>(ceil(static_cast<double>(block_results_width) / (1024*1024)));
        
        unsigned int blocks_per_grid_limit = static_cast<unsigned int>(ceil(static_cast<double>(num_eval_pts) / (xdim*MAX_BLOCKS_PER_GRID)));
        x_chunks = (x_chunks > blocks_per_grid_limit) ? x_chunks : blocks_per_grid_limit;

        return x_chunks;
    }
    
    static unsigned int calc_y_chunksize(size_t num_eval_pts, size_t num_charges, size_t max_malloc_size)
    {
        unsigned int charges_limit = static_cast<unsigned int>(ceil(static_cast<double>(num_charges * sizeof(cl_float4)) / (16*1024*1024)));
        unsigned int blocks_per_grid_limit = static_cast<unsigned int>(ceil(static_cast<double>(num_charges)  / (ydim*100*MAX_BLOCKS_PER_GRID))); // the opencl code also loops over the y direction


        size_t block_results_width = sizeof(CL_VECTORIZED_FLOAT_TYPE)*BLOCK_RESULTS_COLLECTION_WIDTH*(static_cast<unsigned int>(ceil(static_cast<double>(num_eval_pts)/(BLOCK_RESULTS_COLLECTION_WIDTH))));
        unsigned int x_chunks = static_cast<unsigned int>(ceil(static_cast<double>(block_results_width) / (1024*1024)));
        block_results_width /= x_chunks;
        unsigned int block_results_limit = (16*1024*1024) / static_cast<unsigned int>(block_results_width);
        unsigned int block_chunks = static_cast<unsigned int>(ceil(static_cast<double>(num_charges) / (ydim*100*block_results_limit))); // the opencl code also loops over the y direction

        unsigned int y_chunks = (charges_limit > blocks_per_grid_limit) ? charges_limit : blocks_per_grid_limit;
        y_chunks = (y_chunks > block_chunks) ? y_chunks : block_chunks;
        return y_chunks;
    }

    size_t device_mem_reqd() const;

private:

    void allocate_device_mem(const cl_context& ocl_context);
    void init_kernel(const cl_context& ocl_context,
                    const cl_program& ocl_program);

    cl_float kappa;
    cl_uint num_eval_pts;
    cl_uint num_eval_pts_rounded;
    cl_uint num_charges;
    cl_uint num_charges_rounded;
    cl_uint num_loops;
    size_t block_results_size;
    size_t results_size;

    // host memory
    boost::shared_array<cl_float4> host_charges;
    boost::shared_array<cl_float4> host_pts;
    boost::shared_array<CL_VECTORIZED_FLOAT_TYPE> host_results;
    std::vector<EvalPtType*> eval_pts_copy;

    // device memory
    cl_mem cmDevCharges;
    cl_mem cmDevEvalPts;
    cl_mem cmDevBlockResults;
    cl_mem cmDevFinalResults;

    // OpenCL control stuff
    cl_kernel m_evaluation_kernel;
    cl_kernel m_blk_results_kernel;
    cl_event m_event;

    size_t szLocalWorkSize[2];
    size_t szGlobalWorkSize[2];

};

template <typename EvalPtType> 
FMM_Resources<EvalPtType>::FMM_Resources(const std::vector<EvalPtType*>& eval_pts,
                            const std::vector<Charge*>& charges,
                            double kappa_scaled) :
                            kappa(static_cast<cl_float>(kappa_scaled)),
                            cmDevCharges(NULL),
                            cmDevEvalPts(NULL),
                            cmDevFinalResults(NULL),
                            cmDevBlockResults(NULL),
                            m_evaluation_kernel(NULL),
                            m_blk_results_kernel(NULL),
                            m_event(NULL)
{
    OpenCL_WorkBlob::auto_delete = true;

    szLocalWorkSize[0] = static_cast<size_t>(xdim);
    szLocalWorkSize[1] = static_cast<size_t>(ydim);
 
    num_eval_pts = static_cast<cl_uint>(eval_pts.size());
    num_charges = static_cast<cl_uint>(charges.size());

    // sanity check
    assert(num_eval_pts > 0);
    assert(num_charges > 0);
    
    // assuming that we want a limited number of rows of thread blocks
    // we might need to loop over charges in the kernels rather than one
    // per thread.
    //const unsigned int num_rows_max = 1 + (16*1024*1024 / (num_eval_pts * 16));

    num_eval_pts_rounded = (BLOCK_RESULTS_COLLECTION_WIDTH) * static_cast<size_t>(ceil(static_cast<double>(num_eval_pts) / (BLOCK_RESULTS_COLLECTION_WIDTH)));
    num_charges_rounded  = szLocalWorkSize[1] * static_cast<size_t>(ceil(static_cast<double>(num_charges) / szLocalWorkSize[1]));

    size_t block_results_width = sizeof(CL_VECTORIZED_FLOAT_TYPE)*xdim*(static_cast<unsigned int>(ceil(static_cast<double>(num_eval_pts)/xdim)));
    unsigned int block_results_limit = (16*1024*1024) / static_cast<unsigned int>(block_results_width);

    size_t num_rows = num_charges_rounded / szLocalWorkSize[1];
    num_loops = static_cast<cl_uint>(ceil(static_cast<double>(num_rows) / block_results_limit));
    assert(num_loops <= 100); // otherwise the caller should have pre-split the work into separate units
    
    // (actually I think this is redundant - if too many blocks per grid we will have split it up already...)
    unsigned int group_size_limit = MAX_BLOCKS_PER_GRID*ydim;
    cl_uint group_limited_num_loops = static_cast<cl_uint>(ceil(static_cast<double>(num_charges_rounded) / group_size_limit));
    num_loops = (num_loops > group_limited_num_loops) ? num_loops : group_limited_num_loops;
    
    num_charges_rounded =  num_loops * szLocalWorkSize[1] * static_cast<size_t>(ceil(static_cast<double>(num_charges) / (num_loops*szLocalWorkSize[1])));

    host_pts = boost::shared_array<cl_float4>(new cl_float4[num_eval_pts_rounded]);
    host_charges = boost::shared_array<cl_float4>(new cl_float4[num_charges_rounded]);
    host_results = boost::shared_array<CL_VECTORIZED_FLOAT_TYPE>(new CL_VECTORIZED_FLOAT_TYPE[num_eval_pts_rounded]);
    assert(host_pts.get()); // check the memory allocation worked
    assert(host_charges.get()); // check the memory allocation worked
    assert(host_results.get()); // check the memory allocation worked

    // copy the evaluation points
    eval_pts_copy.reserve(num_eval_pts);
    eval_pts_copy.insert(eval_pts_copy.begin(), eval_pts.begin(), eval_pts.end());
    for (unsigned int ii=0; ii < num_eval_pts; ++ii)
    {
        EvalPtType& ep = *(eval_pts[ii]);
        cl_float4& ocl_ep = host_pts[ii];
        ocl_ep.x = static_cast<float>(ep.x);
        ocl_ep.y = static_cast<float>(ep.y);
        ocl_ep.z = static_cast<float>(ep.z);
    }
    for (unsigned int blanks=num_eval_pts; blanks < num_eval_pts_rounded; ++blanks)
    {
        cl_float4& ocl_ep = host_pts[blanks];
        ocl_ep.x = 0;
        ocl_ep.y = 0;
        ocl_ep.z = 0;
    }

    // copy the charges
    for (unsigned int ii=0; ii < num_charges; ++ii)
    {
        cl_float4& ocl_ch = host_charges[ii];
        const Charge& fmm_ch = *(charges[ii]);
        ocl_ch.w = static_cast<float>(fmm_ch.get_charge());
        ocl_ch.x = static_cast<float>(fmm_ch.x);
        ocl_ch.y = static_cast<float>(fmm_ch.y);
        ocl_ch.z = static_cast<float>(fmm_ch.z);
    }
    for (unsigned int blanks=num_charges; blanks < num_charges_rounded; ++blanks)
    {
        cl_float4& ocl_ch = host_charges[blanks];
        ocl_ch.w = 0;        
        ocl_ch.x = 0;
        ocl_ch.y = 0;
        ocl_ch.z = 0;
    }

    szGlobalWorkSize[0] = num_eval_pts_rounded;
    szGlobalWorkSize[1] = num_charges_rounded / num_loops;
    //szGlobalWorkSize[0] = szLocalWorkSize[0] * static_cast<size_t>(ceil(static_cast<double>(num_eval_pts) / szLocalWorkSize[0]));
    //szGlobalWorkSize[1] = szLocalWorkSize[1] * static_cast<size_t>(ceil(static_cast<double>(num_charges) / szLocalWorkSize[1]));

    assert(szGlobalWorkSize[0] > 0);
    assert(szGlobalWorkSize[1] > 0);
    assert(szGlobalWorkSize[0]/szLocalWorkSize[0] <= MAX_BLOCKS_PER_GRID);
    assert(szGlobalWorkSize[1]/szLocalWorkSize[1] <= MAX_BLOCKS_PER_GRID);

    block_results_size = sizeof(CL_VECTORIZED_FLOAT_TYPE)*szGlobalWorkSize[0]*szGlobalWorkSize[1]/szLocalWorkSize[1];
    results_size = sizeof(CL_VECTORIZED_FLOAT_TYPE)*num_eval_pts_rounded;
//     std::cout << "Global worksizes [x/y] = " << szGlobalWorkSize[0] << "," << szGlobalWorkSize[1] << "\n";
//     std::cout << "Local worksizes [x/y] = " << szLocalWorkSize[0] << "," << szLocalWorkSize[1] << "\n";
//     std::cout << "block results size: " << block_results_size << "\n";
//     std::cout << "results size: "<< results_size << "\n";
//     std::cout << "eval_pts and num_charges rounded: " << num_eval_pts_rounded << ", " << num_charges_rounded << "\n";
//     std::cout << "num loops: " << num_loops << "\n";
    size_t mem_req = device_mem_reqd();
    if (mem_req >= 512*1024*1024) {
        std::cout << "Massive FMM Chunk: " << mem_req << std::endl;
        throw std::exception();
    }

}

template <typename EvalPtType> 
FMM_Resources<EvalPtType>::~FMM_Resources() {

    host_pts.reset();
    host_charges.reset();
    host_results.reset();

    if(m_blk_results_kernel != NULL) { int err = clReleaseKernel(m_blk_results_kernel); assert(err == CL_SUCCESS); }
    if(m_event != NULL) { int err = clReleaseEvent(m_event); assert(err == CL_SUCCESS); }
    if(m_evaluation_kernel != NULL) { int err = clReleaseKernel(m_evaluation_kernel); assert(err == CL_SUCCESS); }

    if(cmDevCharges != NULL) { int err = clReleaseMemObject(cmDevCharges); assert(err == CL_SUCCESS); }
    if(cmDevEvalPts != NULL) { int err = clReleaseMemObject(cmDevEvalPts); assert(err == CL_SUCCESS); }
    if(cmDevBlockResults != NULL) { int err = clReleaseMemObject(cmDevBlockResults); assert(err == CL_SUCCESS); }
    if(cmDevFinalResults != NULL) { int err = clReleaseMemObject(cmDevFinalResults); assert(err == CL_SUCCESS); }

    //std::cout << "Destroyed " << this << "(" << qid << ")" << std::endl;
}

template <typename EvalPtType> 
void FMM_Resources<EvalPtType>::allocate(const cl_context& ocl_gpu_context, const cl_program& ocl_program)
{
    try {
        allocate_device_mem(ocl_gpu_context);
        init_kernel(ocl_gpu_context, ocl_program);
    }
    catch (...)
    {
        // cleanup any memory allocated so far
        if(m_blk_results_kernel != NULL) { int err = clReleaseKernel(m_blk_results_kernel); assert(err == CL_SUCCESS); }
        if(m_event != NULL) { int err = clReleaseEvent(m_event); assert(err == CL_SUCCESS); }
        if(m_evaluation_kernel != NULL) { int err = clReleaseKernel(m_evaluation_kernel); assert(err == CL_SUCCESS); }
        
        if(cmDevCharges != NULL) { int err = clReleaseMemObject(cmDevCharges); assert(err == CL_SUCCESS); }
        if(cmDevEvalPts != NULL) { int err = clReleaseMemObject(cmDevEvalPts); assert(err == CL_SUCCESS); }
        if(cmDevBlockResults != NULL) { int err = clReleaseMemObject(cmDevBlockResults); assert(err == CL_SUCCESS); }
        if(cmDevFinalResults != NULL) { int err = clReleaseMemObject(cmDevFinalResults); assert(err == CL_SUCCESS); }
        // rethrow- handler should stop trying to enqueue things for a while!
        // (should retry to allocate this object in next iteration of do_work loop.)
       
        throw;
    }
}

template <typename EvalPtType> 
bool FMM_Resources<EvalPtType>::finished() const
{
    // find out if the kernel event has completed
    cl_int retval;
    cl_int err = clGetEventInfo(m_event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), static_cast<void*>(&retval), NULL);
    if (err != CL_SUCCESS) { 
        std::cerr << "template <typename EvalPtType> FMM_Resources<EvalPtType>::finished(): clGetEventInfo return error code: " << err << std::endl; 
        throw std::exception();
    }
    if (retval < 0)
    {
        std::cerr << "FMM_Resources<EvalPtType>::finished() detected an error whilst waiting for event on a kernel: " << retval << std::endl; 
        throw std::exception();
    }
    return retval == CL_COMPLETE;
}

template <typename EvalPtType> 
void FMM_Resources<EvalPtType>::copy_out_results()
{
    for (unsigned int ii=0; ii < num_eval_pts; ++ii)
    {
        EvalPtType& ep = *(eval_pts_copy[ii]);
        const CL_VECTORIZED_FLOAT_TYPE& ocl_result = host_results[ii];

        ep.add_raw(static_cast<const float*>(&ocl_result.s[0]));
    }
}

template <typename EvalPtType> 
void FMM_Resources<EvalPtType>::enqueue(cl_command_queue& ocl_queue)
{
    int err;

    // Launch kernel
    //printf("clEnqueueNDRangeKernel (evaluation_kernel [%d x %d]) ... ", szGlobalWorkSize[0], szGlobalWorkSize[1], max_group_size);
    //std::cout << std::flush;
    err = clEnqueueNDRangeKernel(ocl_queue, m_evaluation_kernel, 2, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error %d in clEnqueueNDRangeKernel, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw EnqueueException();
    }
    //printf("done\n");

    size_t blk_results_kernel_width = BLOCK_RESULTS_COLLECTION_WIDTH;
    size_t blk_results_total_threads = blk_results_kernel_width * static_cast<size_t>(ceil(static_cast<double>(num_eval_pts_rounded) / blk_results_kernel_width));
    assert(blk_results_total_threads == num_eval_pts_rounded);
    //printf("clEnqueueNDRangeKernel (blk_results_kernel [%d]) ... ", blk_results_total_threads);
    //std::cout << std::flush;

    err = clEnqueueNDRangeKernel(ocl_queue, m_blk_results_kernel, 1, NULL, &blk_results_total_threads, &blk_results_kernel_width, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error %d in clEnqueueNDRangeKernel, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw EnqueueException();
    }
    //printf("done\n");

    //printf("clEnqueueReadBuffer...");
    err = clEnqueueReadBuffer(ocl_queue, cmDevFinalResults, CL_FALSE, 0, results_size, host_results.get(), 0, NULL, &m_event);
    if (err != CL_SUCCESS)
    {
        printf("Error %d in clEnqueueReadBuffer, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw EnqueueException();
    }
    //printf("done\n");

    return;
}

template <typename EvalPtType> 
const unsigned char* FMM_Resources<EvalPtType>::get_opencl_source() 
{ 
    return static_cast<const unsigned char*>(opencl_fmm_kernels_cl); 
}

template <typename EvalPtType> 
size_t FMM_Resources<EvalPtType>::device_mem_reqd() const
{
    size_t total=0;
    total += sizeof(cl_float4)*num_eval_pts_rounded;
    total += sizeof(cl_float4)*num_charges_rounded;
    total += block_results_size;
    total += results_size;
    return total;
}

template <typename EvalPtType> 
void FMM_Resources<EvalPtType>::allocate_device_mem(const cl_context& ocl_context)
{

    //printf("Allocating device memory...");
    int error;
    cmDevEvalPts = clCreateBuffer(ocl_context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, sizeof(cl_float4)*num_eval_pts_rounded, host_pts.get(), &error);
    if (error != CL_SUCCESS) {
        printf("Failed to create a device buffer, Error %d @ Line %u in file %s !!!\n\n", error, __LINE__, __FILE__);
        throw AllocationException();
    }

    cmDevCharges = clCreateBuffer(ocl_context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, sizeof(cl_float4)*num_charges_rounded, host_charges.get(), &error);
    if (error != CL_SUCCESS) {
        printf("Failed to create a device buffer, Error %d @ Line %u in file %s !!!\n\n", error, __LINE__, __FILE__);
        throw AllocationException();
    }
    cmDevBlockResults = clCreateBuffer(ocl_context, CL_MEM_READ_WRITE, block_results_size, NULL, &error);
    if (error != CL_SUCCESS) {
        printf("Failed to create a device buffer, Error %d @ Line %u in file %s !!!\n\n", error, __LINE__, __FILE__);
        throw AllocationException();
    }

    cmDevFinalResults = clCreateBuffer(ocl_context, CL_MEM_READ_WRITE, results_size, NULL, &error);
    if (error != CL_SUCCESS) {
        printf("Failed to create a device buffer, Error %d @ Line %u in file %s !!!\n\n", error, __LINE__, __FILE__);
        throw AllocationException();
    }

    //printf("done\n");

    return;
}

template <typename EvalPtType> 
void FMM_Resources<EvalPtType>::init_kernel(const cl_context& ocl_context,
                                const cl_program& ocl_program)
{
    const std::string& fmm_kernel = FMM_Derivatives_Level<EvalPtType>::fmm_kernel_name;
    const std::string& block_results_kernel = FMM_Derivatives_Level<EvalPtType>::block_results_kernel_name;

    int err;
    // Create the evaluation kernel
    //printf("clCreateKernel (FMM_explicit_kernel)...");
    m_evaluation_kernel = clCreateKernel(ocl_program, fmm_kernel.c_str(), &err);
    if (err != CL_SUCCESS)
    {
        printf("Error %d in clCreateKernel, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw std::exception();
    }
    //printf("done\n");

    // Set the Argument values
    //printf("clSetKernelArg...");
    err  = clSetKernelArg(m_evaluation_kernel, 0, sizeof(cl_mem), (void*)&cmDevEvalPts);
    err |= clSetKernelArg(m_evaluation_kernel, 1, sizeof(cl_mem), (void*)&cmDevCharges);
    err |= clSetKernelArg(m_evaluation_kernel, 2, sizeof(cl_mem), (void*)&cmDevBlockResults);
    err |= clSetKernelArg(m_evaluation_kernel, 3, sizeof(cl_uint), (void*)&num_eval_pts_rounded);
    err |= clSetKernelArg(m_evaluation_kernel, 4, sizeof(cl_uint), (void*)&num_charges_rounded);
    err |= clSetKernelArg(m_evaluation_kernel, 5, sizeof(cl_uint), (void*)&num_loops);
    err |= clSetKernelArg(m_evaluation_kernel, 6, sizeof(cl_float), (void*)&kappa);
    err |= clSetKernelArg(m_evaluation_kernel, 7, sizeof(CL_VECTORIZED_FLOAT_TYPE)*xdim*ydim, NULL);

    if (err != CL_SUCCESS)
    {
        printf("Error %d in clSetKernelArg, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw std::exception();
    }
    //printf("done\n");

    // Create the block-results collection kernel
    //printf("clCreateKernel (processBlockResults)...");
    m_blk_results_kernel = clCreateKernel(ocl_program, block_results_kernel.c_str(), &err);
    if (err != CL_SUCCESS)
    {
        printf("Error %d in clCreateKernel, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw std::exception();
    }
    //printf("done\n");

    // Set the Argument values
    //printf("clSetKernelArg...");
    cl_uint blk_xdim = szGlobalWorkSize[0];
    cl_uint blk_ydim = szGlobalWorkSize[1] / szLocalWorkSize[1];

    err  = clSetKernelArg(m_blk_results_kernel, 0, sizeof(cl_mem), (void*)&cmDevBlockResults);
    err |= clSetKernelArg(m_blk_results_kernel, 1, sizeof(cl_mem), (void*)&cmDevFinalResults);
    err |= clSetKernelArg(m_blk_results_kernel, 2, sizeof(cl_uint), (void*)&blk_xdim);
    err |= clSetKernelArg(m_blk_results_kernel, 3, sizeof(cl_uint), (void*)&blk_ydim);

    if (err != CL_SUCCESS)
    {
        printf("Error %d in clSetKernelArg, Line %u in file %s !!!\n\n", err, __LINE__, __FILE__);
        throw std::exception();
    }
    //printf("done\n");

    return;
}

#endif /* ifdef OPENCL */
#endif /* OPENCL_RESOURCES_H_ */
