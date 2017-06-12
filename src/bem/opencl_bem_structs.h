#ifndef __OPENCL_BEM_STRUCTS_H_
#define __OPENCL_BEM_STRUCTS_H_

//
// WARNING: do not change these without changing the corresponding 
// definitions within opencl_bem_kernels.cl !!!
// NB: the OpenCL code is included in the application via a static 
// char array, created using xxd -i and resulting in the hex-dump 
// opencl_bem_kernels.cl.h. So changing stuff here means changing it
// in opencl_bem_kernels.cl, then regenerating opencl_bem_kernels.cl.h
// using xxd (unix tool for hex dumping files).
//

#define OPENCL_LOCAL_GROUP_DIM_SOURCES 4
#define OPENCL_LOCAL_GROUP_DIM_TARGETS 8
#define BEM_BLOCK_RESULTS_COLLECTION_WIDTH 128

typedef struct ocl_qp
{
    cl_float ptx;
    cl_float pty;
    cl_float ptz;
    cl_float nx;
    cl_float ny;
    cl_float nz;
    cl_float weight;
    cl_float dummy;

} OpenCL_QuadPoint;


typedef struct ocl_np
{
    cl_uint idx;
    cl_float inverse_dielectric_ratio;
    cl_uint first_quad_idx;
    cl_uint num_quad_points;

} OpenCL_NodePatch;

typedef struct ocl_lint
{
    cl_uint src_global_idx;
    cl_uint targ_global_idx;
    cl_float Apt;
    cl_float Bpt;
    cl_float Cpt;
    cl_float Dpt;

} OpenCL_LocalIntegrationResult;

#endif // __OPENCL_BEM_STRUCTS_H_
