/* opencl_bem_kernels.cl

    part of BEEP: Boundary Element Electrostatics Program
    (c) David Fallaize 2011

    This file is the OpenCL code for the BEM.
    In the application code this entire file is
    read in as a static char array, produced
    by running xxd over this file. So if you
    change this then make sure you re-generate
    the corresponding opencl_bem_kernels.cl.h 
    hex dump of this file.
    
    (By the way you might need to append 0x00 to 
    the array in the header file, as xxd does not
    null-terminate the char* array.)

*/

#define OPENCL_LOCAL_GROUP_DIM_SOURCES 4
#define OPENCL_LOCAL_GROUP_DIM_TARGETS 8

typedef struct ocl_qp
{
    float ptx;
    float pty;
    float ptz;
    float nx;
    float ny;
    float nz;
    float weight;
    float dummy;

} OpenCL_QuadPoint;

typedef struct ocl_np
{
    uint idx;
    float inverse_dielectric_ratio;
    uint first_quad_idx;
    uint num_quad_points;

} OpenCL_NodePatch;

typedef struct ocl_lint
{
    uint src_global_idx;
    uint targ_global_idx;
    float Apt;
    float Bpt;
    float Cpt;
    float Dpt;

} OpenCL_LocalIntegrationResult;

#define ONE_OVER_4PI 0.079577471545947673f;

float length2(const float4* a)
{
    return ( ((*a).x * (*a).x) + ((*a).y * (*a).y) + ((*a).z * (*a).z) );
}

float dot_product(const float4* a, const float4* b)
{
    // calculate: a.dot.b and put result in result
    return ((*a).x * (*b).x) + ((*a).y * (*b).y) + ((*a).z * (*b).z);
}

void evaluate_local_BEM_kernels(float kappa,
                                float inv_epsilon,
                                __global const OpenCL_QuadPoint* src_qp,
                                __global const OpenCL_QuadPoint* qp,
                                float* A,
                                float* B,
                                float* C,
                                float* D)
{
    const float sigm = 1e-9;

    float4 dx;
    dx.x = qp->ptx;
    dx.y = qp->pty;
    dx.z = qp->ptz;
    dx.x -= src_qp->ptx;
    dx.y -= src_qp->pty;
    dx.z -= src_qp->ptz;

    float4 normal;
    normal.x = qp->nx;
    normal.y = qp->ny;
    normal.z = qp->nz;

    float4 n0;
    n0.x = src_qp->nx;
    n0.y = src_qp->ny;
    n0.z = src_qp->nz;

    float wt = qp->weight;

    // use inverse multiplication instead of division
    const float r2 = length2(&dx);
    if (r2 <= sigm) { return; }
    const float ir = rsqrt(r2);
    const float r = native_recip(ir);
    const float ir2 = native_recip(r2); // reciprocal native func?
    const float ir3 = ir2*ir;
    const float ir4 = ir2*ir2;
    const float ir5 = ir4*ir;

    const float gr0 = ONE_OVER_4PI;
    const float gr1 = gr0*ir;
    const float gr3 = gr0*ir3;
    const float ur0=native_exp(-kappa*r) * ONE_OVER_4PI;
    const float ur1=ur0*ir;
    const float ur2=ur0*ir2;
    const float ur3=ur0*ir3;

    *A += wt * (gr1 - ur1);
    
    {
        float4 gd;
        gd.x = -gr3 * dx.x;
        gd.y = -gr3 * dx.y;
        gd.z = -gr3 * dx.z;

        float4 ud;
        ud.x = -dx.x * (ur3 + kappa*ur2);
        ud.y = -dx.y * (ur3 + kappa*ur2);
        ud.z = -dx.z * (ur3 + kappa*ur2);

        float4 tmp;
        tmp.x = gd.x * inv_epsilon - ud.x;
        tmp.y = gd.y * inv_epsilon - ud.y;
        tmp.z = gd.z * inv_epsilon - ud.z;

        *B += wt * dot_product(&tmp, &normal);

        tmp.x = ud.x * inv_epsilon - gd.x;
        tmp.y = ud.y * inv_epsilon - gd.y;
        tmp.z = ud.z * inv_epsilon - gd.z;
        *C += wt * dot_product(&tmp, &n0);
    }
    
    float Dtmp = 0.0;
    if (kappa != 0.0)
    {

        const float ur3ur2=ur3+kappa*ur2;
        const float gr5 = gr0*ir5;
        const float ur4=ur0*ir4;
        const float ur5=ur0*ir5;
   
        float pur4ur3=kappa*(3.0*ur4+kappa*ur3)+3.0*ur5;

        Dtmp += n0.x * normal.x * ((3.0*dx.x*dx.x*gr5 - gr3) - (dx.x*dx.x*pur4ur3 - ur3ur2));
        Dtmp += n0.y * normal.x * ((3.0*dx.x*dx.y*gr5)       - (dx.x*dx.y*pur4ur3));
        Dtmp += n0.z * normal.x * ((3.0*dx.x*dx.z*gr5)       - (dx.x*dx.z*pur4ur3));
        Dtmp += n0.x * normal.y * ((3.0*dx.y*dx.x*gr5)       - (dx.y*dx.x*pur4ur3));
        Dtmp += n0.y * normal.y * ((3.0*dx.y*dx.y*gr5 - gr3) - (dx.y*dx.y*pur4ur3 - ur3ur2));
        Dtmp += n0.z * normal.y * ((3.0*dx.y*dx.z*gr5)       - (dx.y*dx.z*pur4ur3));
        Dtmp += n0.x * normal.z * ((3.0*dx.z*dx.x*gr5)       - (dx.z*dx.x*pur4ur3));
        Dtmp += n0.y * normal.z * ((3.0*dx.z*dx.y*gr5)       - (dx.z*dx.y*pur4ur3));
        Dtmp += n0.z * normal.z * ((3.0*dx.z*dx.z*gr5 - gr3) - (dx.z*dx.z*pur4ur3 - ur3ur2));

    }
    *D -= (wt * Dtmp * inv_epsilon);
}

// OpenCL Kernel Function for BEM integrations
__kernel void BEM_kernels(__global const OpenCL_NodePatch* srcs,
                        __global const OpenCL_NodePatch* targs,
                        __global const OpenCL_QuadPoint* quad_points,
                        __global OpenCL_LocalIntegrationResult* results,
                        const uint num_source_patches,
                        const uint num_target_patches,
                        const uint num_quad_points,
                        const float kappa)
{

    // get index into global data array
    uint iGID_src = get_global_id(0);
    uint iGID_targ = get_global_id(1);
    if (iGID_src >= num_source_patches) { return; }
    if (iGID_targ >= num_target_patches) { return; }

    __global const OpenCL_NodePatch* s = srcs + iGID_src;
    __global const OpenCL_NodePatch* t = targs + iGID_targ;
    __global OpenCL_LocalIntegrationResult* res = results + iGID_src*num_target_patches + iGID_targ;

    __global const OpenCL_QuadPoint* src_qp = quad_points + s->first_quad_idx;

    res->src_global_idx = s->idx;
    res->targ_global_idx = t->idx;
    float A_total=0,B_total=0,C_total=0,D_total=0;
    
    if (s->idx != t->idx) {
    
        const uint max_ii = s->num_quad_points;
        for (uint ii=0; ii < max_ii; ++ii)
        {

            float A,B,C,D;
            A=0;
            B=0;
            C=0;
            D=0;

            __global const OpenCL_QuadPoint* qp = quad_points + t->first_quad_idx;
            const uint max_jj = t->num_quad_points;
            const float inv_epsilon = t->inverse_dielectric_ratio;
            for (uint jj=0; jj < max_jj; ++jj)
            {
                evaluate_local_BEM_kernels(kappa, inv_epsilon, src_qp, qp, &A, &B, &C, &D);
                ++qp;
            }
            
            const float wt = src_qp->weight;
            A_total += A * wt;
            B_total += B * wt;
            C_total += C * wt;
            D_total += D * wt;

            ++src_qp;
        
        }
    }
    
    res->Apt = A_total;
    res->Bpt = B_total;
    res->Cpt = C_total;
    res->Dpt = D_total;


}

// OpenCL Kernel Function for BEM integrations
__kernel void BEM_kernels_instant(__global const float *vals,
                                __global const OpenCL_NodePatch* srcs,
                                __global const OpenCL_NodePatch* targs,
                                __global const OpenCL_QuadPoint* quad_points,
                                __global float2* block_results,
                                const uint num_source_patches,
                                const uint num_target_patches,
                                const uint num_quad_points,
                                const float kappa,
                                __local float2 shr[OPENCL_LOCAL_GROUP_DIM_SOURCES*OPENCL_LOCAL_GROUP_DIM_TARGETS])
{

    // get index into global data array
    uint iGID_src = get_global_id(0);
    uint iGID_targ = get_global_id(1);
    bool working_thread = (iGID_src < num_source_patches && iGID_targ < num_target_patches);

    // for thread-block shared memory storage
    uint x_local = get_local_id(0);
    uint y_local = get_local_id(1);

    __local float2* thread_result = &(shr[x_local*OPENCL_LOCAL_GROUP_DIM_TARGETS + y_local]);
    (*thread_result).x = 0;
    (*thread_result).y = 0;

    if (working_thread)
    {
        __global const OpenCL_NodePatch* s = srcs + iGID_src;
        __global const OpenCL_NodePatch* t = targs + iGID_targ;

        if (s->idx != t->idx) {

        __global const OpenCL_QuadPoint* src_qp = quad_points + s->first_quad_idx;

        __global const float *fvals = vals;
        __global const float *hvals = vals + num_target_patches;

        const float f = fvals[iGID_targ];
        const float h = hvals[iGID_targ];

        const uint max_ii = s->num_quad_points;
        for (uint ii=0; ii < max_ii; ++ii)
        {
            __global const OpenCL_QuadPoint* qp = quad_points + t->first_quad_idx;

            float A=0,B=0,C=0,D=0;

            const uint max_jj = t->num_quad_points;
            const float inv_epsilon = t->inverse_dielectric_ratio;
            
            for (uint jj=0; jj < max_jj; ++jj)
            {
                evaluate_local_BEM_kernels(kappa, inv_epsilon, src_qp, qp, &A, &B, &C, &D);
                ++qp;
            }

            // use shared_memory to hold per-thread-block results
            (*thread_result).x += (f*B - h*A)*src_qp->weight;
            (*thread_result).y += (f*D - h*C)*src_qp->weight;

            ++src_qp;
        }
        }
    }

    // Perform parallel reduction to add partial results
    for (uint stride = get_local_size(1)/2; stride > 0; stride /= 2)
    {
        // Synchronize to make sure each work-item is done updating
        // shared memory
        barrier(CLK_LOCAL_MEM_FENCE);

        // Only the first work-items in the work-group add elements
        // together
        if (y_local < stride) {
            (*thread_result) += shr[x_local*OPENCL_LOCAL_GROUP_DIM_TARGETS + y_local + stride];
        }
    }

    // Write the result of the reduction to global memory
    if (y_local == 0)
    {
        block_results[get_group_id(1)*get_global_size(0) + get_global_id(0)] = *thread_result; // grrrrr
    }


}

// Reduction of block results
__kernel void processBEMBlockResults(__global float2* block_results,
                                __global float2* final_results,
                                const unsigned int num_eval_pts,
                                const unsigned int xdim,
                                const unsigned int ydim)
{
    // xdim is the number of block results in the x-direction -- gives the pitch
    // of the block results in memory
    // ydim is the number of block results in the y-direction -- loop over these

    const unsigned int xctr = get_global_id(0);
    if (xctr >= num_eval_pts) { return; }

    float2 accum;
    accum.x=0;
    accum.y=0;

    for (unsigned int yctr=0; yctr < ydim; ++yctr)
    {
        // get result block for this eval pt, for the y'th work group
        accum += block_results[yctr*xdim + xctr];
    }

    final_results[xctr] = accum;

    return;
}

// OpenCL Kernel Function for Singular (same-patch) BEM integrations
__kernel void Singular_BEM_Kernels(__global const OpenCL_NodePatch* srcs,
                                    __global const OpenCL_QuadPoint* quad_points,
                                    __global OpenCL_LocalIntegrationResult* results,
                                    const uint num_source_patches,
                                    const uint num_quad_points,
                                    const float kappa)
{

    // get index into global data array
    uint iGID_src = get_global_id(0);
    if (iGID_src >= num_source_patches) { return; }

    __global const OpenCL_NodePatch* s = srcs + iGID_src;
    __global OpenCL_LocalIntegrationResult* res = results + iGID_src;

    __global const OpenCL_QuadPoint* src_qp = quad_points + s->first_quad_idx;
    __global const OpenCL_QuadPoint* targ_qp = quad_points + s->first_quad_idx;

    res->src_global_idx = s->idx;
    res->targ_global_idx = s->idx;

    float A_total=0,B_total=0,C_total=0,D_total=0,wt_normaliser=0;
    const float inv_epsilon = s->inverse_dielectric_ratio;
    
    const uint max_ii = s->num_quad_points;
    for (uint ii=0; ii < max_ii; ++ii)
    {

        float A,B,C,D;
        A=0;
        B=0;
        C=0;
        D=0;

        __global const OpenCL_QuadPoint* qp = targ_qp;
        for (int jj=0; jj < max_ii; ++jj)
        {
            if (ii == jj) { continue; }
            evaluate_local_BEM_kernels(kappa, inv_epsilon, src_qp, qp, &A, &B, &C, &D);
            ++qp;
        }
        
        const float wt = src_qp->weight;
        wt_normaliser += wt;
        A_total += A * wt;
        B_total += B * wt;
        C_total += C * wt;
        D_total += D * wt;

        ++src_qp;
    
    }
    res->src_global_idx = s->idx;
    res->targ_global_idx = s->idx;
    res->Apt = A_total / wt_normaliser;
    res->Bpt = B_total / wt_normaliser;
    res->Cpt = C_total / wt_normaliser;
    res->Dpt = D_total / wt_normaliser;


}

