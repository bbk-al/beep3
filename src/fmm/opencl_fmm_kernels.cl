// This is the OpenCL file for the FMM explicit interactions
// You should run xxd -i over this file to create a C-compatible
// header file, from whence it can be loaded into OpenCL.  NB:
// xxd produces non-null terminated arrays, so you might need
// to stick "0x00" onto the end of the array in the header file.

#define OPENCL_LOCAL_GROUP_DIM_EVAL_PTS 8
#define OPENCL_LOCAL_GROUP_DIM_CHARGES  8

#define PRECIS 1e-12

// OpenCL Kernel Function for FMM explicit neighbour interactions
__kernel void FMM_explicit_kernel(__global const float4* eval_pts,
                                __global const float4* charges,
                                __global float4* block_results,
                                const uint num_eval_pts,
                                const uint num_charges,
                                const uint num_loops,
                                const float beta,
                                __local float4 shr[OPENCL_LOCAL_GROUP_DIM_EVAL_PTS*OPENCL_LOCAL_GROUP_DIM_CHARGES])

{

    // get index into global data array
    uint x_id = get_global_id(0);
    uint y_id = get_global_id(1);

    uint x_local = get_local_id(0);
    uint y_local = get_local_id(1);

    // use shared_memory to hold per-thread-block results
    __local float4* thread_result = &(shr[x_local*OPENCL_LOCAL_GROUP_DIM_CHARGES + y_local]);

    // global data (don't worry, we promise not to dereference these
    // for the threads which are out of bounds.... honest)
    __global const float4* ch = charges + (y_id*num_loops);
    __global const float4* ep = eval_pts + x_id;
    
    float4 partial_results;
    partial_results.s0123 = 0;
    
    for (uint ctr=0; ctr < num_loops; ++ctr)
    {
        float rx = ((*ep).x - (*ch).x);
        float ry = ((*ep).y - (*ch).y);
        float rz = ((*ep).z - (*ch).z);
        float d_squared = rx*rx + ry*ry + rz*rz;

        // don't bother if this (point) charge is at the evaluation position
        float charge_mag = (*ch).w;
        if (d_squared > PRECIS) 
        {
            // add screened coulomb potential
            float d = native_sqrt(d_squared);
            float pot_term = charge_mag*native_exp(-beta*d) / d;
            float field_term = -pot_term*(1.0 + d*beta) / d_squared;
            partial_results.s0 += (pot_term);
            partial_results.s1 += (field_term * rx);
            partial_results.s2 += (field_term * ry);
            partial_results.s3 += (field_term * rz);
        }
        ++ch;// += OPENCL_LOCAL_GROUP_DIM_CHARGES;
    }
    (*thread_result) = partial_results;
    
    // Perform parallel reduction to add partial results
    for (uint stride = get_local_size(1)/2; stride > 0; stride /= 2)
    {
        // Synchronize to make sure each work-item is done updating
        // shared memory
        barrier(CLK_LOCAL_MEM_FENCE);

        // Only the first work-items in the work-group add elements
        // together
        if (y_local < stride) {
            // Add two elements from the "partialDotProduct" array
            // and store the result in partialDotProduct[index]
            *thread_result += shr[x_local*OPENCL_LOCAL_GROUP_DIM_CHARGES + y_local + stride];
        }
    }
    
    // Write the result of the reduction to global memory
    if (y_local == 0)
    {
        block_results[get_group_id(1)*get_global_size(0) + get_global_id(0)] = *thread_result; // grrrrr
    }

    return;
}

// Reduction of block results
__kernel void processBlockResults(__global const float4* block_results,
                                __global float4* final_results,
                                const unsigned int xdim,
                                const unsigned int ydim)
{
    // xdim is the number of block results in the x-direction -- gives the pitch
    // of the block results in memory
    // ydim is the number of block results in the y-direction -- loop over these

    const unsigned int xctr = get_global_id(0);

    float4 accum;
    accum.s0123 = 0;
    
    // TODO:: pre-load data into shared memory?
    for (unsigned int yctr=0; yctr < ydim; ++yctr)
    {
        // get result block for this eval pt, for the y'th work group
        accum += block_results[yctr*xdim + xctr];
    }

    final_results[xctr] = accum;

    return;
}

// OpenCL Kernel Function for FMM explicit neighbour interactions
__kernel void FMM_explicit_kernel_2nd_derivs(__global const float4* eval_pts,
                                             __global const float4* charges,
                                             __global float16* block_results,
                                             const uint num_eval_pts,
                                             const uint num_charges,
                                             const uint num_loops,
                                             const float beta,
                                             __local float16 shr[OPENCL_LOCAL_GROUP_DIM_EVAL_PTS*OPENCL_LOCAL_GROUP_DIM_CHARGES])

{

    // get index into global data array
    uint x_id = get_global_id(0);
    uint y_id = get_global_id(1);

    uint x_local = get_local_id(0);
    uint y_local = get_local_id(1);

    // use shared_memory to hold per-thread-block results
    __local float16* thread_result = &(shr[x_local*OPENCL_LOCAL_GROUP_DIM_CHARGES + y_local]);

    // global data (don't worry, we promise not to dereference these
    // for the threads which are out of bounds.... honest)
    __global const float4* ch = charges + (y_id*num_loops);
    __global const float4* ep = eval_pts + x_id;

    float16 partial_results;
    partial_results.lo = 0;
    partial_results.hi = 0;
    
    for (uint ctr=0; ctr < num_loops; ++ctr)
    {
        float rx = ((*ep).x - (*ch).x);
        float ry = ((*ep).y - (*ch).y);
        float rz = ((*ep).z - (*ch).z);
        float d_squared = rx*rx + ry*ry + rz*rz;

        // don't bother if this (point) charge is at the evaluation position
        float charge_mag = (*ch).w;
        if (d_squared > PRECIS) 
        {
            // add screened coulomb potential
            float d = native_sqrt(d_squared);
            float pot_term = charge_mag*native_exp(-beta*d) / d;
            float field_term = -pot_term*(1.0 + d*beta) / d_squared;
            
            partial_results.s0 += (pot_term);
            partial_results.s1 += (field_term * rx);
            partial_results.s2 += (field_term * ry);
            partial_results.s3 += (field_term * rz);
            
            // SECOND DERIVS
            // diagonal terms
            const float dbeta = d*beta;
            const float d_beta_squared = dbeta*dbeta;
            const float d4 = d_squared * d_squared;
            
            partial_results.s4 += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*rx*rx) - rx*rx*d_beta_squared) / d4;
            partial_results.s8 += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*ry*ry) - ry*ry*d_beta_squared) / d4;
            partial_results.sC += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*rz*rz) - rz*rz*d_beta_squared) / d4;

            // off-diagonals
            const float f2_mixed_term = pot_term*(3.0 + 3.0*dbeta + d_beta_squared)/d4;
            
            partial_results.s5  += f2_mixed_term * rx*ry;
            partial_results.s6  += f2_mixed_term * rx*rz;
            partial_results.s7  += f2_mixed_term * ry*rx;
            partial_results.s9  += f2_mixed_term * ry*rz;
            partial_results.sA += f2_mixed_term * rz*rx;
            partial_results.sB += f2_mixed_term * rz*ry;
            
        }
        ++ch;// += OPENCL_LOCAL_GROUP_DIM_CHARGES;
    }
    (*thread_result) = partial_results;

    // Perform parallel reduction to add partial results
    for (uint stride = get_local_size(1)/2; stride > 0; stride /= 2)
    {
        // Synchronize to make sure each work-item is done updating
        // shared memory
        barrier(CLK_LOCAL_MEM_FENCE);

        // Only the first work-items in the work-group add elements
        // together
        if (y_local < stride) {
            // Add two elements from the "partialDotProduct" array
            // and store the result in partialDotProduct[index]
            *thread_result += shr[x_local*OPENCL_LOCAL_GROUP_DIM_CHARGES + y_local + stride];
        }
    }
    
    // Write the result of the reduction to global memory
    if (y_local == 0)
    {
        block_results[get_group_id(1)*get_global_size(0) + get_global_id(0)] = *thread_result; // grrrrr
    }

    return;
}

// Reduction of block results
__kernel void processBlockResults_2nd_derivs(__global const float16* block_results,
                                            __global float16* final_results,
                                            const unsigned int xdim,
                                            const unsigned int ydim)
{
    // xdim is the number of block results in the x-direction -- gives the pitch
    // of the block results in memory
    // ydim is the number of block results in the y-direction -- loop over these

    const unsigned int xctr = get_global_id(0);

    float16 accum;
    accum.lo = 0;
    accum.hi = 0;

    // TODO:: pre-load data into shared memory?
    for (unsigned int yctr=0; yctr < ydim; ++yctr)
    {
        // get result block for this eval pt, for the y'th work group
        accum += block_results[yctr*xdim + xctr];
    }

    final_results[xctr] = accum;

    return;
}
