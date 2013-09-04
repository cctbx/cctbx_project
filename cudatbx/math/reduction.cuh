#ifndef CUDATBX_REDUCTION_CUH
#define CUDATBX_REDUCTION_CUH

#include <cudatbx/cuda_base.cuh>

namespace cudatbx {
namespace math {

  /* ==========================================================================
     Weighted Sum Kernel

     This kernel calculates the weighted sum for an array of values by
     breaking the array up and doing pairwise additions within a block. The 
     first thread will add up the sums from each block to calculate the final
     sum.

     This implementation is not optimal, but is easy to follow, and follows
     some of the basic ideas from the reduction example from the NVIDIA SDK.

       values - pointer to array of values
       weights - pointer to array of weights
       n_values - number of values, equals number of weights
       sum - pointer for storing final sum
     --------------------------------------------------------------------------
  */

  template <typename floatType>
  __global__ void weighted_sum_kernel
  (const floatType * values, const floatType * weights, const int n_values,
   floatType * sum) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    // set initial sum to zero
    if (i == 0) {
      *sum = floatType(0.0);
    }
    
    // dynamically allocate shared memory (# of threads * sizeof(floatType))
    extern __shared__ floatType workspace[];
    
    // transfer chunk from global memory to shared memory and apply weight
    if (i < n_values) {
      workspace[tid] = weights[i] * values[i];
    } else {
      workspace[tid] = floatType(0.0);
    }
    __syncthreads();
    
    // do pairwise summing until the first element contains the block sum
    for (int cycle=blockDim.x/2; cycle>0; cycle=cycle/2) {
      if (tid < cycle) {
        workspace[tid] += workspace[tid + cycle];
      }
      __syncthreads();
    }
    
    // add sum from block to total sum
    if (tid == 0) {
      atomicAdd(sum,workspace[0]);
    }
  }
  // ==========================================================================
}
}
#endif // CUDATBX_REDUCTION_CUH
