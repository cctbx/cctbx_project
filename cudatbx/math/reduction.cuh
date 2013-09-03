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
    
    int cycle, current_j;
    floatType s = floatType(0.0);
    
    // dynamically allocate shared memory (# of threads * sizeof(floatType))
    extern __shared__ floatType workspace[];
    
    for (int j=0; j<n_values; j+=blockDim.x) {
      // transfer chunk from global memory to shared memory and apply weight
      current_j = j + threadIdx.x;
      if (current_j < n_values) {
        workspace[threadIdx.x] = weights[current_j] * values[current_j];
      } else {
        workspace[threadIdx.x] = floatType(0.0);
      }
      __syncthreads();
      
      // do pairwise summing until the first element contains the block sum
      cycle = blockDim.x/2;
      while (cycle >= 1) {
        if (threadIdx.x < cycle) {
          workspace[threadIdx.x] = workspace[threadIdx.x] +
            workspace[threadIdx.x + cycle];
        }
        cycle = cycle/2;
        __syncthreads();
      }
      
      // add sum from block to total sum
      if (threadIdx.x == 0) {
        s += workspace[0];
      }
      __syncthreads();
    }

    // transfer total sum back to global memory (array or single value)
    if (threadIdx.x == 0) {
      *sum = s;
    }
  }
  // ==========================================================================
}
}
#endif // CUDATBX_REDUCTION_CUH
