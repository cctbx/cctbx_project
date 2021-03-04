#ifndef SIMTBX_GPU_DETECTOR_CUH
#define SIMTBX_GPU_DETECTOR_CUH

#include <simtbx/nanoBragg/nanotypes.h>

__global__ void scale_array_CUDAKernel(double scale_factor, double * lhs, int array_size){
  const int total_pixels = array_size;
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
       pixIdx < total_pixels; pixIdx += stride) {
    /* position in pixel array */
    const int j = pixIdx;
    lhs[j] = lhs[j] * scale_factor;
    }
}

#endif // SIMTBX_GPU_DETECTOR_CUH
