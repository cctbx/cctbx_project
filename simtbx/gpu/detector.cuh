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

__global__ void get_active_pixel_selection_CUDAKernel(double * dst, const std::size_t * __restrict__ indices,
                                                   const double * __restrict__ src, int total_pixels){
  /* select 'total_pixels' from array 'src' selecting the 'indices' in that order, and deposit the results in 'dst'
  */
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
       pixIdx < total_pixels; pixIdx += stride) {
    /* position in pixel array */
    const int j = pixIdx;
    dst[j] = src[indices[j]];
    }
}

#endif // SIMTBX_GPU_DETECTOR_CUH
