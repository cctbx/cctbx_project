#ifndef SIMTBX_GPU_DETECTOR_CUH
#define SIMTBX_GPU_DETECTOR_CUH

#include <curand_kernel.h>
#include <curand.h>
#include <simtbx/nanoBragg/nanotypes.h>

__global__ void setup_random_states_CUDAKernel(curandState * state, int total_pixels){
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
       pixIdx < total_pixels; pixIdx += stride) {
    // different seed, same offset and sequence number
    curand_init(pixIdx, 0, 0, &state[pixIdx]);
  }
}

__global__ void noisify_CUDAKernel(curandState* states, double * lhs,
      double flicker_noise, double calibration_noise, double readout_noise,
      double quantum_gain, double adc_offset, int array_size, int global_shot_idx){

  const int total_pixels = array_size;
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
       pixIdx < total_pixels; pixIdx += stride) {

    /* position in pixel array */
    const int j = pixIdx;

    // the first `total_pixels` rand calls are dedicated to the calibration (should be the same for every shot)
    // the next 3 calls (2 normal devs, 1 poisson dev) should be different for every shot/pixel
    curandState localState = states[j];
    double calib_dev = (double) curand_normal(&localState);
    unsigned long long nskip = total_pixels + 3*total_pixels*global_shot_idx;
    skipahead(nskip, &localState);

    /* ideal photons/pixel */
    double photons0 = lhs[j];

    /* simulate 1/f noise in source */
    double flick_dev = (double) curand_normal(&localState);
    photons0 *= ( 1.0 + flicker_noise * flick_dev );

    /* simulate photon-counting error (assume calibration error is loss of photons, not electrons) */
    photons0 = (double) curand_poisson(&localState, photons0);

    photons0 *= ( 1.0 + calibration_noise * calib_dev );
    double adu = photons0*quantum_gain + adc_offset;
    double readout_dev = (double) curand_normal(&localState);
    adu += readout_noise * readout_dev;
    // end random code

    // TODO: dont copy state back, just skipahead instead ..
    //states[j] = localState;
    lhs[j] = adu;
  }
}

__global__ void offset_array_CUDAKernel(double* offset, double * lhs, int array_size){
  const int total_pixels = array_size;
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
       pixIdx < total_pixels; pixIdx += stride) {
    /* position in pixel array */
    const int j = pixIdx;
    lhs[j] = lhs[j] + offset[j];
    }
}

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
