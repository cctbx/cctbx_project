#ifndef SIMTBX_GPU_SIMULATION_CUH
#define SIMTBX_GPU_SIMULATION_CUH

#include <simtbx/nanoBragg/nanotypes.h>
using simtbx::nanoBragg::shapetype;
using simtbx::nanoBragg::hklParams;

__global__ void nanoBraggSpotsCUDAKernel(int spixels, int fpixels, int roi_xmin, int roi_xmax,
    int roi_ymin, int roi_ymax, int oversample, int point_pixel,
    CUDAREAL pixel_size, CUDAREAL subpixel_size, int steps, CUDAREAL detector_thickstep,
    int detector_thicksteps, CUDAREAL detector_thick, CUDAREAL detector_mu,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector,
    const CUDAREAL * __restrict__ pix0_vector, int curved_detector, CUDAREAL distance,
    CUDAREAL close_distance, const CUDAREAL * __restrict__ beam_vector,
    CUDAREAL Xbeam, CUDAREAL Ybeam, CUDAREAL dmin, CUDAREAL phi0, CUDAREAL phistep,
    int phisteps, const CUDAREAL * __restrict__ spindle_vector, int sources,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_I, const CUDAREAL * __restrict__ source_lambda,
    const CUDAREAL * __restrict__ a0, const CUDAREAL * __restrict__ b0,
    const CUDAREAL * __restrict c0, shapetype xtal_shape, CUDAREAL mosaic_spread,
    int mosaic_domains, const CUDAREAL * __restrict__ mosaic_umats,
    CUDAREAL Na, CUDAREAL Nb,
    CUDAREAL Nc, CUDAREAL V_cell,
    CUDAREAL water_size, CUDAREAL water_F, CUDAREAL water_MW, CUDAREAL r_e_sqr,
    CUDAREAL fluence, CUDAREAL Avogadro, CUDAREAL spot_scale, int integral_form,
    CUDAREAL default_F,
    int interpolate, const CUDAREAL * __restrict__ Fhkl,
    const hklParams * __restrict__ Fhklparams, int nopolar,
    const CUDAREAL * __restrict__ polar_vector, CUDAREAL polarization, CUDAREAL fudge,
    const int unsigned short * __restrict__ maskimage, float * floatimage /*out*/,
    float * omega_reduction/*out*/, float * max_I_x_reduction/*out*/,
    float * max_I_y_reduction /*out*/, bool * rangemap);

__global__ void nanoBraggSpotsInitCUDAKernel(int spixels, int fpixesl, float * floatimage, float * omega_reduction,
                float * max_I_x_reduction,
                float * max_I_y_reduction, bool * rangemap);

__global__ void add_array_CUDAKernel(double * lhs, float * rhs, int array_size){
  const int total_pixels = array_size;
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx =
      (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
      pixIdx < total_pixels; pixIdx += stride) {
    const int j = pixIdx; /* position in pixel array */
    lhs[j] = lhs[j] + (double)rhs[j]; // specifically add low precision to high precision array
    }
  }

__global__ void add_background_CUDAKernel(int sources, int nanoBragg_oversample,
    CUDAREAL pixel_size, int spixels, int fpixels, int detector_thicksteps,
    CUDAREAL detector_thickstep, CUDAREAL detector_attnlen,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector, const CUDAREAL * __restrict__ pix0_vector,
    CUDAREAL close_distance, int point_pixel, CUDAREAL detector_thick,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_lambda, const CUDAREAL * __restrict__ source_I,
    int stols, const CUDAREAL * stol_of, const CUDAREAL * Fbg_of,
    int nopolar, CUDAREAL polarization, const CUDAREAL * __restrict__ polar_vector,
    CUDAREAL r_e_sqr, CUDAREAL fluence, CUDAREAL amorphous_molecules,
    float * floatimage);

#endif // SIMTBX_GPU_SIMULATION_CUH
