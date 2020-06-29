#ifndef CUDA_STRUCT_H
#define CUDA_STRUCT_H

/*
  Header for defining a struct to house pointers on GPU
*/

#include <simtbx/nanoBragg/nanotypes.h>

#ifndef CUDAREAL
#define CUDAREAL double
#endif

#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

struct hklParams {
        int hkls;
        int h_min;
        int h_max;
        int h_range;
        int k_min;
        int k_max;
        int k_range;
        int l_min;
        int l_max;
        int l_range;
};

struct cudaPointers {

  int cu_spixels, cu_fpixels, cu_roi_xmin, cu_roi_xmax, cu_roi_ymin,
    cu_roi_ymax, cu_oversample, cu_point_pixel;

  CUDAREAL cu_pixel_size, cu_subpixel_size;
  int cu_steps;
  CUDAREAL cu_detector_thickstep, cu_detector_thicksteps, cu_detector_thick,
    cu_detector_mu;
  CUDAREAL * cu_sdet_vector, * cu_fdet_vector, * cu_odet_vector, * cu_pix0_vector;
  int cu_curved_detector;
  CUDAREAL cu_distance, cu_close_distance;
  CUDAREAL * cu_beam_vector;
  CUDAREAL cu_Xbeam, cu_Ybeam, cu_dmin, cu_phi0, cu_phistep;
  int cu_phisteps;
  CUDAREAL * cu_spindle_vector;
  int cu_sources;
  CUDAREAL * cu_source_X, * cu_source_Y, * cu_source_Z, * cu_source_I, * cu_source_lambda,
    * cu_a0, * cu_b0, * cu_c0;
  shapetype cu_xtal_shape;
  CUDAREAL cu_mosaic_spread;
  int cu_mosaic_domains;
  CUDAREAL * cu_mosaic_umats;
  CUDAREAL cu_Na, cu_Nb, cu_Nc, cu_V_cell, cu_water_size, cu_water_F,
    cu_water_MW, cu_r_e_sqr, cu_fluence, cu_Avogadro;
  int cu_integral_form;
  CUDAREAL cu_default_F;
  int cu_interpolate;
  CUDAREAL * cu_Fhkl;
  hklParams * cu_FhklParams;
  int cu_nopolar;
  CUDAREAL * cu_polar_vector;
  CUDAREAL cu_polarization, cu_fudge;
  int unsigned short * cu_maskimage;
  float * cu_floatimage /*out*/,  * cu_omega_reduction/*out*/,
    * cu_max_I_x_reduction/*out*/, * cu_max_I_y_reduction /*out*/;
  bool * cu_rangemap /*out*/;
};

#endif
