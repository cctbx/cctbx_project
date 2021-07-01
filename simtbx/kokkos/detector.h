#ifndef SIMTBX_KOKKOS_DETECTOR_H
#define SIMTBX_KOKKOS_DETECTOR_H

// class is intended to capture several concepts:
//   1) a dxtbx multipanel detector
//   2) the associated experimental data
//   3) the associated simulated data, in process of being accumulated by kernel calls
//   4) mask data
//   5) possibly other metadata


#include "scitbx/array_family/shared.h"
#include "scitbx/array_family/flex_types.h"
#include "dxtbx/model/detector.h"
#include "dxtbx/model/beam.h"
#include "simtbx/nanoBragg/nanoBragg.h"
#include "kokkos_types.h"

#include <iostream>

namespace simtbx { namespace Kokkos {

namespace af = scitbx::af;

struct packed_metrology{
  packed_metrology(){/*printf("NO OPERATION")*/;};
  packed_metrology(dxtbx::model::Detector const &,dxtbx::model::Beam const &);
  packed_metrology(const simtbx::nanoBragg::nanoBragg& nB);
  void show() const;
  af::shared<double>sdet;
  af::shared<double>fdet;
  af::shared<double>odet;
  af::shared<double>pix0;
  af::shared<double>dists;
  af::shared<double>Xbeam;
  af::shared<double>Ybeam;
};

struct kokkos_detector{
  inline kokkos_detector(){printf("NO OPERATION, DEVICE NUMBER IS NEEDED");};
  kokkos_detector(const simtbx::nanoBragg::nanoBragg& nB);
  kokkos_detector(dxtbx::model::Detector const &, dxtbx::model::Beam const &);
  vector_double_t construct_detail(dxtbx::model::Detector const &);

  inline void show_summary(){
    std::cout << "Detector size: " << m_panel_count << " panel" << (m_panel_count>1)? "s" : "" <<std::endl;
    metrology.show();
  }
  void each_image_allocate_cuda();
  void scale_in_place_cuda(const double&);
  void write_raw_pixels_cuda(simtbx::nanoBragg::nanoBragg&);
  af::flex_double get_raw_pixels_cuda();
  void set_active_pixels_on_KOKKOS(af::shared<int>);
  af::shared<double> get_whitelist_raw_pixels_cuda(af::shared<std::size_t>);

  const dxtbx::model::Detector detector;

  // detector size and dimensions in pixels
  int m_panel_count = -1;
  int m_slow_dim_size = -1;
  int m_fast_dim_size = -1;
  int m_total_pixel_count = -1;

  // kokkos arrays
  vector_double_t m_accumulate_floatimage = vector_double_t("m_accumulate_floatimage", 0);
  double * cu_accumulate_floatimage; // pointer to GPU memory

  // per image input variables, pointers to GPU memory
  vector_ushort_t m_maskimage = vector_ushort_t("m_maskimage", 0);
  int unsigned short * cu_maskimage; // nanoBragg only manages single-panel, multipanel case must mask elsewhere

  // per kernel-call output variables, pointers to GPU memory
  vector_float_t m_omega_reduction = vector_float_t("m_omega_reduction", 0);
  vector_float_t m_max_I_x_reduction = vector_float_t("m_max_I_x_reduction", 0);
  vector_float_t m_max_I_y_reduction = vector_float_t("m_max_I_y_reduction", 0);
  float * cu_omega_reduction, * cu_max_I_x_reduction, * cu_max_I_y_reduction;

  vector_bool_t m_rangemap = vector_bool_t("m_rangemap", 0);
  bool * cu_rangemap;

  vector_float_t m_floatimage = vector_float_t("m_floatimage", 0);
  float * cu_floatimage; //most important output, stores the simulated pixel data

  // all-panel packed GPU representation of the multi-panel metrology
  vector_cudareal_t m_sdet_vector = vector_cudareal_t("m_sdet_vector", 0);
  vector_cudareal_t m_fdet_vector = vector_cudareal_t("m_fdet_vector", 0);
  vector_cudareal_t m_odet_vector = vector_cudareal_t("m_odet_vector", 0);
  vector_cudareal_t m_pix0_vector = vector_cudareal_t("m_pix0_vector", 0);
  vector_cudareal_t m_distance = vector_cudareal_t("m_distance", 0);
  vector_cudareal_t m_Xbeam = vector_cudareal_t("m_Xbeam", 0);
  vector_cudareal_t m_Ybeam = vector_cudareal_t("m_Ybeam", 0);
  CUDAREAL * cu_sdet_vector, * cu_fdet_vector;
  CUDAREAL * cu_odet_vector, * cu_pix0_vector;
  CUDAREAL * cu_distance, * cu_Xbeam, * cu_Ybeam;

  // all-panel whitelist of active pixels on GPU
  af::shared<int> active_pixel_list;
  int m_active_pixel_size;
  vector_int_t m_active_pixel_list = vector_int_t("m_active_pixel_list", 0);
  int * cu_active_pixel_list;

  packed_metrology const metrology;
};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_DETECTOR_H

