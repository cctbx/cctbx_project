#ifndef SIMTBX_KOKKOS_DETECTOR_H
#define SIMTBX_KOKKOS_DETECTOR_H

// class is intended to capture several concepts:
//   1) a dxtbx multipanel detector
//   2) the associated experimental data
//   3) the associated simulated data, in process of being accumulated by kernel calls
//   4) mask data
//   5) possibly other metadata
#include <iostream>

#include "scitbx/array_family/shared.h"
#include "scitbx/array_family/flex_types.h"
#include "dxtbx/model/detector.h"
#include "dxtbx/model/beam.h"
#include "simtbx/nanoBragg/nanoBragg.h"
#include "kokkostbx/kokkos_types.h"
#include "kokkostbx/kokkos_vector3.h"
#include "kokkostbx/kokkos_matrix3.h"

using vec3 = kokkostbx::vector3<CUDAREAL>;
using mat3 = kokkostbx::matrix3<CUDAREAL>;


namespace simtbx { namespace Kokkos {

namespace af = scitbx::af;

struct packed_metrology{
  packed_metrology(){/*printf("NO OPERATION")*/;};
  packed_metrology(dxtbx::model::Detector const &,dxtbx::model::Beam const &);
  packed_metrology(const simtbx::nanoBragg::nanoBragg& nB);
  void show() const;
  af::shared<vec3>sdet;
  af::shared<vec3>fdet;
  af::shared<vec3>odet;
  af::shared<vec3>pix0;
  af::shared<double>dists;
  af::shared<double>Xbeam;
  af::shared<double>Ybeam;
};

struct kokkos_detector{
  inline kokkos_detector(){printf("NO OPERATION, DEVICE NUMBER IS NEEDED");};
  kokkos_detector(int const&, const simtbx::nanoBragg::nanoBragg& nB);
  kokkos_detector(int const&, dxtbx::model::Detector const &, dxtbx::model::Beam const &);
  vector_double_t construct_detail(dxtbx::model::Detector const &);

  inline void show_summary(){
    std::cout << "Detector size: " << m_panel_count << " panel" << ( (m_panel_count>1)? "s" : "" ) << std::endl;
    metrology.show();
  }
  void each_image_allocate();
  void scale_in_place(const double&);
  void write_raw_pixels(simtbx::nanoBragg::nanoBragg&);
  af::flex_double get_raw_pixels();
  void set_active_pixels_on_GPU(af::shared<std::size_t>);
  af::shared<double> get_whitelist_raw_pixels(af::shared<std::size_t>);
  inline void each_image_free(){} //no op in Kokkos
  int h_deviceID;

  //const dxtbx::model::Detector detector;
  packed_metrology const metrology;

  // detector size and dimensions in pixels
  int m_panel_count = -1;
  int m_slow_dim_size = -1;
  int m_fast_dim_size = -1;
  int m_total_pixel_count = -1;

  // kokkos arrays
  vector_double_t m_accumulate_floatimage = vector_double_t("m_accumulate_floatimage", 0);

  // per image input variables, pointers to GPU memory

  // nanoBragg only manages single-panel, multipanel case must mask elsewhere
  vector_ushort_t m_maskimage = vector_ushort_t("m_maskimage", 0);

  // per kernel-call output variables, pointers to GPU memory
  vector_float_t m_omega_reduction = vector_float_t("m_omega_reduction", 0);
  vector_float_t m_max_I_x_reduction = vector_float_t("m_max_I_x_reduction", 0);
  vector_float_t m_max_I_y_reduction = vector_float_t("m_max_I_y_reduction", 0);

  vector_bool_t m_rangemap = vector_bool_t("m_rangemap", 0);

  //most important output, stores the simulated pixel data
  vector_float_t m_floatimage = vector_float_t("m_floatimage", 0);

  // all-panel packed GPU representation of the multi-panel metrology
  view_1d_t<vec3> m_sdet_vector = view_1d_t<vec3>("m_sdet_vector", 0);
  view_1d_t<vec3> m_fdet_vector = view_1d_t<vec3>("m_fdet_vector", 0);
  view_1d_t<vec3> m_odet_vector = view_1d_t<vec3>("m_odet_vector", 0);
  view_1d_t<vec3> m_pix0_vector = view_1d_t<vec3>("m_pix0_vector", 0);
  vector_cudareal_t m_distance = vector_cudareal_t("m_distance", 0);
  vector_cudareal_t m_Xbeam = vector_cudareal_t("m_Xbeam", 0);
  vector_cudareal_t m_Ybeam = vector_cudareal_t("m_Ybeam", 0);

  // all-panel whitelist of active pixels on GPU
  af::shared<std::size_t> active_pixel_list;
  std::size_t m_active_pixel_size;
  vector_size_t m_active_pixel_list = vector_size_t("m_active_pixel_list", 0);
};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_DETECTOR_H

