#ifndef SIMTBX_KOKKOS_DETECTOR_H
#define SIMTBX_KOKKOS_DETECTOR_H
// the *.h file can be included for all compiles regardless of the availability of CUDA

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

#include <Kokkos_Core.hpp>
#include <iostream>

using vector_view_t = Kokkos::View<double*>;
using host_view_t = vector_view_t::HostMirror;

namespace simtbx {
namespace Kokkos {

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
  //kokkos_detector(const simtbx::nanoBragg::nanoBragg& nB);
  kokkos_detector(dxtbx::model::Detector const &, dxtbx::model::Beam const &);
  vector_view_t construct_detail(dxtbx::model::Detector const &);

  //inline void show_summary(){
  //  std::cout << "Detector size" << m_panel_count <<std::endl;
  //  metrology.show();
  //}
  //void each_image_allocate_cuda();
  void scale_in_place_cuda(const double&);
  void write_raw_pixels_cuda(simtbx::nanoBragg::nanoBragg&);
  //af::flex_double get_raw_pixels_cuda();
  //void set_active_pixels_on_KOKKOS(af::shared<int>);
  //af::shared<double> get_whitelist_raw_pixels_cuda(af::shared<std::size_t>);
  //void each_image_free_cuda();

  //void free_detail();
  //inline ~kokkos_detector(){ {free_detail();} }

  const dxtbx::model::Detector detector;

  int m_panel_count = -1;
  int m_slow_dim_size = -1;
  int m_fast_dim_size = -1;
  int m_total_pixel_count = -1;
  vector_view_t m_accumulate_floatimage;
  double * cu_accumulate_floatimage; // pointer to GPU memory

  // per image input variables, pointers to GPU memory
  int unsigned short * cu_maskimage; // nanoBragg only manages single-panel, multipanel case must mask elsewhere

  // per kernel-call output variables, pointers to GPU memory
  float * cu_omega_reduction, * cu_max_I_x_reduction, * cu_max_I_y_reduction;
  bool * cu_rangemap;
  float * cu_floatimage; //most important output, stores the simulated pixel data

  // all-panel packed GPU representation of the multi-panel metrology
  CUDAREAL * cu_sdet_vector, * cu_fdet_vector;
  CUDAREAL * cu_odet_vector, * cu_pix0_vector;
  CUDAREAL * cu_distance, * cu_Xbeam, * cu_Ybeam;

  // all-panel whitelist of active pixels on GPU
  af::shared<int> active_pixel_list;
  int * cu_active_pixel_list;

  packed_metrology const metrology;
};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_DETECTOR_H

