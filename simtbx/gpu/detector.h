#ifndef SIMTBX_GPU_DETECTOR_H
#define SIMTBX_GPU_DETECTOR_H
// the *.h file can be included for all compiles regardless of the availability of CUDA

/* class is intended to capture several concepts:
   1) a dxtbx multipanel detector
   2) the associated experimental data
   3) the associated simulated data, in process of being accumulated by kernel calls
   4) mask data
   5) possibly other metadata
*/

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/beam.h>
#include <simtbx/nanoBragg/nanoBragg.h>
#include <iostream>

namespace simtbx {
namespace gpu {

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

struct gpu_detector{
  inline gpu_detector(){printf("NO OPERATION, DEVICE NUMBER IS NEEDED");};
  gpu_detector(int const&, const simtbx::nanoBragg::nanoBragg& nB);
  gpu_detector(int const&, dxtbx::model::Detector const &, dxtbx::model::Beam const &);
  void construct_detail(int const&, dxtbx::model::Detector const &);

  inline int get_deviceID(){return h_deviceID;}
  inline void show_summary(){
    std::cout << "Detector size" << cu_n_panels <<std::endl;
    metrology.show();
  }
  void each_image_allocate();
  void scale_in_place(const double&);
  void write_raw_pixels(simtbx::nanoBragg::nanoBragg&);
  af::flex_double get_raw_pixels();
  void set_active_pixels_on_GPU(af::shared<std::size_t>);
  af::shared<double> get_whitelist_raw_pixels(af::shared<std::size_t>);
  void each_image_free();

  void free_detail();
  inline ~gpu_detector(){ {free_detail();} }

  int h_deviceID;
  const dxtbx::model::Detector detector;

  int  cu_n_panels, cu_slow_pixels, cu_fast_pixels; /* variables on host only */
  double * cu_accumulate_floatimage; /* pointer to GPU memory */

  /* per image input variables, pointers to GPU memory */
  int unsigned short * cu_maskimage; // nanoBragg only manages single-panel, multipanel case must mask elsewhere

  /* per kernel-call output variables, pointers to GPU memory */
  float * cu_omega_reduction, * cu_max_I_x_reduction, * cu_max_I_y_reduction;
  bool * cu_rangemap;
  float * cu_floatimage; //most important output, stores the simulated pixel data

  /* all-panel packed GPU representation of the multi-panel metrology */
  CUDAREAL * cu_sdet_vector, * cu_fdet_vector;
  CUDAREAL * cu_odet_vector, * cu_pix0_vector;
  CUDAREAL * cu_distance, * cu_Xbeam, * cu_Ybeam;

  /* all-panel whitelist of active pixels on GPU */
  af::shared<std::size_t> active_pixel_list;
  std::size_t * cu_active_pixel_list;

  packed_metrology const metrology;
  private:
  int _image_size;
};
} // gpu
} // simtbx
#endif // SIMTBX_GPU_DETECTOR_H
