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
#include <simtbx/nanoBragg/nanoBragg.h>
#include <iostream>

namespace simtbx {
namespace gpu {

namespace af = scitbx::af;

struct gpu_detector{
  gpu_detector(){printf("NO OPERATION, NO DEVICE NUMBER");};
  gpu_detector(int const&, dxtbx::model::Detector const &);

  inline int get_deviceID(){return h_deviceID;}
  inline void show_summary(){
    std::cout << "Detector size" << detector.size() <<std::endl;
    std::cout << "Detector panel 2\n" << detector[2] << std::endl;
  }
  void each_image_allocate_cuda();
  void scale_in_place_cuda(const double&);
  void write_raw_pixels_cuda(simtbx::nanoBragg::nanoBragg&);
  void each_image_free_cuda();

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

  private:
  int _image_size;
};
} // gpu
} // simtbx
#endif // SIMTBX_GPU_DETECTOR_H
