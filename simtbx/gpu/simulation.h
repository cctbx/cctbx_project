#ifndef SIMTBX_GPU_SIMULATION_H
#define SIMTBX_GPU_SIMULATION_H

#include <scitbx/array_family/shared.h>
#include <simtbx/nanoBragg/nanoBragg.h>
#include <simtbx/gpu/structure_factors.h>
#include <simtbx/gpu/detector.h>

namespace simtbx {
namespace gpu {

namespace af = scitbx::af;

struct exascale_api {
  inline
  exascale_api(const simtbx::nanoBragg::nanoBragg& nB):
    SIM(nB){
  }

  void show();
  void add_energy_channel_from_gpu_amplitudes_cuda(int const&,
    simtbx::gpu::gpu_energy_channels &,
    simtbx::gpu::gpu_detector &
  );
  void add_energy_channel_mask_allpanel_cuda(int const&,
    simtbx::gpu::gpu_energy_channels &,
    simtbx::gpu::gpu_detector &,
    af::shared<bool>
  );
  void add_background_cuda(simtbx::gpu::gpu_detector &);
  void allocate_cuda();
  ~exascale_api();

  const simtbx::nanoBragg::nanoBragg& SIM;
  CUDAREAL * cu_current_channel_Fhkl;

  CUDAREAL cu_subpixel_size;
  int cu_steps;

  CUDAREAL * cu_beam_vector;
  CUDAREAL * cu_spindle_vector;
  CUDAREAL * cu_source_X, * cu_source_Y, * cu_source_Z, * cu_source_I;
  CUDAREAL * cu_source_lambda, * cu_a0, * cu_b0, * cu_c0;
  CUDAREAL * cu_mosaic_umats;
  CUDAREAL cu_water_size, cu_water_F, cu_water_MW;
  CUDAREAL * cu_polar_vector;
};
} // gpu
} // simtbx
#endif // SIMTBX_GPU_STRUCTURE_FACTORS_H
