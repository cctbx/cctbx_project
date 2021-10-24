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
    SIM(nB),
    sim_use_diffuse(false),diff_gam_a(100), diff_gam_b(100), diff_gam_c(100), diff_sig_a(1), diff_sig_b(1), diff_sig_c(1)
        {
  }
  inline void use_diffuse(af::shared<double>diffuse){
    sim_use_diffuse=true;
    diff_gam_a=diffuse[0]; diff_gam_b=diffuse[1]; diff_gam_c=diffuse[2]; diff_sig_a=diffuse[3]; diff_sig_b=diffuse[4]; diff_sig_c=diffuse[5];
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
  void add_energy_channel_mask_allpanel_cuda(int const&,
    simtbx::gpu::gpu_energy_channels &,
    simtbx::gpu::gpu_detector &,
    af::shared<int> const
  );
  void add_background_cuda(simtbx::gpu::gpu_detector &, int const&);
  void allocate_cuda();
  ~exascale_api();

  const simtbx::nanoBragg::nanoBragg& SIM;
  bool sim_use_diffuse;
  double diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c;
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
