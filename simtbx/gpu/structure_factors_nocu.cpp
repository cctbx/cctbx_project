#include <simtbx/gpu/structure_factors.h>

namespace simtbx {
namespace gpu {

  gpu_energy_channels::gpu_energy_channels(int const& deviceId){
    h_deviceID = deviceId;
  }

  void
  gpu_energy_channels::structure_factors_to_GPU_detail(af::shared<double> linear_amplitudes){
    double * raw_ptr = linear_amplitudes.begin();
    CUDAREAL * cu_Fhkl = NULL;
    d_channel_Fhkl.push_back(cu_Fhkl);
  }

  void gpu_energy_channels::free_detail(){
  }
} // gpu
} // simtbx
