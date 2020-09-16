//#include <cudatbx/cuda_base.cuh>
#include "hip/hip_runtime.h"
#include <simtbx/gpu/structure_factors.h>

namespace simtbx {
namespace gpu {

  gpu_energy_channels::gpu_energy_channels(int const& deviceId){
    h_deviceID = deviceId;
    hipSetDevice(deviceId);
  }

  void 
  gpu_energy_channels::structure_factors_to_GPU_detail(af::shared<double> linear_amplitudes){
    double * raw_ptr = linear_amplitudes.begin();
    CUDAREAL * cu_Fhkl = NULL;
    hipMalloc((void ** )&cu_Fhkl, 
                     sizeof(*cu_Fhkl) * linear_amplitudes.size());
    hipMemcpy(cu_Fhkl, raw_ptr, 
                     sizeof(*cu_Fhkl) * linear_amplitudes.size(), hipMemcpyHostToDevice);

    d_channel_Fhkl.push_back(cu_Fhkl);
  }

  void gpu_energy_channels::free_detail(){
        hipSetDevice(h_deviceID);
        for (int i_cu_ptr=0; i_cu_ptr < d_channel_Fhkl.size(); ++i_cu_ptr){
          hipFree(d_channel_Fhkl[i_cu_ptr]);
        }
  }
} // gpu
} // simtbx
