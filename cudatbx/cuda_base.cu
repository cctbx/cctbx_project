#include <cudatbx/cuda_base.h>

namespace cudatbx {

  int number_of_gpus() {
    int device_count = 0;
    cudaSafeCall( cudaGetDeviceCount(&device_count) );
    return device_count;
  }

  void reset_gpu(const int& gpu_id) {
    int device_count = number_of_gpus();
    if (gpu_id < device_count) {
      cudaSafeCall( cudaSetDevice(gpu_id) );
      cudaSafeCall( cudaDeviceReset() );
    }
    else {
      std::cerr << "WARNING: Device " << gpu_id << " does not exist.\n";
    }
  }

}
