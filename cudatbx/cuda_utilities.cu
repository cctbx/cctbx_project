#include <cudatbx/cuda_utilities.h>
#include <cudatbx/cuda_base.cuh>

namespace cudatbx {

  int number_of_gpus() {
    int device_count = 0;
    cudaSafeCall( cudaGetDeviceCount(&device_count) );
    return device_count;
  }

  void reset_gpu(const int& gpu_id) {
    cudaSafeCall( cudaSetDevice(gpu_id) );
    cudaSafeCall( cudaDeviceReset() );
  }

}
