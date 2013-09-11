#include <cudatbx/cuda_utilities.cuh>

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

  int calculate_padded_size(const int& size, const int& padding) {
    int padded_size = int(std::floor(size/padding + 1.0)) * padding;
    return padded_size;
  }

  int calculate_blocks_per_grid(const int& size, const int& threads_per_block) {
    int blocks_per_grid = (size + threads_per_block - 1)/threads_per_block;
    return blocks_per_grid;
  }

  /* ==========================================================================
     Basic timer for CUDA using events, one stream only, no checks

     Usage:

       cuda_timer t;
       t.start();

       < run CUDA stuff >

       t.stop();
       float elapsed_time = t.get_elapsed_time();
       std::cout << elapsed_time << "\n";
  */
  cudatbx::cuda_timer::cuda_timer() {
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
  }

  cudatbx::cuda_timer::~cuda_timer() {
    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);
  }

  void cudatbx::cuda_timer::start() {
    cudaEventRecord(start_event);
  }

  void cudatbx::cuda_timer::stop() {
    cudaEventRecord(stop_event);
    cudaEventSynchronize(stop_event);
    cudaEventElapsedTime(&elapsed_time, start_event, stop_event);
  }

  float cudatbx::cuda_timer::get_elapsed_time() {
    return elapsed_time;
  }

  // ==========================================================================

}
