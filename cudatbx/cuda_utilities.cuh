#ifndef CUDA_UTILITIES_CUH
#define CUDA_UTILITIES_CUH

#include <cudatbx/cuda_base.cuh>
#include <cudatbx/cuda_utilities.h>

namespace cudatbx {

  // ==========================================================================

  class cuda_timer {

  public:
    cuda_timer();
    ~cuda_timer();
    void start();
    void stop();
    float get_elapsed_time();

  private:
    cudaEvent_t start_event, stop_event;
    float elapsed_time;
  };

// ==========================================================================

}
#endif // CUDA_UTILITIES_CUH
