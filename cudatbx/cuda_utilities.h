#ifndef CUDA_UTILITIES_H
#define CUDA_UTILITIES_H

namespace cudatbx {

  int number_of_gpus();
  void reset_gpu(const int&);

}
#endif // CUDA_UTILITIES_H
