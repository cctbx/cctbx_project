#ifndef CUDA_UTILITIES_H
#define CUDA_UTILITIES_H

namespace cudatbx {

  int number_of_gpus();
  void reset_gpu(const int&);
  int calculate_padded_size(const int&, const int&);
  int calculate_blocks_per_grid(const int&, const int&);

}
#endif // CUDA_UTILITIES_H
