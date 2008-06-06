#include "sampled_model_density.h"
#include <scitbx/math/basic_statistics.h>

int main() {
  sampled_model_density_test_case test;
  std::cout
    << "sampled model density dynamically loaded by main"
    << std::endl;
  scitbx::math::basic_statistics<double> stats(
    test.sampled_model_density->real_map().as_1d().const_ref());
  std::cout << "map points: " << stats.n << std::endl;
  std::cout << "    min:    " << stats.min << std::endl;
  std::cout << "    max:    " << stats.max << std::endl;
  std::cout << "    mean:   " << stats.mean << std::endl;
}
