#include "sampled_model_density.h"

int main() {
  sampled_model_density_test_case test;
  std::cout
    << "sampled model density dynamically loaded by main"
    << std::endl;
  std::cout << test.sampled_model_density->real_map().as_1d().const_ref()
            << std::endl;

}
