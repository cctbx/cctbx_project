#ifndef CUDATBX_FACTORIALS_CUH
#define CUDATBX_FACTORIALS_CUH

#include <cudatbx/cuda_base.cuh>

namespace cudatbx {
namespace math {

  // f = n!
  template <typename IntType, typename FloatType>
  __device__ FloatType factorial(const IntType& n) {
    FloatType f = 1.0;
    for (IntType i=2; i<(n+1); i++) {
      f *= i;
    }
    return f;
  }

  // f = n!!
  template <typename IntType, typename FloatType>
  __device__ FloatType double_factorial(const IntType& n) {
    if (n < 2) {
      return 1.0;
    }
    else if (n == 2) {
      return 2.0;
    }
    else {
      FloatType f = 1.0;
      for (IntType i=n; i>1; i -= 2) {
        f *= n;
      }
      return f;
    }
  }

}
}
#endif // CUDATBX_FACTORIALS_CUH
