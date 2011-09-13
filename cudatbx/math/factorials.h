#ifndef CUDATBX_FACTORIALS_H
#define CUDATBX_FACTORIALS_H

#include <cudatbx/cuda_base.h>

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
    FloatType f = 1.0;
    for (IntType i=n; i>2; i -= 2) {
      f *= n;
    }
    return f;
  }

}
}
#endif // CUDATBX_FACTORIALS_H
