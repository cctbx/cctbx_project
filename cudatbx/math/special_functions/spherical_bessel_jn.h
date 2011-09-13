#ifndef CUDATBX_SPHERICAL_BESSEL_JN_H
#define CUDATBX_SPHERICAL_BESSEL_JN_H

#include <scitbx/array_family/shared.h>

#include <cudatbx/cuda_base.h>
#include <cudatbx/math/factorials.h>

namespace cudatbx {
namespace math {
namespace special_functions {

  template <typename IntType, typename FloatType>
  __device__ FloatType a_k(const IntType& k, const IntType& n) {
    if (k > n) {
      return 0.0;
    }
    else {
      FloatType f_npk_k = 1.0;  // (n + k)! / k! = (k + 1) * (k + 2) ... (k + n)
      for (IntType i=1.0; i<(n+1); i++) {
        f_npk_k *= k + i;
      }
      return f_npk_k/(pow(2.0,k)*factorial<IntType,FloatType>(n-k));
    }
  }

  template <typename IntType, typename FloatType>
  __device__ FloatType spherical_bessel_jn
    (const IntType& n, const FloatType& z) {
    // handle z ~ 0.0
    if (z < 1.0e-10) {
      if (n == 0) {
        return 1.0;
      }
      return 0.0;
    }
    // handle low z region
    // interface will become discontinuous above order 50
    else if (z < 0.014*n*n) {
      return pow(z,n) / double_factorial<IntType,FloatType>(2*n + 1);
    }
    // all other values for z
    else {
      FloatType sin_sum, cos_sum, sin_z, cos_z;
      // sin and cos terms
      sincos(z - CUDART_PIO2*n,&sin_z,&cos_z);

      // sums
      sin_sum = 0.0;
      for (IntType k=0; k<IntType(floor(n/2.0))+1; k++) {
        sin_sum += pow(-1.0,k) * a_k<IntType,FloatType>(2*k,n) / pow(z,2*k+1);
      }
      cos_sum = 0.0;
      for (IntType k=0; k<IntType(floor((n-1)/2.0))+1; k++) {
        cos_sum += pow(-1.0,k) * a_k<IntType,FloatType>(2*k+1,n) / pow(z,2*k+2);
      }
      return sin_z*sin_sum + cos_z*cos_sum;
    }
  }

  __global__ void spherical_bessel_jn_kernel
    (const int, const double*, const int, double*);
  __device__ double2 spherical_bessel_j0_j1(double);
  scitbx::af::shared<double> cuda_spherical_bessel_jn
    (const int&, const scitbx::af::const_ref<double>&, const int&);

}
}
}
#endif // CUDATBX_SPHERICAL_BESSEL_JN_H
