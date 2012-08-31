#ifndef FORM_FACTORS_CUH
#define FORM_FACTORS_CUH

#include <cudatbx/cuda_base.cuh>

namespace cudatbx {
namespace scattering {

  // form factors
  const int max_types = 50;          // number of atom types
  const int max_terms = 10;          // maximum number of terms for form factor
  __device__ __constant__ fType d_a[max_types * max_terms];
  __device__ __constant__ fType d_b[max_types * max_terms];
  __device__ __constant__ fType d_c[max_types];
  __device__ __constant__ int d_n_types;
  __device__ __constant__ int d_n_terms;

  template <typename floatType>
  __device__ floatType form_factor(const int type, const floatType stol_sq) {
    floatType ff = 0.0;
    for (int term=0; term<d_n_terms; term++) {
      ff += d_a[type*d_n_terms + term] *
        __expf(-d_b[type*d_n_terms + term] * stol_sq);
    }
    ff += d_c[type];

    return ff;
  }
}
}
#endif // FORM_FACTORS_CUH
