#ifndef FORM_FACTORS_CUH
#define FORM_FACTORS_CUH

#include <cudatbx/cuda_base.cuh>

namespace cudatbx {
namespace scattering {

  // form factors
  const int max_types = 50;          // number of atom types
  const int max_terms = 10;          // maximum number of terms for form factor
  __device__ __constant__ fType dc_a[max_types * max_terms];
  __device__ __constant__ fType dc_b[max_types * max_terms];
  __device__ __constant__ fType dc_c[max_types];
  __device__ __constant__ int dc_n_types;
  __device__ __constant__ int dc_n_terms;
  __device__ __constant__ bool dc_complex_form_factor;

  template <typename floatType>
  __device__ floatType form_factor(const floatType* a, const floatType* b,
                                   const floatType& c, const int& n_terms,
                                   const floatType& stol_sq) {
    floatType ff = 0.0;
    for (int term=0; term<n_terms; term++) {
      ff += a[term] * __expf(-b[term] * stol_sq);
    }
    ff += c;
    return ff;
  }

  template <typename floatType>
  __device__ floatType form_factor(const int& type, const floatType& stol_sq) {
    return form_factor(&dc_a[type*dc_n_terms],&dc_b[type*dc_n_terms],
                       dc_c[type],dc_n_terms,stol_sq);
  }

}
}
#endif // FORM_FACTORS_CUH
