#ifndef SCITBX_MATRIX_CHOLESKY_H
#define SCITBX_MATRIX_CHOLESKY_H

#include <scitbx/matrix/packed.h>
#include <scitbx/array_family/shared.h>
#include <cmath>

namespace scitbx { namespace matrix { namespace cholesky {

  template <typename FloatType>
  af::shared<FloatType>
  decomposition(
    af::const_ref<FloatType> const& a,
    FloatType const& relative_eps=1.e-15)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    FloatType eps;
    if (relative_eps <= 0 || a.size() == 0) {
      eps = 0;
    }
    else {
      eps = relative_eps * af::max_absolute(a);
    }
    af::shared<FloatType> result(a.size(), af::init_functor_null<FloatType>());
    FloatType *c = result.begin();
    std::size_t i_a = 0;
    std::size_t k0 = 0;
    for(unsigned k=0;k<n;k++,k0+=k) {
      std::size_t kj = k0;
      FloatType sum = 0;
      for(unsigned j=0;j<k;j++) {
        FloatType const& c_kj = c[kj++];
        sum += c_kj * c_kj;
      }
      FloatType d = a[i_a++] - sum;
      if (d <= eps) return af::shared<FloatType>();
      d = std::sqrt(d);
      c[kj++] = d;
      std::size_t i0 = kj;
      for(unsigned i=k+1;i<n;i++,i0+=i) {
        std::size_t ij = i0;
        kj = k0;
        sum = 0;
        for(unsigned j=0;j<k;j++) {
          sum += c[ij++] * c[kj++];
        }
        c[ij] = (a[i_a++] - sum) / d;
      }
    }
    return result;
  }

}}} // namespace scitbx::matrix::cholesky

#endif // SCITBX_MATRIX_CHOLESKY_H
