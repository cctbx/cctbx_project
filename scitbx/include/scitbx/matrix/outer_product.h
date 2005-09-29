#ifndef SCITBX_MATRIX_OUTER_PRODUCT_H
#define SCITBX_MATRIX_OUTER_PRODUCT_H

#include <scitbx/array_family/versa.h>

namespace scitbx { namespace matrix {

  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> >
  outer_product(
    af::const_ref<FloatType> const& lhs,
    af::const_ref<FloatType> const& rhs)
  {
    af::versa<FloatType, af::c_grid<2> >
      result(
        af::c_grid<2>(lhs.size(), rhs.size()),
        af::init_functor_null<FloatType>());
    FloatType* r = result.begin();
    for(unsigned i=0;i<lhs.size();i++) {
      FloatType li = lhs[i];
      for(unsigned j=0;j<rhs.size();j++) {
        *r++ = li * rhs[j];
      }
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_OUTER_PRODUCT_H
