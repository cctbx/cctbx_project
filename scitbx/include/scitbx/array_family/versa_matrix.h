#ifndef SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
#define SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H

#include <scitbx/array_family/versa.h>
#include <scitbx/mat_ref.h>

namespace scitbx { namespace af {

  template <typename FloatType>
  versa<FloatType, c_grid<2> >
  matrix_multiply(
    const_ref<FloatType, c_grid<2> > const& a,
    const_ref<FloatType, c_grid<2> > const& b)
  {
    versa<FloatType, c_grid<2> > ab(
      c_grid<2>(a.accessor()[0], b.accessor()[1]),
      af::init_functor_null<FloatType>());
    mat_const_ref<FloatType> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<FloatType> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<FloatType> ab_(ab.begin(), ab.accessor()[0], ab.accessor()[1]);
    multiply(a_, b_, ab_);
    return ab;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
