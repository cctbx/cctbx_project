#ifndef SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
#define SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
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
      init_functor_null<FloatType>());
    mat_const_ref<FloatType> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<FloatType> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<FloatType> ab_(ab.begin(), ab.accessor()[0], ab.accessor()[1]);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename FloatType>
  shared<FloatType>
  matrix_multiply(
    const_ref<FloatType, c_grid<2> > const& a,
    const_ref<FloatType> const& b)
  {
    shared<FloatType> ab(a.accessor()[0], init_functor_null<FloatType>());
    mat_const_ref<FloatType> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<FloatType> b_(b.begin(), b.size(), 1);
    mat_ref<FloatType> ab_(ab.begin(), a.accessor()[0], 1);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename FloatType>
  shared<FloatType>
  matrix_multiply(
    const_ref<FloatType> const& a,
    const_ref<FloatType, c_grid<2> > const& b)
  {
    shared<FloatType> ab(b.accessor()[1], init_functor_null<FloatType>());
    mat_const_ref<FloatType> a_(a.begin(), 1, a.size());
    mat_const_ref<FloatType> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<FloatType> ab_(ab.begin(), 1, b.accessor()[1]);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename FloatType>
  FloatType
  matrix_multiply(
    const_ref<FloatType> const& a,
    const_ref<FloatType> const& b)
  {
    FloatType ab;
    mat_const_ref<FloatType> a_(a.begin(), 1, a.size());
    mat_const_ref<FloatType> b_(b.begin(), b.size(), 1);
    mat_ref<FloatType> ab_(&ab, 1, 1);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename FloatType, typename FlexGridIndexType>
  void
  transpose_in_place(versa<FloatType, flex_grid<FlexGridIndexType> >& a)
  {
    SCITBX_ASSERT(a.accessor().nd() == 2);
    SCITBX_ASSERT(a.accessor().is_0_based());
    SCITBX_ASSERT(!a.accessor().is_padded());
    typedef typename FlexGridIndexType::value_type index_value_type;
    index_value_type n_rows = a.accessor().all()[0];
    index_value_type n_columns = a.accessor().all()[1];
    mat_ref<FloatType> a_(a.begin(), n_rows, n_columns);
    a_.transpose_in_place();
    a.resize(flex_grid<FlexGridIndexType>(n_columns, n_rows));
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
