#ifndef SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
#define SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/mat_ref.h>
#include <scitbx/matrix/lu_decomposition.h>
#include <scitbx/matrix/diagonal.h>

namespace scitbx { namespace af {

  template <typename NumType>
  shared<NumType>
  matrix_diagonal(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor()[0] == a.accessor()[1]);
    shared<NumType> result(a.accessor()[0], init_functor_null<NumType>());
    matrix::diagonal(
      a.begin(), a.accessor()[0], result.begin());
    return result;
  }

  template <typename NumType>
  NumType
  matrix_diagonal_sum(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor()[0] == a.accessor()[1]);
    return matrix::diagonal_sum(a.begin(), a.accessor()[0]);
  }

  template <typename NumType>
  NumType
  matrix_diagonal_product(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor()[0] == a.accessor()[1]);
    return matrix::diagonal_product(a.begin(), a.accessor()[0]);
  }

  template <typename NumType>
  versa<NumType, c_grid<2> >
  matrix_multiply(
    const_ref<NumType, c_grid<2> > const& a,
    const_ref<NumType, c_grid<2> > const& b)
  {
    versa<NumType, c_grid<2> > ab(
      c_grid<2>(a.accessor()[0], b.accessor()[1]),
      init_functor_null<NumType>());
    mat_const_ref<NumType> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<NumType> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<NumType> ab_(ab.begin(), ab.accessor()[0], ab.accessor()[1]);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename NumType>
  shared<NumType>
  matrix_multiply(
    const_ref<NumType, c_grid<2> > const& a,
    const_ref<NumType> const& b)
  {
    shared<NumType> ab(a.accessor()[0], init_functor_null<NumType>());
    mat_const_ref<NumType> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<NumType> b_(b.begin(), b.size(), 1);
    mat_ref<NumType> ab_(ab.begin(), a.accessor()[0], 1);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename NumType>
  shared<NumType>
  matrix_multiply(
    const_ref<NumType> const& a,
    const_ref<NumType, c_grid<2> > const& b)
  {
    shared<NumType> ab(b.accessor()[1], init_functor_null<NumType>());
    mat_const_ref<NumType> a_(a.begin(), 1, a.size());
    mat_const_ref<NumType> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<NumType> ab_(ab.begin(), 1, b.accessor()[1]);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename NumType>
  NumType
  matrix_multiply(
    const_ref<NumType> const& a,
    const_ref<NumType> const& b)
  {
    NumType ab;
    mat_const_ref<NumType> a_(a.begin(), 1, a.size());
    mat_const_ref<NumType> b_(b.begin(), b.size(), 1);
    mat_ref<NumType> ab_(&ab, 1, 1);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename NumType, typename FlexGridIndexType>
  void
  transpose_in_place(versa<NumType, flex_grid<FlexGridIndexType> >& a)
  {
    SCITBX_ASSERT(a.accessor().nd() == 2);
    SCITBX_ASSERT(a.accessor().is_0_based());
    SCITBX_ASSERT(!a.accessor().is_padded());
    typedef typename FlexGridIndexType::value_type index_value_type;
    index_value_type n_rows = a.accessor().all()[0];
    index_value_type n_columns = a.accessor().all()[1];
    mat_ref<NumType> a_(a.begin(), n_rows, n_columns);
    a_.transpose_in_place();
    a.resize(flex_grid<FlexGridIndexType>(n_columns, n_rows));
  }

  template <typename FloatType>
  shared<std::size_t>
  matrix_lu_decomposition_in_place(
    ref<FloatType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor()[0] == a.accessor()[1]);
    shared<std::size_t>
      pivot_indices(a.accessor()[0]+1, init_functor_null<std::size_t>());
    matrix::lu_decomposition_in_place(
      a.begin(), a.accessor()[0], pivot_indices.begin());
    return pivot_indices;
  }

  template <typename FloatType>
  shared<FloatType>
  matrix_lu_back_substitution(
    const_ref<FloatType, c_grid<2> > const& a,
    const_ref<std::size_t> const& pivot_indices,
    const_ref<FloatType> const& b)
  {
    SCITBX_ASSERT(a.accessor()[0] == a.accessor()[1]);
    SCITBX_ASSERT(pivot_indices.size() == a.accessor()[0]+1);
    SCITBX_ASSERT(b.size() == a.accessor()[0]);
    shared<FloatType> x(b.begin(), b.end());
    matrix::lu_back_substitution(
      a.begin(), a.accessor()[0], pivot_indices.begin(), x.begin());
    return x;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
