#ifndef SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
#define SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/matrix/lu_decomposition.h>
#include <scitbx/matrix/inversion.h>
#include <scitbx/matrix/diagonal.h>
#include <scitbx/mat_ref.h>
#include <boost/optional.hpp>

namespace scitbx { namespace af {

  template <typename NumType>
  shared<NumType>
  matrix_diagonal(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    shared<NumType> result(a.accessor()[0], init_functor_null<NumType>());
    matrix::diagonal(
      a.begin(), a.accessor()[0], result.begin());
    return result;
  }

  template <typename NumType>
  void
  matrix_diagonal_set_in_place(
    ref<NumType, c_grid<2> > const& a,
    NumType const& value)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typedef typename c_grid<2>::index_value_type ivt;
    ivt n = a.accessor()[0];
    ivt n_sq = n*n;
    ivt n_plus_1 = n + 1;
    for(ivt i=0;i<n_sq;i+=n_plus_1) {
      a[i] = value;
    }
  }

  template <typename NumType>
  void
  matrix_diagonal_add_in_place(
    ref<NumType, c_grid<2> > const& a,
    NumType const& value)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typedef typename c_grid<2>::index_value_type ivt;
    ivt n = a.accessor()[0];
    ivt n_sq = n*n;
    ivt n_plus_1 = n + 1;
    for(ivt i=0;i<n_sq;i+=n_plus_1) {
      a[i] += value;
    }
  }

  template <typename NumType>
  NumType
  matrix_diagonal_sum(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    return matrix::diagonal_sum(a.begin(), a.accessor()[0]);
  }

  template <typename NumType>
  NumType
  matrix_diagonal_product(
    const_ref<NumType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    return matrix::diagonal_product(a.begin(), a.accessor()[0]);
  }

  template <typename NumTypeA, typename NumTypeB>
  versa<
    typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic,
    c_grid<2> >
  matrix_multiply(
    const_ref<NumTypeA, c_grid<2> > const& a,
    const_ref<NumTypeB, c_grid<2> > const& b)
  {
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_ab;
    versa<numtype_ab, c_grid<2> > ab(
      c_grid<2>(a.accessor()[0], b.accessor()[1]),
      init_functor_null<numtype_ab>());
    mat_const_ref<NumTypeA> a_(a.begin(), a.accessor()[0], a.accessor()[1]);
    mat_const_ref<NumTypeB> b_(b.begin(), b.accessor()[0], b.accessor()[1]);
    mat_ref<numtype_ab> ab_(ab.begin(), ab.accessor()[0], ab.accessor()[1]);
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

  template <typename NumType>
  shared<NumType>
  matrix_transpose_multiply_as_packed_u(
    const_ref<NumType, c_grid<2> > const& a)
  {
    unsigned na = a.accessor()[1];
    shared<NumType> ata(na*(na+1)/2, init_functor_null<NumType>());
    matrix::transpose_multiply_as_packed_u(
      a.begin(), a.accessor()[0], na, ata.begin());
    return ata;
  }

  template <typename NumType>
  versa<NumType, c_grid<2> >
  matrix_transpose(const_ref<NumType, c_grid<2> > const& a)
  {
    typedef typename c_grid<2>::value_type index_value_type;
    index_value_type n_rows = a.accessor()[0];
    index_value_type n_columns = a.accessor()[1];
    versa<NumType, c_grid<2> > result(
      c_grid<2>(n_columns, n_rows), init_functor_null<NumType>());
    NumType* r = result.begin();
    for (index_value_type ic=0;ic<n_columns;ic++) {
      std::size_t ir_nc_ic = ic;
      for (index_value_type ir=0;ir<n_rows;ir++,ir_nc_ic+=n_columns) {
        *r++ = a[ir_nc_ic];
      }
    }
    return result;
  }

  template <typename NumType, typename FlexGridIndexType>
  void
  matrix_transpose_in_place(versa<NumType, flex_grid<FlexGridIndexType> >& a)
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
    SCITBX_ASSERT(a.accessor().is_square());
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
    SCITBX_ASSERT(a.accessor().is_square());
    SCITBX_ASSERT(pivot_indices.size() == a.accessor()[0]+1);
    SCITBX_ASSERT(b.size() == a.accessor()[0]);
    shared<FloatType> x(b.begin(), b.end());
    matrix::lu_back_substitution(
      a.begin(), a.accessor()[0], pivot_indices.begin(), x.begin());
    return x;
  }

  template <typename FloatType>
  FloatType
  matrix_determinant_via_lu(
    const_ref<FloatType, c_grid<2> > const& a,
    const_ref<std::size_t> const& pivot_indices)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    SCITBX_ASSERT(pivot_indices.size() == a.accessor()[0]+1);
    FloatType result = matrix_diagonal_product(a);
    if (pivot_indices[a.accessor()[0]] % 2) result = -result;
    return result;
  }

  template <typename FloatType>
  FloatType
  matrix_determinant_via_lu(
    const_ref<FloatType, c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    std::vector<FloatType> a_;
    a_.reserve(a.accessor().size_1d());
    std::copy(a.begin(), a.end(), a_.begin());
    FloatType result;
    try {
      shared<std::size_t>
        pivot_indices = matrix_lu_decomposition_in_place(
          ref<FloatType, c_grid<2> >(&*a_.begin(), a.accessor()));
      result = matrix_diagonal_product(
        const_ref<FloatType, c_grid<2> >(&*a_.begin(), a.accessor()));
      if (pivot_indices[a.accessor()[0]] % 2) result = -result;
    }
    catch (std::runtime_error const& e) {
      if (std::string(e.what())
          != "lu_decomposition_in_place: singular matrix") throw;
      result = 0;
    }
    return result;
  }

  template <typename FloatType>
  void
  matrix_inversion_in_place(
    ref<FloatType, c_grid<2> > const& a,
    ref<FloatType, c_grid<2> > const& b)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    if (   b.accessor()[0] != 0
        && b.accessor()[1] != a.accessor()[0]) {
      throw std::runtime_error(
        "matrix_inversion_in_place: if a is a (n*n) matrix b must be (m*n)");
    }
    matrix::inversion_in_place(
      a.begin(),
      static_cast<std::size_t>(a.accessor()[0]),
      b.begin(),
      static_cast<std::size_t>(b.accessor()[0]));
  }

  template <typename FloatType>
  void
  matrix_inversion_in_place(
    ref<FloatType, c_grid<2> > const& a)
  {
    matrix_inversion_in_place(
      a, ref<FloatType, c_grid<2> >(0, c_grid<2>(0,0)));
  }

  template <typename FloatType>
  boost::optional<FloatType>
  cos_angle(
    const_ref<FloatType> const& a,
    const_ref<FloatType> const& b)
  {
    SCITBX_ASSERT(b.size() == a.size());
    FloatType a_sum_sq = 0;
    FloatType b_sum_sq = 0;
    FloatType a_dot_b = 0;
    for(std::size_t i=0;i<a.size();i++) {
      const FloatType& ai = a[i];
      a_sum_sq += ai * ai;
      const FloatType& bi = b[i];
      b_sum_sq += bi * bi;
      a_dot_b +=  ai * bi;
    }
    if (a_sum_sq == 0 || b_sum_sq == 0) {
      return boost::optional<FloatType>();
    }
    FloatType d = a_sum_sq * b_sum_sq;
    if (d == 0) return boost::optional<FloatType>();
    return boost::optional<FloatType>(a_dot_b / std::sqrt(d));
  }

  template <typename FloatType>
  FloatType
  cos_angle(
    const_ref<FloatType> const& a,
    const_ref<FloatType> const& b,
    FloatType const& value_if_undefined)
  {
    boost::optional<FloatType> result = cos_angle(a, b);
    if (result) return *result;
    return value_if_undefined;
  }

  template <typename FloatType>
  boost::optional<FloatType>
  angle(
    const_ref<FloatType> const& a,
    const_ref<FloatType> const& b)
  {
    boost::optional<FloatType> c = cos_angle(a, b);
    if (!c) return c;
    FloatType result = std::acos(*c);
    return boost::optional<FloatType>(result);
  }

  template <typename FloatType>
  boost::optional<FloatType>
  angle(
    const_ref<FloatType> const& a,
    const_ref<FloatType> const& b,
    bool deg)
  {
    boost::optional<FloatType> rad = angle(a, b);
    if (!rad || !deg) return rad;
    return boost::optional<FloatType>((*rad) / constants::pi_180);
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
