#ifndef SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
#define SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/ref_matrix.h>
#include <scitbx/matrix/lu_decomposition.h>
#include <scitbx/matrix/inversion.h>
#include <scitbx/matrix/diagonal.h>
#include <scitbx/matrix/packed.h>
#include <scitbx/matrix/triangular_systems.h>
#include <scitbx/constants.h>
#include <boost/optional.hpp>
#include <boost/scoped_array.hpp>

namespace scitbx { namespace af {

  template <typename ElementType>
  versa<ElementType, mat_grid>
  versa_mat_grid(
    ElementType const* values,
    unsigned n_rows,
    unsigned n_columns)
  {
    return versa<ElementType, mat_grid>(
      shared_plain<ElementType>(values, values+n_rows*n_columns),
      mat_grid(n_rows, n_columns));
  }

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

  /// The pair (d,f) where d is the diagonal of a and f is the superdiagonal
  template <typename NumType>
  std::pair<shared<NumType>, shared<NumType> >
  matrix_upper_bidiagonal(const_ref<NumType, c_grid<2> > const& a)
  {
    int p = std::min(a.n_rows(), a.n_columns());
    shared<NumType> d(p,   init_functor_null<NumType>()),
                    f(p-1, init_functor_null<NumType>());
    for (int i=0; i<p; ++i) {
      d[i] = a(i, i);
      if (i < p-1) f[i] = a(i, i + 1);
    }
    return std::make_pair(d, f);
  }

  /// The pair (d,f) where d is the diagonal of a and f is the subdiagonal
  template <typename NumType>
  std::pair<shared<NumType>, shared<NumType> >
  matrix_lower_bidiagonal(const_ref<NumType, c_grid<2> > const& a)
  {
    int p = std::min(a.n_rows(), a.n_columns());
    shared<NumType> d(p,   init_functor_null<NumType>()),
                    f(p-1, init_functor_null<NumType>());
    for (int i=0; i<p; ++i) {
      d[i] = a(i, i);
      if (i < p-1) f[i] = a(i+1, i);
    }
    return std::make_pair(d, f);
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

  template <typename NumType, typename NumTypeD>
  void
  matrix_diagonal_set_in_place(
    ref<NumType, c_grid<2> > const& a,
    const_ref<NumTypeD> const& diagonal)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    SCITBX_ASSERT(diagonal.size() == a.accessor()[0]);
    typedef typename c_grid<2>::index_value_type ivt;
    ivt n = a.accessor()[0];
    ivt n_sq = n*n;
    ivt n_plus_1 = n + 1;
    for(ivt i=0,j=0;i<n_sq;i+=n_plus_1) {
      a[i] = diagonal[j++];
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

  //! a * b
  template <typename NumTypeA, typename NumTypeB>
  versa<
    typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic,
    mat_grid>
  matrix_multiply(
    const_ref<NumTypeA, mat_grid> const& a,
    const_ref<NumTypeB, mat_grid> const& b)
  {
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_ab;
    versa<numtype_ab, mat_grid> ab(
      c_grid<2>(a.accessor()[0], b.accessor()[1]),
      init_functor_null<numtype_ab>());
    multiply(a, b, ab.ref());
    return ab;
  }

  //! a.transpose() * b
  template <typename NumTypeA, typename NumTypeB>
  versa<
    typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic,
    mat_grid>
  matrix_transpose_multiply(
    const_ref<NumTypeA, mat_grid> const& a,
    const_ref<NumTypeB, mat_grid> const& b)
  {
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_atb;
    versa<numtype_atb, mat_grid> atb(
      c_grid<2>(a.accessor()[1], b.accessor()[1]),
      init_functor_null<numtype_atb>());
    transpose_multiply(a, b, atb.ref());
    return atb;
  }

  //! a * b.transpose()
  template <typename NumTypeA, typename NumTypeB>
  versa<
    typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic,
    mat_grid>
  matrix_multiply_transpose(
    const_ref<NumTypeA, mat_grid> const& a,
    const_ref<NumTypeB, mat_grid> const& b)
  {
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_atb;
    versa<numtype_atb, mat_grid> atb(
      c_grid<2>(a.accessor()[0], b.accessor()[0]),
      init_functor_null<numtype_atb>());
    multiply_transpose(a, b, atb.ref());
    return atb;
  }

  template <typename NumType>
  shared<NumType>
  matrix_multiply(
    const_ref<NumType, c_grid<2> > const& a,
    const_ref<NumType> const& b)
  {
    shared<NumType> ab(a.accessor()[0], init_functor_null<NumType>());
    af::const_ref<NumType, af::mat_grid> a_(a.begin(),
                                            a.accessor()[0], a.accessor()[1]);
    af::const_ref<NumType, af::mat_grid> b_(b.begin(), b.size(), 1);
    af::ref<NumType, af::mat_grid> ab_(ab.begin(), a.accessor()[0], 1);
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
    af::const_ref<NumType, af::mat_grid> a_(a.begin(), 1, a.size());
    af::const_ref<NumType, af::mat_grid> b_(b.begin(),
                                            b.accessor()[0], b.accessor()[1]);
    af::ref<NumType, af::mat_grid> ab_(ab.begin(), 1, b.accessor()[1]);
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
    af::const_ref<NumType, af::mat_grid> a_(a.begin(), 1, a.size());
    af::const_ref<NumType, af::mat_grid> b_(b.begin(), b.size(), 1);
    af::ref<NumType, af::mat_grid> ab_(&ab, 1, 1);
    multiply(a_, b_, ab_);
    return ab;
  }

  template <typename NumTypeA, typename NumTypeB>
  versa<
    typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic,
    c_grid<2> >
  matrix_multiply_packed_u(
    const_ref<NumTypeA, c_grid<2> > const& a,
    const_ref<NumTypeB> const& b)
  {
    unsigned a_n_rows = a.accessor()[0];
    unsigned a_n_columns = a.accessor()[1];
    SCITBX_ASSERT(dimension_from_packed_size(b.size())
               == a_n_columns);
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_ab;
    versa<numtype_ab, c_grid<2> > ab(
      c_grid<2>(a_n_rows, a_n_columns),
      init_functor_null<numtype_ab>());
    matrix::multiply_packed_u(
      a.begin(), b.begin(), a_n_rows, a_n_columns, ab.begin());
    return ab;
  }

  template <typename NumTypeA, typename NumTypeB>
  shared<typename binary_operator_traits<NumTypeA, NumTypeB>::arithmetic>
  matrix_multiply_packed_u_multiply_lhs_transpose(
    const_ref<NumTypeA, c_grid<2> > const& a,
    const_ref<NumTypeB> const& b)
  {
    unsigned a_n_rows = a.accessor()[0];
    unsigned a_n_columns = a.accessor()[1];
    SCITBX_ASSERT(dimension_from_packed_size(b.size())
               == a_n_columns);
    typedef typename
      binary_operator_traits<NumTypeA, NumTypeB>::arithmetic
        numtype_ab;
    shared<numtype_ab> abat(
      a_n_rows*(a_n_rows+1)/2, init_functor_null<numtype_ab>());
    boost::scoped_array<numtype_ab> ab(new numtype_ab[a_n_rows * a_n_columns]);
    matrix::multiply_packed_u_multiply_lhs_transpose(
      a.begin(),
      b.begin(),
      a_n_rows,
      a_n_columns,
      ab.get(),
      abat.begin());
    return abat;
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
  shared<NumType>
  matrix_transpose_multiply_diagonal_multiply_as_packed_u(
    const_ref<NumType, c_grid<2> > const& a,
    const_ref<NumType> const& diagonal_elements)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    unsigned n = a.accessor()[0];
    shared<NumType> atda(n*(n+1)/2, init_functor_null<NumType>());
    matrix::transpose_multiply_diagonal_multiply_as_packed_u(
      a.begin(), diagonal_elements.begin(), n, atda.begin());
    return atda;
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
    af::ref<NumType, af::mat_grid> a_(a.begin(), n_rows, n_columns);
    a_.transpose_in_place();
    a.resize(flex_grid<FlexGridIndexType>(n_columns, n_rows));
  }

  template <typename NumType>
  versa<NumType, c_grid<2> >
  matrix_rot90(const_ref<NumType, c_grid<2> > const& a, int k)
  {
    typedef typename c_grid<2>::value_type index_value_type;
    index_value_type n_rows = a.accessor()[0];
    index_value_type n_columns = a.accessor()[1];
    versa<NumType, c_grid<2> > result(
      k%2==0?c_grid<2>(n_rows, n_columns):c_grid<2>(n_columns, n_rows),
      init_functor_null<NumType>());
    NumType* r = result.begin();
    std::size_t ir_nc_ic;
    switch (k%4) {
    case 0:
      if (a.begin() != 0) {
        std::copy(a.begin(), a.end(), result.begin());
      }
      break;
    case -3:
    case +1:
      for (index_value_type ic=0;ic<n_columns;ic++) {
        ir_nc_ic = n_columns - 1 - ic;
        for (index_value_type ir=0;ir<n_rows;ir++,ir_nc_ic+=n_columns) {
          *r++ = a[ir_nc_ic];
        }
      }
      break;
    case -2:
    case +2:
      ir_nc_ic = n_columns * n_rows - 1;
      for (index_value_type ic=0;ic<n_columns;ic++) {
        for (index_value_type ir=0;ir<n_rows;ir++,ir_nc_ic--) {
          *r++ = a[ir_nc_ic];
        }
      }
      break;
    case -1:
    case +3:
      for (index_value_type ic=0;ic<n_columns;ic++) {
        ir_nc_ic = n_columns * (n_rows - 1) + ic;
        for (index_value_type ir=0;ir<n_rows;ir++,ir_nc_ic-=n_columns) {
          *r++ = a[ir_nc_ic];
        }
      }
      break;
    }
    return result;
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
  af::shared<FloatType>
  matrix_forward_substitution(const_ref<FloatType> const &l,
                              ref<FloatType> const &b,
                              bool unit_diag=false)
  {
    SCITBX_ASSERT(dimension_from_packed_size(l.size()) == b.size());
    af::shared<FloatType> x(b.begin(), b.end());
    matrix::forward_substitution(b.size(), l.begin(), x.begin(), unit_diag);
    return x;
  }

  template <typename FloatType>
  af::shared<FloatType>
  matrix_back_substitution(const_ref<FloatType> const &u,
                           ref<FloatType> const &b,
                           bool unit_diag=false)
  {
    SCITBX_ASSERT(dimension_from_packed_size(u.size()) == b.size());
    af::shared<FloatType> x(b.begin(), b.end());
    matrix::back_substitution(b.size(), u.begin(), x.begin(), unit_diag);
    return x;
  }

  template <typename FloatType>
  af::shared<FloatType>
  matrix_forward_substitution_given_transpose(const_ref<FloatType> const &u,
                                              ref<FloatType> const &b,
                                              bool unit_diag=false)
  {
    SCITBX_ASSERT(dimension_from_packed_size(u.size()) == b.size());
    af::shared<FloatType> x(b.begin(), b.end());
    matrix::forward_substitution_given_transpose(b.size(), u.begin(), x.begin(),
                                                 unit_diag);
    return x;
  }

  template <typename FloatType>
  af::shared<FloatType>
  matrix_back_substitution_given_transpose(const_ref<FloatType> const &l,
                                           ref<FloatType> const &b,
                                           bool unit_diag=false)
  {
    SCITBX_ASSERT(dimension_from_packed_size(l.size()) == b.size());
    af::shared<FloatType> x(b.begin(), b.end());
    matrix::back_substitution_given_transpose(b.size(), l.begin(), x.begin(),
                                              unit_diag);
    return x;
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
    boost::scoped_array<FloatType> a_(new FloatType[a.accessor().size_1d()]);
    std::copy(a.begin(), a.end(), a_.get());
    FloatType result;
    try {
      shared<std::size_t>
        pivot_indices = matrix_lu_decomposition_in_place(
          ref<FloatType, c_grid<2> >(a_.get(), a.accessor()));
      result = matrix_diagonal_product(
        const_ref<FloatType, c_grid<2> >(a_.get(), a.accessor()));
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
    FloatType cv = *c;
    if      (cv >  1) cv = static_cast<FloatType>(1);
    else if (cv < -1) cv = static_cast<FloatType>(-1);
    FloatType result = std::acos(cv);
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

  template <typename ElementType>
  versa<ElementType, c_grid<2> >
  mat_const_ref_as_versa(
    af::const_ref<ElementType, af::mat_grid> const& m)
  {
    versa<ElementType, c_grid<2> > result(
      c_grid<2>(m.n_rows(), m.n_columns()),
      init_functor_null<ElementType>());
    if (m.begin() != 0) {
      std::copy(m.begin(), m.end(), result.begin());
    }
    else {
      SCITBX_ASSERT(m.size() == 0);
    }
    return result;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_MATRIX_H
