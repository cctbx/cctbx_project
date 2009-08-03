#ifndef SCITBX_RIGID_BODY_MATRIX_HELPERS_H
#define SCITBX_RIGID_BODY_MATRIX_HELPERS_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/tiny_reductions.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace scitbx { namespace rigid_body {

  template <typename FloatType>
  af::tiny<FloatType, 4>
  vec4_normalize(
    af::tiny<FloatType, 4> const& v)
  {
    FloatType den = std::sqrt(af::sum_sq(v));
    SCITBX_ASSERT(den != 0);
    return v / den;
  }

  template <typename FloatType>
  af::tiny<FloatType, 4>
  mat4x3_mul_vec3(
    af::tiny<FloatType, 4*3> const& m,
    vec3<FloatType> const& v)
  {
    af::tiny<FloatType, 4> result;
    matrix::multiply(
      m.begin(),
      v.begin(),
      /*ar*/ 4,
      /*ac*/ 3,
      /*bc*/ 1,
      result.begin());
    return result;
  }

  template <typename FloatType>
  af::tiny<FloatType, 4*3>
  mat4x4_mul_mat4x3(
    af::tiny<FloatType, 4*4> const& a,
    af::tiny<FloatType, 4*3> const& b)
  {
    af::tiny<FloatType, 4*3> result;
    matrix::multiply(
      a.begin(),
      b.begin(),
      /*ar*/ 4,
      /*ac*/ 4,
      /*bc*/ 3,
      result.begin());
    return result;
  }

  template <typename FloatType>
  af::tiny<FloatType, 6>
  mat_6xn_mul_vec_n(
    af::const_ref<FloatType, af::mat_grid> const& a,
    af::const_ref<FloatType> const& b)
  {
    SCITBX_ASSERT(a.accessor().n_rows() == 6);
    unsigned ac = a.accessor().n_columns();
    SCITBX_ASSERT(b.size() == ac);
    af::tiny<FloatType, 6> result;
    matrix::multiply(
      a.begin(),
      b.begin(),
      /*ar*/ 6,
      ac,
      /*bc*/ 1,
      result.begin());
    return result;
  }

  template <typename FloatType>
  af::tiny<FloatType, 6>
  mat_6x6_transpose_mul_vec6(
    af::const_ref<FloatType, af::mat_grid> const& a,
    af::const_ref<FloatType> const& b)
  {
    SCITBX_ASSERT(a.accessor().n_rows() == 6);
    SCITBX_ASSERT(a.accessor().n_columns() == 6);
    SCITBX_ASSERT(b.size() == 6);
    af::tiny<FloatType, 6> result;
    matrix::transpose_multiply(
      a.begin(),
      b.begin(),
      /*ar*/ 6,
      /*ac*/ 6,
      /*bc*/ 1,
      result.begin());
    return result;
  }

  template <typename FloatType>
  af::small<FloatType, 6>
  mat_mxn_mul_vec_n(
    af::const_ref<FloatType, af::mat_grid> const& a,
    af::const_ref<FloatType> const& b)
  {
    unsigned ar = a.accessor().n_rows();
    unsigned ac = a.accessor().n_columns();
    SCITBX_ASSERT(ar <= 6);
    SCITBX_ASSERT(b.size() == ac);
    af::small<FloatType, 6> result(ar);
    matrix::multiply(
      a.begin(),
      b.begin(),
      ar,
      ac,
      /*bc*/ 1,
      result.begin());
    return result;
  }

  template <typename FloatType>
  af::small<FloatType, 6>
  mat_mxn_transpose_mul_vec_n(
    af::const_ref<FloatType, af::mat_grid> const& a,
    af::const_ref<FloatType> const& b)
  {
    unsigned ar = a.accessor().n_rows();
    unsigned ac = a.accessor().n_columns();
    SCITBX_ASSERT(ac <= 6);
    SCITBX_ASSERT(b.size() == ar);
    af::small<FloatType, 6> result(ac);
    matrix::transpose_multiply(
      a.begin(),
      b.begin(),
      ar,
      ac,
      /*bc*/ 1,
      result.begin());
    return result;
  }

  template <typename ElementType, std::size_t N>
  ElementType
  dot_product(
    af::tiny<ElementType, N> const& a,
    af::tiny<ElementType, N> const& b)
  {
    ElementType result = a[0] * b[0];
    for(std::size_t i=1;i<N;i++) result += a[i] * b[i];
    return result;
  }

  template <
    typename FloatType,
    std::size_t ResultSize>
  void
  matrix_mul(
    af::tiny<FloatType, ResultSize>& result,
    af::const_ref<FloatType, af::mat_grid> const& lhs,
    af::const_ref<FloatType> const& rhs)
  {
    SCITBX_ASSERT(ResultSize == lhs.n_rows());
    SCITBX_ASSERT(lhs.n_columns() == rhs.size());
    matrix::multiply(
      lhs.begin(),
      rhs.begin(),
      /*ar*/ lhs.n_rows(),
      /*ac*/ lhs.n_columns(),
      /*bc*/ 1,
      result.begin());
  }

  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  a_transpose_mul_b_mul_a(
    af::const_ref<FloatType, af::mat_grid> const& a,
    af::const_ref<FloatType, af::mat_grid> const& b)
  {
    return af::matrix_multiply(
      af::matrix_transpose_multiply(a, b).const_ref(),
      a);
  }

}} // namespace scitbx::rigid_body

#endif // GUARD
