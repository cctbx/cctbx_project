#ifndef SCITBX_RIGID_BODY_MATRIX_HELPERS_H
#define SCITBX_RIGID_BODY_MATRIX_HELPERS_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/tiny_reductions.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/matrix/multiply.h>

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

  template <typename ElementType, std::size_t N>
  ElementType
  dot_product(
    af::tiny<ElementType, N> const& a,
    af::tiny<ElementType, N> const& b)
  {
    ElementType result = a[0] * a[0];
    for(std::size_t i=1;i<N;i++) result += a[i] * a[i];
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

}} // namespace scitbx::rigid_body

#endif // GUARD
