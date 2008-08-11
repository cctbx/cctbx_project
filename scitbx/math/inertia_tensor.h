#ifndef SCITBX_MATH_INERTIA_TENSOR_H
#define SCITBX_MATH_INERTIA_TENSOR_H

#include <scitbx/sym_mat3.h>

namespace scitbx { namespace math {

  template <typename FloatType>
  void
  inertia_tensor(
    af::const_ref<vec3<FloatType> > const& points,
    vec3<FloatType> const& pivot,
    sym_mat3<FloatType>& result)
  {
    for(std::size_t i=0;i<points.size();i++) {
      vec3<FloatType> p = points[i] - pivot;
      vec3<FloatType> pp(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
      result(0,0) += pp[1] + pp[2];
      result(1,1) += pp[0] + pp[2];
      result(2,2) += pp[0] + pp[1];
      result(0,1) -= p[0] * p[1];
      result(0,2) -= p[0] * p[2];
      result(1,2) -= p[1] * p[2];
    }
  }

  template <typename FloatType>
  void
  inertia_tensor(
    af::const_ref<vec3<FloatType> > const& points,
    af::const_ref<FloatType> const& weights,
    vec3<FloatType> const& pivot,
    sym_mat3<FloatType>& result)
  {
    SCITBX_ASSERT(weights.size() == points.size());
    for(std::size_t i=0;i<points.size();i++) {
      vec3<FloatType> p = points[i] - pivot;
      vec3<FloatType> pp(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
      FloatType w = weights[i];
      result(0,0) += w * (pp[1] + pp[2]);
      result(1,1) += w * (pp[0] + pp[2]);
      result(2,2) += w * (pp[0] + pp[1]);
      result(0,1) -= w * p[0] * p[1];
      result(0,2) -= w * p[0] * p[2];
      result(1,2) -= w * p[1] * p[2];
    }
  }

  template <typename FloatType>
  sym_mat3<FloatType>
  inertia_tensor(
    af::const_ref<vec3<FloatType> > const& points,
    vec3<FloatType> const& pivot)
  {
    sym_mat3<FloatType> result(0,0,0,0,0,0);
    inertia_tensor(points, pivot, result);
    return result;
  }

  template <typename FloatType>
  sym_mat3<FloatType>
  inertia_tensor(
    af::const_ref<vec3<FloatType> > const& points,
    af::const_ref<FloatType> const& weights,
    vec3<FloatType> const& pivot)
  {
    sym_mat3<FloatType> result(0,0,0,0,0,0);
    inertia_tensor(points, weights, pivot, result);
    return result;
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_INERTIA_TENSOR_H
