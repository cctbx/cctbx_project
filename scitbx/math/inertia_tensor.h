#ifndef SCITBX_MATH_INERTIA_TENSOR_H
#define SCITBX_MATH_INERTIA_TENSOR_H

#include <scitbx/sym_mat3.h>
#include <scitbx/math/accumulators.h>

namespace scitbx { namespace math {

  template <typename FloatType>
  void
  inertia_tensor(
    af::const_ref<vec3<FloatType> > const& points,
    vec3<FloatType> const& pivot,
    sym_mat3<FloatType>& result)
  {
    scitbx::math::accumulator::inertia_accumulator<FloatType>
      accumulate;
    for(std::size_t i_p=0;i_p<points.size();i_p++) {
      accumulate(points[i_p]);
    }
    result = accumulate.inertia_tensor(pivot);
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
    scitbx::math::accumulator::inertia_accumulator<FloatType>
      accumulate;
    for(std::size_t i_p=0;i_p<points.size();i_p++) {
    FloatType w = weights[i_p];
      accumulate(points[i_p], w);
    }
    result = accumulate.inertia_tensor(pivot);
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
