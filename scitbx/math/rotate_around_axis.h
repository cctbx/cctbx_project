#ifndef SCITBX_MATH_ROTATE_AROUND_AXIS_H
#define SCITBX_MATH_ROTATE_AROUND_AXIS_H

#include <scitbx/math/fast_approx_math.h>

namespace scitbx { namespace math {

// Approximate - see code below
template <typename FloatType>
vec3<FloatType>
rotate_point_around_axis(
  vec3<FloatType> const& axis_point_1,
  vec3<FloatType> const& axis_point_2,
  vec3<FloatType> const& point,
  FloatType angle_rad,
  af::const_ref<double> const& sin_table,
  af::const_ref<double> const& cos_table,
  double const& step,
  int const& n)
{
  FloatType xa = axis_point_1[0];
  FloatType ya = axis_point_1[1];
  FloatType za = axis_point_1[2];
  FloatType xb = axis_point_2[0];
  FloatType yb = axis_point_2[1];
  FloatType zb = axis_point_2[2];
  FloatType xl = xb-xa;
  FloatType yl = yb-ya;
  FloatType zl = zb-za;
  FloatType xlsq = xl*xl;
  FloatType ylsq = yl*yl;
  FloatType zlsq = zl*zl;
  FloatType dlsq = xlsq + ylsq + zlsq;
  // Use faster approximate functions instead. Standard functions can be enabled
  // by overload.
  FloatType dl  = std::sqrt(dlsq);
  //FloatType ca  = std::cos(angle_rad);
  //FloatType dsa = std::sin(angle_rad)/dl;
  //FloatType dl  = math::approx_sqrt(dlsq);
  FloatType ca  = math::cos_table(cos_table,angle_rad,step,n,false);
  FloatType dsa = math::sin_table(sin_table,angle_rad,step,n,false)/dl;
  FloatType oca = (1-ca)/dlsq;
  FloatType xlylo = xl*yl*oca;
  FloatType xlzlo = xl*zl*oca;
  FloatType ylzlo = yl*zl*oca;
  FloatType xma = point[0]-xa;
  FloatType yma = point[1]-ya;
  FloatType zma = point[2]-za;
  return vec3<FloatType>(
    xma*(xlsq*oca+ca) +yma*(xlylo-zl*dsa)+zma*(xlzlo+yl*dsa)+xa,
    xma*(xlylo+zl*dsa)+yma*(ylsq*oca+ca) +zma*(ylzlo-xl*dsa)+ya,
    xma*(xlzlo-yl*dsa)+yma*(ylzlo+xl*dsa)+zma*(zlsq*oca+ca)+za);
}

// exact. Copy-paste from the code above. Here we don't need sin/cos tables.
template <typename FloatType>
vec3<FloatType>
rotate_point_around_axis(
  vec3<FloatType> const& axis_point_1,
  vec3<FloatType> const& axis_point_2,
  vec3<FloatType> const& point,
  FloatType angle_rad)
{
  FloatType xa = axis_point_1[0];
  FloatType ya = axis_point_1[1];
  FloatType za = axis_point_1[2];
  FloatType xb = axis_point_2[0];
  FloatType yb = axis_point_2[1];
  FloatType zb = axis_point_2[2];
  FloatType xl = xb-xa;
  FloatType yl = yb-ya;
  FloatType zl = zb-za;
  FloatType xlsq = xl*xl;
  FloatType ylsq = yl*yl;
  FloatType zlsq = zl*zl;
  FloatType dlsq = xlsq + ylsq + zlsq;
  // Use faster approximate functions instead. Standard functions can be enabled
  // by overload.
  FloatType dl  = std::sqrt(dlsq);
  FloatType ca  = std::cos(angle_rad);
  FloatType dsa = std::sin(angle_rad)/dl;
  //FloatType dl  = math::approx_sqrt(dlsq);
  // FloatType ca  = math::cos_table(cos_table,angle_rad,step,n,false);
  // FloatType dsa = math::sin_table(sin_table,angle_rad,step,n,false)/dl;
  FloatType oca = (1-ca)/dlsq;
  FloatType xlylo = xl*yl*oca;
  FloatType xlzlo = xl*zl*oca;
  FloatType ylzlo = yl*zl*oca;
  FloatType xma = point[0]-xa;
  FloatType yma = point[1]-ya;
  FloatType zma = point[2]-za;
  return vec3<FloatType>(
    xma*(xlsq*oca+ca) +yma*(xlylo-zl*dsa)+zma*(xlzlo+yl*dsa)+xa,
    xma*(xlylo+zl*dsa)+yma*(ylsq*oca+ca) +zma*(ylzlo-xl*dsa)+ya,
    xma*(xlzlo-yl*dsa)+yma*(ylzlo+xl*dsa)+zma*(zlsq*oca+ca)+za);
}

template <typename FloatType>
void
rotate_points_around_axis(
  std::size_t const& i_axis_point_1,
  std::size_t const& i_axis_point_2,
  af::ref<vec3<FloatType> > const& all_points,
  af::const_ref<std::size_t> const& rotatable_point_indices,
  FloatType angle_rad,
  af::const_ref<double> const& sin_table,
  af::const_ref<double> const& cos_table,
  double const& step,
  int const& n)
{
  for(std::size_t i=0;i<rotatable_point_indices.size();i++) {
    std::size_t j = rotatable_point_indices[i];
    all_points[j] = rotate_point_around_axis(
      all_points[i_axis_point_1],
      all_points[i_axis_point_2],
      all_points[j],
      angle_rad,
      sin_table,
      cos_table,
      step,
      n);
  }
}

}} // namespace scitbx::math

#endif // SCITBX_MATH_ROTATE_AROUND_AXIS_H
