/*
 * angle_derivative.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef SCITBX_MATH_ANGLE_DERIVATIVE_H
#define SCITBX_MATH_ANGLE_DERIVATIVE_H

#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/error.h>

namespace scitbx { namespace math {

  using scitbx::vec3;

  /**
   * A function to calculate the derivative of the angle between two vectors
   * u and v with respect to each of the elements of u and v. The result is
   * returned in the form of a pair of vectors.
   */
  inline
  af::tiny<vec3<double>, 2>
  angle_derivative_wrt_vectors(
      vec3<double> u,
      vec3<double> v){

    double u_dot_v = u * v;
    double u_norm = u.length();
    double v_norm = v.length();
    double one_over_u2 = 1./(u_norm * u_norm);
    double one_over_v2 = 1./(v_norm * v_norm);
    double one_over_uv = 1./(u_norm * v_norm);
    double t = std::acos(u_dot_v * one_over_uv);
    double sin_t = std::sin(t);
    SCITBX_ASSERT(sin_t > 0.0);
    double cos_t = std::cos(t);
    double one_over_sin_t = 1./sin_t;
    double cos_t_over_sin_t = cos_t * one_over_sin_t;
    double fac = one_over_sin_t * one_over_uv;

    vec3<double> dtheta_du = cos_t_over_sin_t * one_over_u2 * u - fac * v;
    vec3<double> dtheta_dv = cos_t_over_sin_t * one_over_v2 * v - fac * u;

    return af::tiny<vec3<double>, 2>(dtheta_du, dtheta_dv);
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_ANGLE_DERIVATIVE_H
