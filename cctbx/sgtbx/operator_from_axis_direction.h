#ifndef CCTBX_SGTBX_OPERATOR_FROM_AXIS_DIRECTION_H
#define CCTBX_SGTBX_OPERATOR_FROM_AXIS_DIRECTION_H

#include <cctbx/uctbx.h>
#include <scitbx/constants.h>

namespace cctbx { namespace sgtbx {

  inline
  uc_mat3
  two_fold_operator_from_axis_direction(uc_vec3 const& ev_cart)
  {
    double f = 2. / ev_cart.length_sq();
    double x = ev_cart[0];
    double y = ev_cart[1];
    double z = ev_cart[2];
    return uc_mat3(f*x*x-1., f*x*y,    f*x*z,
                   f*y*x,    f*y*y-1., f*y*z,
                   f*z*x,    f*z*y,    f*z*z-1.);
  }

  inline
  uc_mat3
  n_fold_operator_from_axis_direction(
    uc_vec3 const& ev_cart,
    int n,
    int sense=1)
  {
    if (n == 1) return uc_mat3(1,0,0,0,1,0,0,0,1);
    if (n == 2) return two_fold_operator_from_axis_direction(ev_cart);
    CCTBX_ASSERT(sense == 1 || sense == -1);
    CCTBX_ASSERT(n == 1 || n == 2 || n == 3 || n == 4 || n == 6);
    uc_vec3 v = ev_cart.normalize();
    double angle = scitbx::constants::two_pi / (sense * n);
    double c = std::cos(angle);
    double s = std::sin(angle);
    double d = 1 - c;
    return uc_mat3(
      c+d*v[0]*v[0],      d*v[0]*v[1]-s*v[2], d*v[0]*v[2]+s*v[1],
      d*v[0]*v[1]+s*v[2], c+d*v[1]*v[1],      d*v[1]*v[2]-s*v[0],
      d*v[0]*v[2]-s*v[1], d*v[1]*v[2]+s*v[0], c+d*v[2]*v[2]);
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_OPERATOR_FROM_AXIS_DIRECTION_H
