#ifndef CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
#define CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H

#include <cctbx/maptbx/eight_point_interpolation.h>

namespace cctbx { namespace maptbx { namespace real_space_refinement {
template <typename FloatType>
FloatType
residual(
  af::const_ref<FloatType, af::flex_grid<> > const& map,
  scitbx::mat3<FloatType> const& gridding_matrix,
  af::const_ref<scitbx::vec3<FloatType> > const& sites_cart)
{
  FloatType val = 0.0;
  std::size_t num = sites_cart.size();

  for(std::size_t i = 0; i < num; ++i) {
    val -= cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, sites_cart[i], false);
  }
  return val / static_cast<FloatType>(num);
}

template <typename FloatType>
af::shared<scitbx::vec3<FloatType> >
gradients(
  af::const_ref<FloatType, af::flex_grid<> > const& map,
  scitbx::mat3<FloatType> const& gridding_matrix,
  af::const_ref<scitbx::vec3<FloatType> > const& sites_cart)
{
  af::shared<scitbx::vec3<FloatType> > grad_vals;
  scitbx::vec3<FloatType> grad_val;
  FloatType delta_x, delta_y, delta_z;
  delta_x = delta_y = delta_z = 1.0;
  std::size_t num = sites_cart.size();

  for(std::size_t i = 0; i < num; ++i) {
    scitbx::vec3<FloatType> temp;
    temp = sites_cart[i];
    temp[0] = sites_cart[i][0] + delta_x;
    grad_val[0] = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    temp[0] = sites_cart[i][0] - delta_x;
    grad_val[0] -= cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    grad_val[0] /= 2.0;
    temp = sites_cart[i];
    temp[1] = sites_cart[i][1] + delta_y;
    grad_val[1] = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    temp[1] = sites_cart[i][1] - delta_y;
    grad_val[1] -= cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    grad_val[1] /= 2.0;
    temp = sites_cart[i];
    temp[2] = sites_cart[i][2] + delta_z;
    grad_val[2] = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    temp[2] = sites_cart[i][2] - delta_z;
    grad_val[2] -= cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, temp, false);
    grad_val[2] /= 2.0;
    grad_vals.push_back(grad_val);
  }

  return grad_vals;
}

}}} // namespace cctbx::maptbx::real_space_refinement

#endif // CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
