#ifndef CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
#define CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
// Done by Erik McKee

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
  scitbx::vec3<FloatType> deltas;
  std::size_t num = sites_cart.size();
  for(std::size_t i = 0; i < 3; ++i) {
    deltas[i] = 1.0;
  }

  for(std::size_t i = 0; i < num; ++i) {
    for(std::size_t j = 0; j < 3; ++j) {
      scitbx::vec3<FloatType> temp;
      temp = sites_cart[i];
      temp[j] = sites_cart[i][j] - deltas[j];
      grad_val[j] = cctbx::maptbx::
        non_crystallographic_eight_point_interpolation<FloatType>(
        map, gridding_matrix, temp, false);
      temp[j] = sites_cart[i][j] + deltas[j];
      grad_val[j] -= cctbx::maptbx::
        non_crystallographic_eight_point_interpolation<FloatType>(
        map, gridding_matrix, temp, false);
      grad_val[j] /= 2.0;
    }
    grad_vals.push_back(grad_val);
  }

  return grad_vals;
}

}}} // namespace cctbx::maptbx::real_space_refinement

#endif // CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
