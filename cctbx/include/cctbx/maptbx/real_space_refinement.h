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
  int num = sites_cart.size();

  for(int i = 0; i < num; ++i) {
    val -= cctbx::maptbx::non_crystallographic_eight_point_interpolation<double>(
	                               map, gridding_matrix, sites_cart[i], false);
  }

  return val / static_cast<FloatType>(num);
}
}}} // namespace cctbx::maptbx::real_space_refinement

#endif // CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
