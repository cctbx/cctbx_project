#ifndef CCTBX_MAPTBX_MASK_UTILS_H
#define CCTBX_MAPTBX_MASK_UTILS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <cstddef>
#include <cctbx/uctbx.h>

namespace cctbx { namespace maptbx {

af::shared<scitbx::vec3<double> >
sample_mask_regions(
  af::const_ref<int, af::flex_grid<> > const& mask,
  int n_zone,
  int volume,
  int sampling_rate,
  cctbx::uctbx::unit_cell const& unit_cell)

{
  CCTBX_ASSERT(mask.accessor().nd() == 3);
  CCTBX_ASSERT(mask.accessor().all().all_gt(0));
  af::shared <  scitbx::vec3<double> > result_cart;
  af::c_grid<3> a = mask.accessor();
  int n_sampled = 0;
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
        if (mask(i,j,k) == n_zone) {
          n_sampled += 1;
          if (n_sampled == 1 ||
              n_sampled == volume ||
              n_sampled % sampling_rate == 0) {
            cctbx::fractional<> grid_frac = cctbx::fractional<>(
                double(i)/a[0],
                double(j)/a[1],
                double(k)/a[2]);
            cctbx::cartesian<> site_cart = unit_cell.orthogonalize(grid_frac);
            // result_cart[mask(i,j,k)].push_back(site_cart);
            result_cart.push_back(site_cart);
            // result_cart.push_back(scitbx::vec3<double>(i,j,k));
          }
        }
      }
    }
  }
  return result_cart;

}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_MASK_UTILS_H
