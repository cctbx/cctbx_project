#ifndef CCTBX_MAPTBX_MASK_UTILS_H
#define CCTBX_MAPTBX_MASK_UTILS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <cstddef>
#include <cctbx/uctbx.h>

namespace cctbx { namespace maptbx {

class sample_all_mask_regions {
protected:
  af::shared<af::shared<scitbx::vec3<double> > > result_cart_;

public:
  sample_all_mask_regions(
    af::const_ref<int, af::flex_grid<> > const& mask,
    af::shared<int> const& volumes,
    af::shared<int> const& sampling_rates,
    cctbx::uctbx::unit_cell const& unit_cell)
  {
    CCTBX_ASSERT(mask.accessor().nd() == 3);
    CCTBX_ASSERT(mask.accessor().all().all_gt(0));
    CCTBX_ASSERT(volumes.size() == sampling_rates.size());

    // initializing result_cart_ array
    for (int i=0; i<volumes.size(); i++)
    {
      af::shared<scitbx::vec3<double> > tmp;
      result_cart_.push_back(tmp);
    }
    af::shared<int> count_dict(volumes.size(), 0);

    af::c_grid<3> a = mask.accessor();

    for(int i = 0; i < a[0]; i++) {
      for(int j = 0; j < a[1]; j++) {
        for(int k = 0; k < a[2]; k++) {
          int mv = mask(i,j,k);
          if (mv > 0) {
            count_dict[mv] += 1;
            if (count_dict[mv] == 1 ||
                count_dict[mv] == volumes[mv] ||
                count_dict[mv] % sampling_rates[mv] == 0) {
              cctbx::fractional<> grid_frac = cctbx::fractional<>(
                  double(i)/a[0],
                  double(j)/a[1],
                  double(k)/a[2]);
              cctbx::cartesian<> site_cart = unit_cell.orthogonalize(grid_frac);
              result_cart_[mv].push_back(site_cart);
            }
          }

        } //for k
      } // for j
    } // for i
  } // constructor

  af::shared<scitbx::vec3<double> >
  get_array(int n)
  {
    CCTBX_ASSERT(n < result_cart_.size());
    return result_cart_[n];
  }

}; // class sample_all_mask_regions

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_MASK_UTILS_H
