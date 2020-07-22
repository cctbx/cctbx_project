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

class binary_filter {

private:
  af::versa<double, af::c_grid<3> > map_new;
  af::tiny<int, 3> map_dimensions;

public:
  binary_filter (
    af::const_ref<double, af::flex_grid<> > const& map,
    float const& threshold)
  {
    CCTBX_ASSERT(map.accessor().nd() == 3);
    CCTBX_ASSERT(map.accessor().all().all_gt(0));

    af::c_grid<3> a = map.accessor();

    map_dimensions = af::adapt(map.accessor().all());
    map_new.resize(af::c_grid<3>(map_dimensions), 0);

    int boundary = 1;
    float value = 0;
    float cutoff = 27. * threshold;  // sum of 27 values must be > cutoff
    // Create new map where value at each point is 1 if
    // average of all 27 points next to this point plus this point > threshold
    int i_min=(boundary);
    int i_max=a[0]-boundary;
    int j_min=(boundary);
    int j_max=a[1]-boundary;
    int k_min=(boundary);
    int k_max=a[2]-boundary;

    for(int i = i_min ; i < i_max; i++) {
      for(int j = j_min ; j < j_max; j++) {
        for(int k = k_min ; k < k_max; k++) {
              value = 0.;
              for ( int ii = -1; ii < 2; ii++) {
                for ( int jj = -1; jj < 2; jj++) {
                  for ( int kk = -1; kk < 2; kk++) {
                    value += map(i+ii,j+jj,k+kk);
                  }
                }
              }
              if (value >= cutoff) {
                map_new(i,j,k) = 1.;
              } else {
                map_new(i,j,k) = 0.;
              }
        } //for k
      } // for j
    } // for i
  } // constructor
  af::versa<double, af::c_grid<3> > result() {return map_new;}

}; // class binary_filter

class zero_boundary_box_map {

private:
  af::versa<double, af::c_grid<3> > map_new;
  af::tiny<int, 3> map_dimensions;

public:
  zero_boundary_box_map(
    af::const_ref<double, af::flex_grid<> > const& mask,
    int const& boundary)
  {
    CCTBX_ASSERT(mask.accessor().nd() == 3);
    CCTBX_ASSERT(mask.accessor().all().all_gt(0));

    af::c_grid<3> a = mask.accessor();
    CCTBX_ASSERT(boundary >= 0);  // if 0, all 1 in map
    CCTBX_ASSERT(2*boundary < a[0]);
    CCTBX_ASSERT(2*boundary < a[1]);
    CCTBX_ASSERT(2*boundary < a[2]);

    map_dimensions = af::adapt(mask.accessor().all());
    map_new.resize(af::c_grid<3>(map_dimensions), 0);

    // set map=1  except within boundary of any side
    int i_min=(boundary);
    int i_max=a[0]-boundary;
    int j_min=(boundary);
    int j_max=a[1]-boundary;
    int k_min=(boundary);
    int k_max=a[2]-boundary;

    for(int i = i_min; i < i_max ; i++) {
      for(int j = j_min; j < j_max; j++) {
        for(int k = k_min; k < k_max; k++) {
          map_new(i,j,k) = 1;
        } //for k
      } // for j
    } // for i
  } // constructor
  af::versa<double, af::c_grid<3> > result() {return map_new;}

}; // class zero_boundary_box_map

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_MASK_UTILS_H
