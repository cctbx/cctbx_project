#ifndef CCTBX_MAPTBX_CONNECTIVITY_H
#define CCTBX_MAPTBX_CONNECTIVITY_H

#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace maptbx {

//! Placeholder for map connectivity analysis code.

class connectivity {
public:
  af::versa<double, af::c_grid<3> > map_new;

  connectivity(
    af::const_ref<double, af::c_grid<3> > const& map_data,
    double const& threshold)
  {
    af::tiny<int, 3> a = map_data.accessor();
    map_new.resize(af::c_grid<3>(a), 0);
    for (int i = 0; i < a[0]; i++) {
      for (int j = 0; j < a[1]; j++) {
        for (int k = 0; k < a[2]; k++) {
          double rho = map_data(i,j,k);
          if(rho < threshold) map_new(i,j,k) = 0;
          else                map_new(i,j,k) = 1;
        }
      }
    }
  }

  af::versa<double, af::c_grid<3> > result() {return map_new;}

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_CONNECTIVITY_H
