#ifndef CCTBX_MARKED_GRID_POINTS_H
#define CCTBX_MARKED_GRID_POINTS_H

#include <scitbx/array_family/accessors/c_grid.h>


namespace cctbx { namespace maptbx {

//! Map marked_grid_points analysis.


class marked_grid_points {

private:
  af::versa<bool, af::c_grid<3> > final_result;
  af::shared<scitbx::vec3<int> > grid_points;
  af::tiny<int, 3> map_dimensions;

public:
  template <typename MapType>
  marked_grid_points(
    af::const_ref<MapType, af::flex_grid<> > const& map_data,
    int const& every_nth_point)
  {
    CCTBX_ASSERT(map_data.accessor().nd() == 3);
    CCTBX_ASSERT(map_data.accessor().all().all_gt(0));
    map_dimensions = af::adapt(map_data.accessor().all());
    int pointer_empty=0, pointer_current=0, cur_reg = 0;
    af::const_ref<MapType, af::c_grid<3> > data_ref(
        map_data.begin(),
        af::c_grid<3>(af::adapt(map_data.accessor().all())));

    // do something result
    int istart = every_nth_point/2;
    for (int i = istart; i < map_dimensions[0]; i+=every_nth_point) {
      for (int j = istart; j < map_dimensions[1]; j+=every_nth_point) {
        for (int k = istart; k < map_dimensions[2]; k+=every_nth_point) {
          if (map_data(i,j,k) ) {
            // save this point
            grid_points.push_back(scitbx::vec3<int>(i,j,k));
          }

        }
      }
    }

  }

  af::shared<scitbx::vec3<int> > result() {return grid_points;}

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MARKED_GRID_POINTS_H
