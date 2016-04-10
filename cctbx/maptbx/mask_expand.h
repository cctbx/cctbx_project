#ifndef CCTBX_MAPTBX_MASK_EXPAND_H
#define CCTBX_MAPTBX_MASK_EXPAND_H

#include <scitbx/array_family/accessors/c_grid.h>


namespace cctbx { namespace maptbx {

//! Map mask_expand analysis.


class marked_grid_points {

private:
  af::versa<bool, af::c_grid<3> > final_result;
  af::shared<scitbx::vec3<int> > grid_points;
  af::versa<int, af::c_grid<3> > map_new;
  af::shared<int> region_vols;
  af::tiny<int, 3> map_dimensions;
  af::shared<double> region_maximum_values;
  af::shared<scitbx::vec3<int> > region_maximum_coors;
  int n_regions;
  bool border_wrapping;

  int
  get_six_neighbours(
    int const& x,
    int const& y,
    int const& z,
    af::shared<scitbx::vec3<int> > neighbours)
  // returns number of neigbours, coordinates are in neighbours, should be
  // at least 6 long
  {
    int n=0;
    if (border_wrapping) {
      // x+-1
      neighbours[0][0] = ((x+1<map_dimensions[0]) ? x+1 : 0);
      neighbours[1][0] = ((x-1>=0) ? x-1 : map_dimensions[0]-1);
      neighbours[0][1] = neighbours[1][1] = y;
      neighbours[0][2] = neighbours[1][2] = z;
      // y+-1
      neighbours[2][1] = ((y+1<map_dimensions[1]) ? y+1 : 0);
      neighbours[3][1] = ((y-1>=0) ? y-1 : map_dimensions[1]-1);
      neighbours[2][0] = neighbours[3][0] = x;
      neighbours[2][2] = neighbours[3][2] = z;
      // z+-1
      neighbours[4][2] = ((z+1<map_dimensions[2]) ? z+1 : 0);
      neighbours[5][2] = ((z-1>=0) ? z-1 : map_dimensions[2]-1);
      neighbours[4][0] = neighbours[5][0] = x;
      neighbours[4][1] = neighbours[5][1] = y;
      n = 6;
    }
    else {
      //x+1
      if (x+1<map_dimensions[0]) {
        neighbours[n] = scitbx::vec3<int>(x+1, y, z);
        n++;
      }
      if (x-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x-1, y, z);
        n++;
      }
      if (y+1<map_dimensions[1]) {
        neighbours[n] = scitbx::vec3<int>(x, y+1, z);
        n++;
      }
      if (y-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y-1, z);
        n++;
      }
      if (z+1<map_dimensions[2]) {
        neighbours[n] = scitbx::vec3<int>(x, y, z+1);
        n++;
      }
      if (z-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y, z-1);
        n++;
      }
    }
    return n;
  }


  inline
  int get_coordinate_in_boundaries(int x, int map_size)
  {
    if (x > map_size-1) return x-map_size;
    if (x < 0) return map_size + x;
    return x;
  }

  scitbx::vec3<int>
  put_coordinates_in_boundaries(int x, int y, int z)
  {
    scitbx::vec3<int> result(0,0,0);
    result[0] = get_coordinate_in_boundaries(x, map_dimensions[0]);
    result[1] = get_coordinate_in_boundaries(y, map_dimensions[1]);
    result[2] = get_coordinate_in_boundaries(z, map_dimensions[2]);
    return result;
  }

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
class mask_expand {

private:
  af::versa<bool, af::c_grid<3> > final_result;
  af::versa<int, af::c_grid<3> > map_new;
  af::shared<int> region_vols;
  af::tiny<int, 3> map_dimensions;
  af::shared<double> region_maximum_values;
  af::shared<scitbx::vec3<int> > region_maximum_coors;
  int n_regions;
  bool border_wrapping;

  int
  get_six_neighbours(
    int const& x,
    int const& y,
    int const& z,
    af::shared<scitbx::vec3<int> > neighbours)
  // returns number of neigbours, coordinates are in neighbours, should be
  // at least 6 long
  {
    int n=0;
    if (border_wrapping) {
      // x+-1
      neighbours[0][0] = ((x+1<map_dimensions[0]) ? x+1 : 0);
      neighbours[1][0] = ((x-1>=0) ? x-1 : map_dimensions[0]-1);
      neighbours[0][1] = neighbours[1][1] = y;
      neighbours[0][2] = neighbours[1][2] = z;
      // y+-1
      neighbours[2][1] = ((y+1<map_dimensions[1]) ? y+1 : 0);
      neighbours[3][1] = ((y-1>=0) ? y-1 : map_dimensions[1]-1);
      neighbours[2][0] = neighbours[3][0] = x;
      neighbours[2][2] = neighbours[3][2] = z;
      // z+-1
      neighbours[4][2] = ((z+1<map_dimensions[2]) ? z+1 : 0);
      neighbours[5][2] = ((z-1>=0) ? z-1 : map_dimensions[2]-1);
      neighbours[4][0] = neighbours[5][0] = x;
      neighbours[4][1] = neighbours[5][1] = y;
      n = 6;
    }
    else {
      //x+1
      if (x+1<map_dimensions[0]) {
        neighbours[n] = scitbx::vec3<int>(x+1, y, z);
        n++;
      }
      if (x-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x-1, y, z);
        n++;
      }
      if (y+1<map_dimensions[1]) {
        neighbours[n] = scitbx::vec3<int>(x, y+1, z);
        n++;
      }
      if (y-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y-1, z);
        n++;
      }
      if (z+1<map_dimensions[2]) {
        neighbours[n] = scitbx::vec3<int>(x, y, z+1);
        n++;
      }
      if (z-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y, z-1);
        n++;
      }
    }
    return n;
  }


  inline
  int get_coordinate_in_boundaries(int x, int map_size)
  {
    if (x > map_size-1) return x-map_size;
    if (x < 0) return map_size + x;
    return x;
  }

  scitbx::vec3<int>
  put_coordinates_in_boundaries(int x, int y, int z)
  {
    scitbx::vec3<int> result(0,0,0);
    result[0] = get_coordinate_in_boundaries(x, map_dimensions[0]);
    result[1] = get_coordinate_in_boundaries(y, map_dimensions[1]);
    result[2] = get_coordinate_in_boundaries(z, map_dimensions[2]);
    return result;
  }

public:
  template <typename MapType>
  mask_expand(
    af::const_ref<MapType, af::flex_grid<> > const& map_data,
    MapType const& id_to_expand,
    bool wrapping=true)
  {
    CCTBX_ASSERT(map_data.accessor().nd() == 3);
    CCTBX_ASSERT(map_data.accessor().all().all_gt(0));
    map_dimensions = af::adapt(map_data.accessor().all());
    border_wrapping=wrapping;
    int pointer_empty=0, pointer_current=0, cur_reg = 0;
    af::const_ref<MapType, af::c_grid<3> > data_ref(
        map_data.begin(),
        af::c_grid<3>(af::adapt(map_data.accessor().all())));

    int expand_size=1;

    final_result.resize(af::c_grid<3>(map_dimensions), false);
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          if (map_data(i,j,k) == id_to_expand)
          {
            for (int ii=i-expand_size; ii<=i+expand_size; ii++) {
              for (int jj=j-expand_size; jj<=j+expand_size; jj++) {
                for (int kk=k-expand_size; kk<=k+expand_size; kk++) {
                  scitbx::vec3<int> adj_coor = put_coordinates_in_boundaries(ii,jj,kk);
                  // result(af::adapt(adj_coor)) = true;
                  final_result(adj_coor) = true;
                }
              }
            }
          }
        }
      }
    }
  }

  af::versa<bool, af::c_grid<3> > result() {return final_result;}

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_MASK_EXPAND_H
