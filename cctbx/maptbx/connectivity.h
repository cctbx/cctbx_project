#ifndef CCTBX_MAPTBX_CONNECTIVITY_H
#define CCTBX_MAPTBX_CONNECTIVITY_H

#include <scitbx/array_family/accessors/c_grid.h>


namespace cctbx { namespace maptbx {

//! Map connectivity analysis.


class connectivity {

private:
  af::versa<int, af::c_grid<3> > map_new;
  af::shared<int> region_vols;
  af::tiny<int, 3> map_dimensions;
  af::shared<double> region_maximum_values;
  af::shared<scitbx::vec3<int> > region_maximum_coors;
  int n_regions;

  void
  get_six_neighbours(
    int const& x,
    int const& y,
    int const& z,
    af::shared<scitbx::vec3<int> > neighbours)
  {
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
  connectivity(
    af::const_ref<MapType, af::flex_grid<> > const& map_data,
    MapType const& threshold)
  {
    CCTBX_ASSERT(map_data.accessor().nd() == 3);
    CCTBX_ASSERT(map_data.accessor().all().all_gt(0));
    map_dimensions = af::adapt(map_data.accessor().all());
    int pointer_empty=0, pointer_current=0, cur_reg = 0;
    af::const_ref<MapType, af::c_grid<3> > data_ref(
        map_data.begin(),
        af::c_grid<3>(af::adapt(map_data.accessor().all())));
    // estimating size of working array tempcoors. If this code fails with
    // segmentation fault or reveal a bug, here is the first place to look.
    // To make sure the cause is not here, just make tempcoors of map_data
    // size and delete boundaries check (~L122, L143, look for
    // if ( ... >= needed_size)
    int maxside = ((map_dimensions[0]>map_dimensions[1]) ?
                        map_dimensions[0] : map_dimensions[1]);
    maxside = maxside>map_dimensions[2] ? maxside : map_dimensions[2];
    int needed_size = 4*4*maxside*maxside;
    af::shared<scitbx::vec3<int> > tempcoors(needed_size);
    // af::shared<scitbx::vec3<int> > tempcoors(needed_size);
    map_new.resize(af::c_grid<3>(map_dimensions), -1);
    af::shared<scitbx::vec3<int> > neighbours(6);
    region_vols.push_back(0);
    region_maximum_values.push_back(-10000000);
    region_maximum_coors.push_back(scitbx::vec3<int>(0,0,0));
    int v0 = 0, cur_reg_vol;

    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          if (map_new(i,j,k)<0) {
            if (data_ref(i,j,k) > threshold) {
              // got a new point, start filling
              cur_reg += 1;
              tempcoors[0] = scitbx::vec3<int> (i,j,k);
              map_new(i,j,k) = cur_reg;
              MapType cur_max_value = data_ref(i,j,k);
              scitbx::vec3<int> cur_max (i,j,k);
              cur_reg_vol = 1;
              pointer_empty = 1;
              pointer_current = 0;
              while (pointer_empty != pointer_current) {
                get_six_neighbours(tempcoors[pointer_current][0],
                                   tempcoors[pointer_current][1],
                                   tempcoors[pointer_current][2],
                                   neighbours);
                for (int l = 0; l<6; l++) {
                  //processing 6 neighbours
                  if (map_new(af::adapt(neighbours[l]))<0) {
                    if (data_ref(af::adapt(neighbours[l])) > threshold) {
                      map_new(af::adapt(neighbours[l])) = cur_reg;
                      cur_reg_vol += 1;
                      tempcoors[pointer_empty] = neighbours[l];
                      pointer_empty += 1;
                      if (pointer_empty >= needed_size) pointer_empty = 0;
                      if (data_ref(af::adapt(neighbours[l])) > cur_max_value)
                      {
                        cur_max_value = data_ref(af::adapt(neighbours[l]));
                        cur_max = neighbours[l];
                      }
                    }
                    else {
                      map_new(af::adapt(neighbours[l])) = 0;
                      v0 += 1;
                      if (data_ref(af::adapt(neighbours[l])) >
                          region_maximum_values[0])
                      {
                        region_maximum_values[0] =
                            data_ref(af::adapt(neighbours[l]));
                        region_maximum_coors[0] = neighbours[l];
                      }
                    }
                  }
                }
                pointer_current += 1;
                if (pointer_current >= needed_size) pointer_current = 0;
              }
              region_vols.push_back(cur_reg_vol);
              region_maximum_values.push_back(cur_max_value);
              region_maximum_coors.push_back(cur_max);
            }
            else {
              map_new(i,j,k) = 0;
              v0 += 1;
              if (data_ref(i,j,k) > region_maximum_values[0])
              {
                region_maximum_values[0] = data_ref(i,j,k);
                region_maximum_coors[0] = scitbx::vec3<int>(i,j,k);
              }
            }
          }
        }
      }
    }
    region_vols[0] = v0;
    n_regions = cur_reg;
  }

  af::versa<int, af::c_grid<3> > result() {return map_new;}
  af::shared<int> regions() {return region_vols;}
  af::shared<double> maximum_values() {return region_maximum_values;}
  af::shared<scitbx::vec3<int> > maximum_coors() {return region_maximum_coors;}
  af::versa<int, af::c_grid<3> > volume_cutoff_mask(int const& volume_cutoff)
  {
    af::versa<int, af::c_grid<3> > vol_mask;
    vol_mask.resize(af::c_grid<3>(map_dimensions), -1);
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          if (map_new(i,j,k)>0 && region_vols[map_new(i,j,k)] > volume_cutoff)
          {
            vol_mask(i,j,k)=1;
          }
          else
          {
            vol_mask(i,j,k)=0;
          }
        }
      }
    }
    return vol_mask;
  }

  af::versa<int, af::c_grid<3> >
  get_blobs_boundaries()
  {
    // boundaries - array[min/max, number_of_blob, x/y/z]
    af::versa<int, af::c_grid<3> > boundaries;
    boundaries.resize(af::c_grid<3>(2, n_regions+1, 3));
    // initialization
    for (int j = 0; j < n_regions+1; j++) {
      for (int k = 0; k < 3; k++) {
        boundaries(0,j,k) = 1000000;
        boundaries(1,j,k) = -1000000;
      }
    }
    // do the cycle
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          int n_reg = map_new(i,j,k);
          if (i < boundaries(0,n_reg,0)) {boundaries(0,n_reg,0) = i;}
          if (j < boundaries(0,n_reg,1)) {boundaries(0,n_reg,1) = j;}
          if (k < boundaries(0,n_reg,2)) {boundaries(0,n_reg,2) = k;}
          if (i > boundaries(1,n_reg,0)) {boundaries(1,n_reg,0) = i;}
          if (j > boundaries(1,n_reg,1)) {boundaries(1,n_reg,1) = j;}
          if (k > boundaries(1,n_reg,2)) {boundaries(1,n_reg,2) = k;}
        }
      }
    }
    return boundaries;
  }

  af::versa<bool, af::c_grid<3> >
  expand_mask(int id_to_expand, int expand_size)
  {
    CCTBX_ASSERT(expand_size > 0);
    CCTBX_ASSERT(id_to_expand >= 0);
    af::versa<bool, af::c_grid<3> > result;
    result.resize(af::c_grid<3>(map_dimensions), false);
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          if (map_new(i,j,k) == id_to_expand)
          {
            for (int ii=i-expand_size; ii<=i+expand_size; ii++) {
              for (int jj=j-expand_size; jj<=j+expand_size; jj++) {
                for (int kk=k-expand_size; kk<=k+expand_size; kk++) {
                  scitbx::vec3<int> adj_coor = put_coordinates_in_boundaries(ii,jj,kk);
                  // result(af::adapt(adj_coor)) = true;
                  result(adj_coor) = true;
                }
              }
            }
          }
        }
      }
    }
    return result;
  }

  af::versa<int, af::c_grid<3> >
  noise_elimination_two_cutoffs(
    connectivity const& connectivity_object_at_t1,
    int const& elimination_volume_threshold_at_t1,
    bool zero_all_interblob_region=false)
  {
    af::versa<int, af::c_grid<3> > result_mask;
    result_mask.resize(af::c_grid<3>(map_dimensions), 0);
    af::shared<int> fill_data(n_regions+1, 0);
    for (int i=1; i < connectivity_object_at_t1.n_regions+1; i++ )
    {
      if (connectivity_object_at_t1.region_vols[i] >
          elimination_volume_threshold_at_t1)
      {
        int good_reg_number_t2 =
            map_new(connectivity_object_at_t1.region_maximum_coors[i]);
        fill_data[good_reg_number_t2] = (good_reg_number_t2>0) ? 1 : 0;
      }
    }
    fill_data[0] = (zero_all_interblob_region) ? 0 : 1;
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          result_mask(i,j,k) = fill_data[map_new(i,j,k)];
        }
      }
    }
    return result_mask;
  }

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_CONNECTIVITY_H
