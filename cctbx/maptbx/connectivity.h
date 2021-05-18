#ifndef CCTBX_MAPTBX_CONNECTIVITY_H
#define CCTBX_MAPTBX_CONNECTIVITY_H


// This needed for grid_symop.h inclusion for some reason
#include <cctbx/sgtbx/direct_space_asu/proto/asymmetric_unit.h>

#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <mmtbx/masks/grid_symop.h>
#include <cctbx/maptbx/asymmetric_map.h>

namespace cctbx { namespace maptbx {

using scitbx::af::int3;
//! Map connectivity analysis.


class connectivity {

private:
  af::versa<int, af::c_grid<3> > map_new;
  af::shared<int> region_vols;
  af::tiny<int, 3> map_dimensions;
  af::shared<double> region_maximum_values;
  af::shared<scitbx::vec3<int> > region_maximum_coors;
  int n_regions;
  bool border_wrapping;
  bool preprocess_shallow;
  std::vector<cctbx::sgtbx::grid_symop> symops;
  int3 uc_dims;


  inline
  bool within_cell(scitbx::vec3<int> &c)
  {
    return c[0]>=0 && c[1]>=0 && c[2]>=0 &&
      c[0]<map_dimensions[0] && c[1]<map_dimensions[1] && c[2]<map_dimensions[2];
  }

  int
  get_six_neighbours_sg(
    int const& x,
    int const& y,
    int const& z,
    af::shared<scitbx::vec3<int> > neighbours)
  {
    // border_wrapping is not used here, because we have a sg defined and
    // always 'wrap' using sg.
    // x+-1
    neighbours[0][0] = x+1;
    neighbours[1][0] = x-1;
    neighbours[0][1] = neighbours[1][1] = y;
    neighbours[0][2] = neighbours[1][2] = z;
    // y+-1
    neighbours[2][1] = y+1;
    neighbours[3][1] = y-1;
    neighbours[2][0] = neighbours[3][0] = x;
    neighbours[2][2] = neighbours[3][2] = z;
    // z+-1
    neighbours[4][2] = z+1;
    neighbours[5][2] = z-1;
    neighbours[4][0] = neighbours[5][0] = x;
    neighbours[4][1] = neighbours[5][1] = y;
    int n = 6;

    // now checking/moving them to an appropriate place
    for (int i=0; i<n; i++)
    {
      size_t symop=0;
      while (!within_cell(neighbours[i]) && symop<symops.size())
      {
        int3 sym_pos = symops[symop].apply_to(neighbours[i]);
        scitbx::int3 pos_in_cell(sym_pos);
        translate_into_cell(pos_in_cell, uc_dims);
        neighbours[i] = (pos_in_cell);
        symop++;
        // int reg_on_map = map_new(pos_in_cell);
        // std::cout << "    sym: " << sym_pos << " -> " << pos_in_cell
        //           << " region " << reg_on_map << "\n";
      }
    }

    return n;
  }

  std::vector<sgtbx::grid_symop> grid_symops(const asymmetric_map &amap) const
  {
    sgtbx::space_group group = amap.space_group();
    unsigned short order = group.order_z();
    CCTBX_ASSERT( order>0 );
    const int3 n = amap.unit_cell_grid_size();
    CCTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
    std::vector<sgtbx::grid_symop> symops;
    symops.reserve(order);
    for(size_t i=0; i<order; ++i)
    {
      sgtbx::grid_symop grsym( group(i), n );
      symops.push_back(grsym);
    }
    CCTBX_ASSERT( symops.size() == order );
    return symops;
  }

  void
  get_six_neighbours_asu(
      int3 pos,
      af::shared<int3> neighbours,
      const asymmetric_map &amap)
  {
    // std::cout << "  pos " << pos << " -> ";
    // border_wrapping is not used here, because we have a sg defined and
    // always 'wrap' using sg.

    int3 asu_begin = amap.box_begin();
    int3 asu_end = amap.box_end();
    int3 uc_d = amap.unit_cell_grid_size();
    std::vector<sgtbx::grid_symop> symops = grid_symops(amap);

    // init the neighbours
    for (int i=0; i<6; i++) neighbours[i] = pos;
    // x+-1
    neighbours[0][0] += 1;
    neighbours[1][0] -= 1;
    // // y+-1
    neighbours[2][1] += 1;
    neighbours[3][1] -= 1;
    // // z+-1
    neighbours[4][2] += 1;
    neighbours[5][2] -= 1;
    // int n = 6;

    for (int i=0; i<6; i++) {
      bool in_asu = neighbours[i].all_ge(asu_begin) && neighbours[i].all_le(asu_end);
      if (!in_asu) {
        // somehow move this one to asu...
        int3 asu_pos = neighbours[i];
        // int nsym = 0;
        // while (nsym < symops.size() && ! (asu_pos.all_ge(asu_begin) && asu_pos.all_le(asu_end))) {
        //   asu_pos = symops[nsym].apply_to(neighbours[i]);
        //   scitbx::int3 pos_in_cell(asu_pos);
        //   translate_into_cell(pos_in_cell,uc_d);
        //   std::cout << "!!! " << asu_pos << "("<< pos_in_cell << ")"<< "!!! ";
        //   nsym++;
        // }
        // CCTBX_ASSERT(asu_pos.all_ge(asu_begin) && asu_pos.all_le(asu_end));
        int n_res = 0;
        for (int nsym=0; nsym<symops.size(); nsym++) {
          asu_pos = symops[nsym].apply_to(neighbours[i]);
          scitbx::int3 pos_in_cell(asu_pos);
          translate_into_cell(pos_in_cell,uc_d);
          int3 pic(pos_in_cell);
          if (pic.all_ge(asu_begin) && pic.all_le(asu_end)) n_res +=1;
          // std::cout << "!!! " << asu_pos << "("<< pos_in_cell << ")"<< "!!! ";
        }
        neighbours[i] = asu_pos;
        // std::cout << "\nn_res " << n_res << "\n";
      }
      // std::cout << "    " << neighbours[i] << " " << in_asu << " ";
    }
    // std::cout << "\n";

    // // now checking/moving them to an appropriate place
    // for (int i=0; i<n; i++)
    // {
    //   size_t symop=0;
    //   while (!within_cell(neighbours[i]) && symop<symops.size())
    //   {
    //     int3 sym_pos = symops[symop].apply_to(neighbours[i]);
    //     scitbx::int3 pos_in_cell(sym_pos);
    //     translate_into_cell(pos_in_cell, uc_dims);
    //     neighbours[i] = (pos_in_cell);
    //     symop++;
    //     // int reg_on_map = map_new(pos_in_cell);
    //     // std::cout << "    sym: " << sym_pos << " -> " << pos_in_cell
    //     //           << " region " << reg_on_map << "\n";
    //   }
    // }

    // return n;
    // return 6;
  }



  int
  get_six_neighbours(
    int const& x,
    int const& y,
    int const& z,
    af::shared<scitbx::vec3<int> > neighbours,
    bool pad_with_negative=false)
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
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
        n++;
      }
      if (x-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x-1, y, z);
        n++;
      }
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
        n++;
      }

      if (y+1<map_dimensions[1]) {
        neighbours[n] = scitbx::vec3<int>(x, y+1, z);
        n++;
      }
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
        n++;
      }
      if (y-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y-1, z);
        n++;
      }
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
        n++;
      }
      if (z+1<map_dimensions[2]) {
        neighbours[n] = scitbx::vec3<int>(x, y, z+1);
        n++;
      }
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
        n++;
      }
      if (z-1>=0) {
        neighbours[n] = scitbx::vec3<int>(x, y, z-1);
        n++;
      }
      else if (pad_with_negative) {
        neighbours[n] = scitbx::vec3<int>(-1, -1, -1);
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
  int preprocessing_changed_voxels;
  int preprocessing_n_passes;

  template <typename MapType>
  connectivity(
    af::ref<MapType, af::c_grid<3> > map_data,
    MapType const& threshold,
    bool wrapping=true,
    bool preprocess_against_shallow=false) //: sg(cctbx::sgtbx::space_group("P 1"))
  {
    //new code
    // std::cout << "Old constructor.\n";
    // cctbx::sgtbx::space_group dummy_sg = cctbx::sgtbx::space_group("P 1");
    // sg = dummy_sg;
    // old code
    map_dimensions = map_data.accessor();
    border_wrapping=wrapping;
    preprocess_shallow=preprocess_against_shallow;
    int pointer_empty=0, pointer_current=0, cur_reg = 0;
    af::shared<scitbx::vec3<int> > neighbours(6);
    preprocessing_changed_voxels = 0;
    preprocessing_n_passes = 0;
    if (preprocess_shallow) {
      int n_changed;
      do {
        n_changed = 0;
        for (int i = 0; i < map_dimensions[0]; i++) {
          for (int j = 0; j < map_dimensions[1]; j++) {
            for (int k = 0; k < map_dimensions[2]; k++) {
              if (map_data(i,j,k) > threshold) {
                int n_neib = get_six_neighbours(i,j,k, neighbours, true);
                CCTBX_ASSERT(n_neib == 6);
                bool keep=true;
                int l=0;
                while (keep && l<3) {
                  MapType v1 = (neighbours[l*2][0]>=0) ? map_data(af::adapt(neighbours[l*2])) : threshold-1;
                  MapType v2 = (neighbours[l*2+1][0]>=0) ? map_data(af::adapt(neighbours[l*2+1])) : threshold-1;
                  if (v1 <= threshold && v2 <= threshold) keep=false;
                  l += 1;
                }
                if (!keep) {
                  map_data(i,j,k) = threshold-1;
                  n_changed += 1;
                }
              }
            }
          }
        }
        preprocessing_changed_voxels += n_changed;
        preprocessing_n_passes += 1;
      } while (n_changed > 0);
    }
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
    region_vols.push_back(0);
    region_maximum_values.push_back(-10000000);
    region_maximum_coors.push_back(scitbx::vec3<int>(0,0,0));
    int v0 = 0, cur_reg_vol;

    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          if (map_new(i,j,k)<0) {
            if (map_data(i,j,k) > threshold) {
              // got a new point, start filling
              cur_reg += 1;
              tempcoors[0] = scitbx::vec3<int> (i,j,k);
              map_new(i,j,k) = cur_reg;
              MapType cur_max_value = map_data(i,j,k);
              scitbx::vec3<int> cur_max (i,j,k);
              cur_reg_vol = 1;
              pointer_empty = 1;
              pointer_current = 0;
              while (pointer_empty != pointer_current) {
                int n_neib = get_six_neighbours(tempcoors[pointer_current][0],
                                   tempcoors[pointer_current][1],
                                   tempcoors[pointer_current][2],
                                   neighbours);
                for (int l = 0; l < n_neib; l++) {
                  //processing neighbours
                  if (map_new(af::adapt(neighbours[l]))<0) {
                    if (map_data(af::adapt(neighbours[l])) > threshold) {
                      map_new(af::adapt(neighbours[l])) = cur_reg;
                      cur_reg_vol += 1;
                      tempcoors[pointer_empty] = neighbours[l];
                      pointer_empty += 1;
                      if (pointer_empty >= needed_size) pointer_empty = 0;
                      if (map_data(af::adapt(neighbours[l])) > cur_max_value)
                      {
                        cur_max_value = map_data(af::adapt(neighbours[l]));
                        cur_max = neighbours[l];
                      }
                    }
                    else {
                      map_new(af::adapt(neighbours[l])) = 0;
                      v0 += 1;
                      if (map_data(af::adapt(neighbours[l])) >
                          region_maximum_values[0])
                      {
                        region_maximum_values[0] =
                            map_data(af::adapt(neighbours[l]));
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
              if (map_data(i,j,k) > region_maximum_values[0])
              {
                region_maximum_values[0] = map_data(i,j,k);
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

  template <typename MapType>
  connectivity(
    asymmetric_map &amap,
    MapType const& threshold,
    bool preprocess_against_shallow=false)

  {
    std::cout << "Asymmetric map constructor.\n";

    typedef scitbx::af::nested_loop<scitbx::int3> grid_iterator_t;
    typedef scitbx::af::ref<double, asymmetric_map::asu_grid_t  >  data_ref_t;
    data_ref_t amap_data=amap.data_ref();

    cctbx::sgtbx::space_group sg = amap.space_group();

    std::cout << "  amap box begin/end    " << amap.box_begin() << ";" << amap.box_end() << "\n";

    // map_dimensions = map_data.accessor();
    // border_wrapping=wrapping;
    preprocess_shallow=preprocess_against_shallow;
    int pointer_empty=0, pointer_current=0, cur_reg = 0;
    af::shared<scitbx::vec3<int> > neighbours(6);
    preprocessing_changed_voxels = 0;
    preprocessing_n_passes = 0;
    int n0=0, n1=0;
    if (preprocess_shallow) {
      int n_changed;
      do {
        n_changed = 0;
        // for (int i = 0; i < map_dimensions[0]; i++) {
        //   for (int j = 0; j < map_dimensions[1]; j++) {
        //     for (int k = 0; k < map_dimensions[2]; k++) {
        for(grid_iterator_t i3=amap.grid_begin(); !i3.over(); i3.incr() ) {
          int3 pos = i3();
          if (amap_data(pos) > threshold) {
            n1 += 1;
            af::shared<int3> neigh(6);

            // std::cout << "    data "<< amap_data(pos) << " > thresh\n";
            get_six_neighbours_asu(pos, neigh, amap);

            // CCTBX_ASSERT(n_neib == 6);
            // bool keep=true;
            // int l=0;
            // while (keep && l<3) {
            //   MapType v1 = (neighbours[l*2][0]>=0) ? map_data(af::adapt(neighbours[l*2])) : threshold-1;
            //   MapType v2 = (neighbours[l*2+1][0]>=0) ? map_data(af::adapt(neighbours[l*2+1])) : threshold-1;
            //   if (v1 <= threshold && v2 <= threshold) keep=false;
            //   l += 1;
            // }
            // if (!keep) {
            //   map_data(i,j,k) = threshold-1;
            //   n_changed += 1;
            // }
          } else n0+=1;
        }
        preprocessing_changed_voxels += n_changed;
        preprocessing_n_passes += 1;
      } while (n_changed > 0);
    }
    std::cout << "  C++ n0 " << n0 << "\n";
    std::cout << "  C++ n1 " << n1 << "\n";

  }

  template <typename MapType>
  connectivity(
    af::ref<MapType, af::c_grid<3> > map_data
    // MapType const& threshold
    // const cctbx::sgtbx::space_group &space_group,
    // const cctbx::sgtbx::space_group_type & space_group,
    // int3 uc_dimensions,
    // bool wrapping=true,
    // bool preprocess_against_shallow=false
    )
  {
    // addition
    std::cout << "Symmetry-aware constructor.\n";

    // uc_dims = uc_dimensions;
    // // sg = space_group;
    // unsigned short order = space_group.order_z();
    // CCTBX_ASSERT( order>0 );
    // // const int3 n = map_dimensions;
    // // CCTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
    // // std::vector<cctbx::sgtbx::grid_symop> symops;
    // symops.reserve(order);
    // for(size_t i=0; i<order; ++i)
    // {
    //   sgtbx::grid_symop grsym( space_group(i), uc_dims );
    //   symops.push_back(grsym);
    // }
    // CCTBX_ASSERT( symops.size() == order );
    // // end addition


    // map_dimensions = map_data.accessor();
    // border_wrapping=wrapping;
    // preprocess_shallow=preprocess_against_shallow;
    // int pointer_empty=0, pointer_current=0, cur_reg = 0;
    // af::shared<scitbx::vec3<int> > neighbours(6);
    // preprocessing_changed_voxels = 0;
    // preprocessing_n_passes = 0;
    // if (preprocess_shallow) {
    //   int n_changed;
    //   do {
    //     n_changed = 0;
    //     for (int i = 0; i < map_dimensions[0]; i++) {
    //       for (int j = 0; j < map_dimensions[1]; j++) {
    //         for (int k = 0; k < map_dimensions[2]; k++) {
    //           if (map_data(i,j,k) > threshold) {
    //             int n_neib = get_six_neighbours_sg(i,j,k, neighbours);
    //             CCTBX_ASSERT(n_neib == 6);
    //             bool keep=true;
    //             int l=0;
    //             while (keep && l<3) {
    //               MapType v1 = (neighbours[l*2][0]>=0) ? map_data(af::adapt(neighbours[l*2])) : threshold-1;
    //               MapType v2 = (neighbours[l*2+1][0]>=0) ? map_data(af::adapt(neighbours[l*2+1])) : threshold-1;
    //               if (v1 <= threshold && v2 <= threshold) keep=false;
    //               l += 1;
    //             }
    //             if (!keep) {
    //               map_data(i,j,k) = threshold-1;
    //               n_changed += 1;
    //             }
    //           }
    //         }
    //       }
    //     }
    //     preprocessing_changed_voxels += n_changed;
    //     preprocessing_n_passes += 1;
    //   } while (n_changed > 0);
    // }
    // // estimating size of working array tempcoors. If this code fails with
    // // segmentation fault or reveal a bug, here is the first place to look.
    // // To make sure the cause is not here, just make tempcoors of map_data
    // // size and delete boundaries check (~L122, L143, look for
    // // if ( ... >= needed_size)
    // int maxside = ((map_dimensions[0]>map_dimensions[1]) ?
    //                     map_dimensions[0] : map_dimensions[1]);
    // maxside = maxside>map_dimensions[2] ? maxside : map_dimensions[2];
    // int needed_size = 4*4*maxside*maxside;
    // af::shared<scitbx::vec3<int> > tempcoors(needed_size);
    // // af::shared<scitbx::vec3<int> > tempcoors(needed_size);
    // map_new.resize(af::c_grid<3>(map_dimensions), -1);
    // region_vols.push_back(0);
    // region_maximum_values.push_back(-10000000);
    // region_maximum_coors.push_back(scitbx::vec3<int>(0,0,0));
    // int v0 = 0, cur_reg_vol;

    // for (int i = 0; i < map_dimensions[0]; i++) {
    //   for (int j = 0; j < map_dimensions[1]; j++) {
    //     for (int k = 0; k < map_dimensions[2]; k++) {
    //       if (map_new(i,j,k)<0) {
    //         if (map_data(i,j,k) > threshold) {
    //           // got a new point, start filling
    //           cur_reg += 1;
    //           tempcoors[0] = scitbx::vec3<int> (i,j,k);
    //           map_new(i,j,k) = cur_reg;
    //           MapType cur_max_value = map_data(i,j,k);
    //           scitbx::vec3<int> cur_max (i,j,k);
    //           cur_reg_vol = 1;
    //           pointer_empty = 1;
    //           pointer_current = 0;
    //           while (pointer_empty != pointer_current) {
    //             int n_neib = get_six_neighbours_sg(tempcoors[pointer_current][0],
    //                                tempcoors[pointer_current][1],
    //                                tempcoors[pointer_current][2],
    //                                neighbours);
    //             for (int l = 0; l < n_neib; l++) {
    //               //processing neighbours
    //               if (map_new(af::adapt(neighbours[l]))<0) {
    //                 if (map_data(af::adapt(neighbours[l])) > threshold) {
    //                   map_new(af::adapt(neighbours[l])) = cur_reg;
    //                   cur_reg_vol += 1;
    //                   tempcoors[pointer_empty] = neighbours[l];
    //                   pointer_empty += 1;
    //                   if (pointer_empty >= needed_size) pointer_empty = 0;
    //                   if (map_data(af::adapt(neighbours[l])) > cur_max_value)
    //                   {
    //                     cur_max_value = map_data(af::adapt(neighbours[l]));
    //                     cur_max = neighbours[l];
    //                   }
    //                 }
    //                 else {
    //                   map_new(af::adapt(neighbours[l])) = 0;
    //                   v0 += 1;
    //                   if (map_data(af::adapt(neighbours[l])) >
    //                       region_maximum_values[0])
    //                   {
    //                     region_maximum_values[0] =
    //                         map_data(af::adapt(neighbours[l]));
    //                     region_maximum_coors[0] = neighbours[l];
    //                   }
    //                 }
    //               }
    //             }
    //             pointer_current += 1;
    //             if (pointer_current >= needed_size) pointer_current = 0;
    //           }
    //           region_vols.push_back(cur_reg_vol);
    //           region_maximum_values.push_back(cur_max_value);
    //           region_maximum_coors.push_back(cur_max);
    //         }
    //         else {
    //           map_new(i,j,k) = 0;
    //           v0 += 1;
    //           if (map_data(i,j,k) > region_maximum_values[0])
    //           {
    //             region_maximum_values[0] = map_data(i,j,k);
    //             region_maximum_coors[0] = scitbx::vec3<int>(i,j,k);
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
    // region_vols[0] = v0;
    // n_regions = cur_reg;
  }


  void experiment_with_symmetry(cctbx::sgtbx::space_group &space_group, int3 uc_dims)
  {
    std::cout << "Start experiment\n";
    unsigned short order = space_group.order_z();
    CCTBX_ASSERT( order>0 );
    const int3 n = map_dimensions;
    CCTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
    std::vector<cctbx::sgtbx::grid_symop> symops;
    symops.reserve(order);
    for(size_t i=0; i<order; ++i)
    {
      sgtbx::grid_symop grsym( space_group(i), uc_dims );
      symops.push_back(grsym);
    }
    std::cout << "SG order:" << order << "\n";
    CCTBX_ASSERT( symops.size() == order );

    std::vector<int3> coords;
    // coords.reserve(3);
    // int3 c0(0,0,17);
    // int3 c1(33,0,0);
    // int3 c2(0,33,0);

    int3 c0(1,0,17);
    int3 c1(2,0,17);
    int3 c2(3,0,17);

    coords.push_back(c0);
    coords.push_back(c1);
    coords.push_back(c2);

    for (int i=0; i<3; i++)
    {
      std::cout << "original coords: " << coords[i] << "\n";
      for(size_t symop=0; symop<symops.size(); symop++)
      {
        int3 sym_pos = symops[symop].apply_to(coords[i]);
        scitbx::int3 pos_in_cell(sym_pos);
        translate_into_cell(pos_in_cell, uc_dims);
        int reg_on_map = map_new(pos_in_cell);
        std::cout << "    sym: " << sym_pos << " -> " << pos_in_cell
                  << " region " << reg_on_map << "\n";
      }
    }
  }

  void merge_symmetry_related_regions(cctbx::sgtbx::space_group &space_group)
  {
    // copy-paste from asymmetric_map.cpp asymmetric_map::grid_symops()
    // sgtbx::space_group group = this->space_group();
    unsigned short order = space_group.order_z();
    CCTBX_ASSERT( order>0 );
    const int3 n = map_dimensions;
    CCTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
    std::vector<cctbx::sgtbx::grid_symop> symops;
    symops.reserve(order);
    for(size_t i=0; i<order; ++i)
    {
      sgtbx::grid_symop grsym( space_group(i), n );
      symops.push_back(grsym);
    }
    CCTBX_ASSERT( symops.size() == order );
    // return symops;
    // end of copy-paste
    int n_regions = region_vols.size();
    af::shared<int> remap_list(n_regions);
    for (int i = 0; i < n_regions; i++) remap_list[i] = -1;
    remap_list[0] = 0;
    int cur_region_to_fill = 0;
    for (int i = 1; i < n_regions; i++)
    {
      // std::cout << "loop # " << i <<"\n";
      // for (int j = 0; j < n_regions; j++) std::cout << remap_list[j] <<", ";
      // std::cout << "\n";
      if (remap_list[i]<0) // not assigned yet
      {
        cur_region_to_fill += 1;
        remap_list[i] = cur_region_to_fill;
        // now remap all symmetry-related regions
        int3 cur_coords = region_maximum_coors[i];
        // cur_coords[2] = 17;
        // std::cout << "  cur_coords " << cur_coords << "\n";
        // If 0th symop is always self, we can start from 1...
        for(size_t symop=0; symop<symops.size(); symop++)
        {
          bool mapped_with_self=false;
          int3 sym_pos = symops[symop].apply_to(cur_coords);
          scitbx::int3 pos_in_cell(sym_pos);
          translate_into_cell(pos_in_cell, n);
          int reg_on_map = map_new(pos_in_cell);
          // std::cout << "    sym: " << sym_pos << " -> " << pos_in_cell
          //           << " region " << reg_on_map << "\n";
          // safeguarding a little
          if (remap_list[reg_on_map] < 0 ) {
            mapped_with_self = true;
            remap_list[reg_on_map] = cur_region_to_fill;
          }
          else {
            if (reg_on_map < i)
            {
              if (mapped_with_self) {
                // This branch is not tested and not clear if needed at all,
                // i.e. if there are such cases.
                CCTBX_ASSERT(false);
                int rl_value = cur_region_to_fill;
                for (int j=0; j < n_regions; j++)
                  if (remap_list[j] == rl_value) {
                    remap_list[j] = reg_on_map;
                  }
              }
              else {
                remap_list[i] = reg_on_map;
              }
            }
            else if (reg_on_map > i)
            {
              mapped_with_self = true;
              remap_list[reg_on_map] = cur_region_to_fill;
            }
          }
        }
      }
    }
    // Here we are done with remap_list, sanity check:
    // std::cout << "remap_list\n";
    for (int i = 0; i < n_regions; i++)
    {
      // std::cout << remap_list[i] << ", ";
      CCTBX_ASSERT(remap_list[i] >=0);
    }
    // std::cout << "\n";
    // std::cout << "Correcting other info\n";
    // Need to refill map_new, cut and adjust region_vols,
    // cut region_maximum_values, region_maximum_coors
    int new_n_regions = -1;
    for (int i = 0; i < n_regions; i++)
      if (remap_list[i]> new_n_regions) {
        new_n_regions = remap_list[i];}
    new_n_regions += 1;
    region_vols.resize(new_n_regions);
    region_maximum_values.resize(0);
    region_maximum_coors.resize(0);
    for (int i=0; i<new_n_regions; i++) region_vols[i] = 0;
    for (int i = 0; i < map_dimensions[0]; i++) {
      for (int j = 0; j < map_dimensions[1]; j++) {
        for (int k = 0; k < map_dimensions[2]; k++) {
          map_new(i,j,k) = remap_list[map_new(i,j,k)];
          region_vols[map_new(i,j,k)] += 1;
        }
      }
    }
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
