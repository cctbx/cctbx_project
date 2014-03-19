#ifndef CCTBX_MAPTBX_CONNECTIVITY_H
#define CCTBX_MAPTBX_CONNECTIVITY_H

#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace maptbx {

//! Map connectivity analysis.

class connectivity {

public:
  af::versa<int, af::c_grid<3> > map_new;

  void
  get_six_neighbours(
    int const& x,
    int const& y,
    int const& z,
    af::tiny<int, 3> const& dimensions,
    af::versa<int, af::c_grid<2> > neighbours)
  {
    // x+-1
    neighbours(0,0) = ((x+1<dimensions[0]) ? x+1 : 0);
    neighbours(1,0) = ((x-1>=0) ? x-1 : dimensions[0]-1);
    neighbours(0,1) = neighbours(1,1) = y;
    neighbours(0,2) = neighbours(1,2) = z;
    // y+-1
    neighbours(2,1) = ((y+1<dimensions[1]) ? y+1 : 0);
    neighbours(3,1) = ((y-1>=0) ? y-1 : dimensions[1]-1);
    neighbours(2,0) = neighbours(3,0) = x;
    neighbours(2,2) = neighbours(3,2) = z;
    // z+-1
    neighbours(4,2) = ((z+1<dimensions[2]) ? z+1 : 0);
    neighbours(5,2) = ((z-1>=0) ? z-1 : dimensions[2]-1);
    neighbours(4,0) = neighbours(5,0) = x;
    neighbours(4,1) = neighbours(5,1) = y;
  }


  connectivity(
    af::const_ref<double, af::c_grid<3> > const& map_data,
    double const& threshold)
  {
    af::tiny<int, 3> a = map_data.accessor();
    int pointer_empty=0, pointer_current=0, cur_reg = 0;

    // estimating size of working array tempcoors. If this code fails with
    // segmentation fault or reveal a bug, here is the first place to look.
    // To make sure the cause is not here, just make tempcoors of map_data
    // size and delete boundaries check (91 and 97th lines).
    int maxside = ((a[0]>a[1]) ? a[0] : a[1]);
    maxside = maxside>a[2] ? maxside : a[2];
    int needed_size = 4*4*maxside*maxside;
    af::versa<int,af::c_grid<2> >  tempcoors(af::c_grid<2>(needed_size,3));

    map_new.resize(af::c_grid<3>(a), -1);
    af::versa<int, af::c_grid<2> > neighbours(af::c_grid<2>(6,3));

    for (int i = 0; i < a[0]; i++) {
      for (int j = 0; j < a[1]; j++) {
        for (int k = 0; k < a[2]; k++) {
          if (map_new(i,j,k)<0) {
            if (map_data(i,j,k) > threshold) {
              // got a new point, start filling
              cur_reg += 1;
              tempcoors(0,0) = i;
              tempcoors(0,1) = j;
              tempcoors(0,2) = k;
              map_new(i,j,k) = cur_reg;
              pointer_empty = 1;
              pointer_current = 0;
              while (pointer_empty != pointer_current) {
                get_six_neighbours(tempcoors(pointer_current,0),
                                   tempcoors(pointer_current,1),
                                   tempcoors(pointer_current,2),
                                   a,
                                   neighbours);
                for (int l = 0; l<6; l++) {
                  //processing 6 neighbours
                  int x = neighbours(l,0);
                  int y = neighbours(l,1);
                  int z = neighbours(l,2);
                  if (map_new(x,y,z)<0) {
                    if (map_data(x,y,z) > threshold) {
                      map_new(x,y,z) = cur_reg;
                      tempcoors(pointer_empty,0) = x;
                      tempcoors(pointer_empty,1) = y;
                      tempcoors(pointer_empty,2) = z;
                      pointer_empty += 1;
                      if (pointer_empty >= needed_size) pointer_empty = 0;
                    }
                    else map_new(x,y,z) = 0;
                  }
                }
                pointer_current += 1;
                if (pointer_current >= needed_size) pointer_current = 0;
              }
            }
            else map_new(i,j,k) = 0;
          }
        }
      }
    }
  }

  af::versa<int, af::c_grid<3> > result() {return map_new;}

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_CONNECTIVITY_H
