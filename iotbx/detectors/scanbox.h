#ifndef SCANBOX_H
#define SCANBOX_H

#include <vector>
#include <boost/shared_ptr.hpp>

namespace Distl {

struct interval {
    int first, last;
    interval(const int& first, const int&last):first(first),last(last){}
    inline size_t size()const{return last-first+1;}
};

typedef std::vector<interval> interval_list;
typedef interval_list::const_iterator interval_ptr;

class scanbox_tiling {
 public:
  interval_list persistent_x_tiles, persistent_y_tiles;
  int firstx, lastx, firsty, lasty;
 public:
  scanbox_tiling(const int& firstx, const int& lastx,
                 const int& firsty, const int& lasty):
                 firstx(firstx), lastx(lastx), firsty(firsty), lasty(lasty),
                 persistent_x_tiles(interval_list()),
                 persistent_y_tiles(interval_list()){}

  interval_list
  generate_normal_spacing(const int& first, const int& last,
                          const int& interv) const {
     //printf("Detector range from %5d to %5d\n",first, last);
     interval_list result;

     //legacy formula.
     /*
     for (int x = first; x < last - interv; x+=interv) {
        result.push_back(interval(x, x+interv-1));
     }
     */

     //expansive formula -- cover the region given by intervals that
     // may not be the same size.
     int n_intervals = (last-first+1)/interv;
     //SCITBX_EXAMINE(n_intervals);

     if (n_intervals == 0){
       // Avoid divide-by-zero;
       // Require that scanbox size >= input box size parameter interv
       return result;
     }

     double window = (double)(1+last-first)/n_intervals;
     for (int x = 0; x < n_intervals; ++x) {
        result.push_back(interval(x*window, int((x+1)*window)-1));
     }

     return result;
  }

  virtual
  interval_ptr
  x_tiles(const int& interval) {
    if (persistent_x_tiles.size()==0){
      persistent_x_tiles = generate_normal_spacing(firstx,lastx,interval);
    }
    return persistent_x_tiles.begin();
  }

  virtual
  interval_ptr
  y_tiles(const int& interval) {
    if (persistent_y_tiles.size()==0){
      persistent_y_tiles = generate_normal_spacing(firsty,lasty,interval);
    }
    return persistent_y_tiles.begin();
  }
  interval_ptr
  x_end() const {
    return persistent_x_tiles.end();
  }
  interval_ptr
  y_end() const {
    return persistent_y_tiles.end();
  }
  void reset(){
    persistent_x_tiles = interval_list();
    persistent_y_tiles = interval_list();
  }


};

typedef boost::shared_ptr<scanbox_tiling> ptr_tiling;

class scanbox_tiling_pilatus6M : public scanbox_tiling {
 public:
  scanbox_tiling_pilatus6M(const int& firstx, const int& lastx,
                 const int& firsty, const int& lasty):
                 scanbox_tiling(firstx,lastx,firsty,lasty){}

  interval_list
  generate_pilatus_spacing(const int& module_width, const int& interspacing,
                           const int& n_modules, const int& interv,
                           const int& allowable_first, const int& allowable_last) const {
     interval_list result;

     //special Pilatus formula -- cover the entire region
     //but ignore the interspacings with inactive pixels

     int n_intervals = (module_width-1)/interv;
     double window = (double)(module_width)/n_intervals;
     for (int mod = 0; mod < n_modules; ++mod) {
       for (int x = 0; x < n_intervals; ++x) {
        int module_origin = mod * (module_width+interspacing);
        int potential_first = module_origin+(int)(x*window);
        int potential_last  = module_origin+(int)((x+1)*window)-1;
        if (potential_first >= allowable_last) {continue;}
        if (potential_last  <= allowable_first) {continue;}
        result.push_back(interval(potential_first,potential_last));
       }
    }
    return result;
  }

  virtual
  interval_ptr
  x_tiles(const int& scanbox_width) {
    if (persistent_x_tiles.size()==0){
      // x axis:  2527 total pixels
      persistent_x_tiles = generate_pilatus_spacing(195, 17, 12, scanbox_width, firstx, lastx);
    }
    return persistent_x_tiles.begin();
  }

  virtual
  interval_ptr
  y_tiles(const int& scanbox_width) {
    if (persistent_y_tiles.size()==0){
      // y axis:  2463 total pixels
      persistent_y_tiles = generate_pilatus_spacing(487, 7, 5, scanbox_width, firsty, lasty);
    }
    return persistent_y_tiles.begin();
  }

};

class scanbox_tiling_pilatus2M : public scanbox_tiling_pilatus6M {
 public:
  scanbox_tiling_pilatus2M(const int& firstx, const int& lastx,
                 const int& firsty, const int& lasty):
                 scanbox_tiling_pilatus6M(firstx,lastx,firsty,lasty){}

  virtual
  interval_ptr
  x_tiles(const int& scanbox_width) {
    if (persistent_x_tiles.size()==0){
      // x axis:  1679 total pixels
      persistent_x_tiles = generate_pilatus_spacing(195, 17, 8, scanbox_width, firstx, lastx);
    }
    return persistent_x_tiles.begin();
  }

  virtual
  interval_ptr
  y_tiles(const int& scanbox_width) {
    if (persistent_y_tiles.size()==0){
      // y axis:  1475 total pixels
      persistent_y_tiles = generate_pilatus_spacing(487, 7, 3, scanbox_width, firsty, lasty);
    }
    return persistent_y_tiles.begin();
  }

};

}

#endif //scanbox_h
