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
  int firstx, lastx, firsty, lasty;
  interval_list persistent_x_tiles, persistent_y_tiles;
 public:
  inline
  scanbox_tiling(){}
  inline
  scanbox_tiling(const int& firstx, const int& lastx,
                 const int& firsty, const int& lasty):
                 firstx(firstx), lastx(lastx), firsty(firsty), lasty(lasty),
                 persistent_x_tiles(interval_list()),
                 persistent_y_tiles(interval_list()){}

  inline
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
        result.push_back(interval(int(x*window), int((x+1)*window)-1));
     }

     return result;
  }

  inline
  virtual
  interval_ptr
  x_tiles(const int& interval) {
    if (persistent_x_tiles.size()==0){
      persistent_x_tiles = generate_normal_spacing(firstx,lastx,interval);
    }
    return persistent_x_tiles.begin();
  }

  inline
  virtual
  interval_ptr
  y_tiles(const int& interval) {
    if (persistent_y_tiles.size()==0){
      persistent_y_tiles = generate_normal_spacing(firsty,lasty,interval);
    }
    return persistent_y_tiles.begin();
  }

  inline
  interval_ptr
  x_end() const {
    return persistent_x_tiles.end();
  }

  inline
  interval_ptr
  y_end() const {
    return persistent_y_tiles.end();
  }

  inline
  void reset(){
    persistent_x_tiles = interval_list();
    persistent_y_tiles = interval_list();
  }

  inline virtual ~scanbox_tiling(){}
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

#include <scitbx/array_family/flex_types.h>

//For CXI CS Pad detector.  Explicitly define the rectangular active areas ahead of time.
class scanbox_tiling_explicit : public scanbox_tiling {
 size_t tile_count,peripheral_margin,internal_tile_state;
 scitbx::af::flex_int tiles;
 std::vector<int> tile_lookup;

 public:
  scanbox_tiling_explicit(scitbx::af::flex_int explicit_tiling,int const& peripheral_margin):
    tile_count(explicit_tiling.size()/4),peripheral_margin(peripheral_margin),
    internal_tile_state(0),tiles(explicit_tiling)
  {}

  interval_list
  generate_defined_x_tiles(int const& width){
     internal_tile_state= 0;
     tile_lookup.clear();
     interval_list result;
     for (int tile = 0 ; tile < tile_count; ++tile){
       int start = tiles[4*tile]+peripheral_margin;
       int finish= tiles[4*tile+2]-peripheral_margin-1;
       interval I = interval(start,finish);
       int available_width = I.size();
       int n_intervals = available_width/width;
       int next = start;
       for (int iival= 0; iival < n_intervals; ++iival){
         int increment = int((double(iival+1)/double(n_intervals))*(finish-start));
         result.push_back( interval(next,start + increment) );
         tile_lookup.push_back(tile);
         next = start + increment + 1;
       }
     }
     return result;
  }

  interval_list
  generate_defined_y_tiles(int const& width){
     interval_list result;
     int tile = tile_lookup[internal_tile_state];
     int start = tiles[4*tile+1]+peripheral_margin;
     int finish= tiles[4*tile+3]-peripheral_margin-1;
     interval I = interval(start,finish);
     int available_width = I.size();
     int n_intervals = available_width/width;
     int next = start;
     for (int iival= 0; iival < n_intervals; ++iival){
       int increment = int((double(iival+1)/double(n_intervals))*(finish-start));
       result.push_back( interval(next,start + increment) );
       next = start + increment + 1;
     }

     internal_tile_state+=1;
     return result;
  }

  virtual
  interval_ptr
  x_tiles(const int& scanbox_width) {
    persistent_x_tiles = generate_defined_x_tiles(scanbox_width);
    return persistent_x_tiles.begin();
  }

  virtual
  interval_ptr
  y_tiles(const int& scanbox_width) {
    persistent_y_tiles = generate_defined_y_tiles(scanbox_width);
    return persistent_y_tiles.begin();
  }

};

}

#endif //scanbox_h
