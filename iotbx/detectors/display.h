/* -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*- */

#ifndef DETECTORS_IMAGE_DISPV_H
#define DETECTORS_IMAGE_DISPV_H

#include <vector>
#include <algorithm>
#include <limits>
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/mat3.h>
#include <scitbx/mat2.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/graphics_utils/colors.h>
#include <iotbx/detectors/context/spot_xy_convention.h>
#include <scitbx/random.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/math/polygon_intersection.h>

namespace af = scitbx::af;
namespace iotbx { namespace detectors { namespace display {

struct Color {
  af::tiny<int,3> RGB;
  inline Color(int r,int g,int b):RGB(r,g,b){}
  inline af::tiny<double,3>
  as_unit_rgb(){
    af::tiny<double,3> rv;
    for (int i=0;i<3;++i){rv[i]=(RGB[i]/255.);}
    return rv;
  }
};

class ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int&,const int&){return true;}
  inline virtual bool is_active_area_by_linear_index(const int&){return true;}
};

typedef boost::shared_ptr<ActiveAreaDefault> ptr_area;

class ActiveAreaPilatus6M: public ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int& x,const int& y){
    /*
    x, vertical, slow, size1: 12 blocks of 195, separated by 11 stripes of 17, giving 2527.
    y  horizont, fast, size2:  5 blocks of 487, separated by 4 stripes of 7, giving 2463.
    Takes 0.02 seconds extra to apply this test to entire Pilatus 6M image.
    */
    return ( (x%212<195) && (y%494<487) );
  }
  inline virtual bool is_active_area_by_linear_index(const int& i){
    return is_active_area(i/2463,i%2463);
  }
};
class ActiveAreaPilatus2M: public ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int& x,const int& y){
    /*
    x, vertical, slow, size1:  8 blocks of 195, separated by 7 stripes of 17, giving 1679.
    y  horizont, fast, size2:  3 blocks of 487, separated by 2 stripes of 7, giving 1475.
    */
    return ( (x%212<195) && (y%494<487) );
  }
  inline virtual bool is_active_area_by_linear_index(const int& i){
    return is_active_area(i/1475,i%1475);
  }
};
class ActiveAreaPilatus300K: public ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int& x,const int& y){
    /*
    x, vertical, slow, size1:  3 blocks of 195, separated by 2 stripes of 17, giving 619.
    y  horizont, fast, size2:  1 blocks of 487, separated by 0 stripes of 7, giving 487.
    */
    return ( (x%212<195) && (y%494<487) );
  }
  inline virtual bool is_active_area_by_linear_index(const int& i){
    return is_active_area(i/487,i%487);
  }
};

template <int FastModules>
class ActiveAreaEigerX: public ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int& x,const int& y){
    /*
    EigerX-16M, FastModules=4
    x, vertical, slow, size1:  8 blocks of 514, separated by 7 stripes of 37, giving 4371.
    y  horizont, fast, size2:  4 blocks of 1030, separated by 3 stripes of 10, giving 4150.

    EigerX-9M, FastModules=3
    x, vertical, slow, size1:  6 blocks of 514, separated by 5 stripes of 37, giving 3269.
    y  horizont, fast, size2:  3 blocks of 1030, separated by 2 stripes of 10, giving 3110.

    EigerX=4M, FastModules=2
    x, vertical, slow, size1:  4 blocks of 514, separated by 3 stripes of 37, giving 2167.
    y  horizont, fast, size2:  2 blocks of 1030, separated by 1 stripes of 10, giving 2070.

    EigerX-1M, FastModules=1
    x, vertical, slow, size1:  2 blocks of 514, separated by 1 stripes of 37, giving 1065.
    y  horizont, fast, size2:  1 blocks of 1030, separated by 0 stripes of 10, giving 1030.

    EigerX-500K, FastModules=1 (no inactive areas, do not instantiate)
    x, vertical, slow, size1:  1 blocks of 514, separated by 0 stripes of 37, giving 514.
    y  horizont, fast, size2:  1 blocks of 1030, separated by 0 stripes of 10, giving 1030.
    */
    return ( (x%551<514) && (y%1040<1030) );
  }
  inline virtual bool is_active_area_by_linear_index(const int& i){
    return is_active_area(i/(FastModules*1030+(FastModules-1)*10),i%(FastModules*1030+(FastModules-1)*10));
  }
};

template <int FastModules>
class ActiveAreaEiger2X: public ActiveAreaDefault {
public:
  inline virtual bool is_active_area(const int& x,const int& y){
    /*
    Eiger2X-16M, FastModules=4
    x, vertical, slow, size1:  8 blocks of 512, separated by 7 stripes of 38, giving 4362.
    y  horizont, fast, size2:  4 blocks of 1028, separated by 3 stripes of 12, giving 4148.

    Eiger2X-9M, FastModules=3
    x, vertical, slow, size1:  6 blocks of 512, separated by 5 stripes of 38, giving 3262.
    y  horizont, fast, size2:  3 blocks of 1028, separated by 2 stripes of 12, giving 3108.

    Eiger2X=4M, FastModules=2
    x, vertical, slow, size1:  4 blocks of 512, separated by 3 stripes of 38, giving 2162.
    y  horizont, fast, size2:  2 blocks of 1028, separated by 1 stripes of 12, giving 2068.

    Eiger2X-1M, FastModules=1
    x, vertical, slow, size1:  2 blocks of 512, separated by 1 stripes of 38, giving 1062.
    y  horizont, fast, size2:  1 blocks of 1028, separated by 0 stripes of 12, giving 1028.

    Eiger2X-500K, FastModules=1 (no inactive areas, do not instantiate)
    x, vertical, slow, size1:  1 blocks of 512, separated by 0 stripes of 38, giving 512.
    y  horizont, fast, size2:  1 blocks of 1028, separated by 0 stripes of 12, giving 1028.
    */
    return ( (x%550<512) && (y%1040<1028) );
  }
  inline virtual bool is_active_area_by_linear_index(const int& i){
    return is_active_area(i/(FastModules*1040+(FastModules-1)*12),i%(FastModules*1040+(FastModules-1)*12));
  }
};


inline int iround(double const& x)
    {if (x < 0) return static_cast<int>(x-0.5);
      return static_cast<int>(x+.5);}

#define COLOR_GRAY 0
#define COLOR_RAINBOW 1
#define COLOR_HEAT 2
#define COLOR_INVERT 3
template <typename DataType = int>
class FlexImage {

public:
  typedef af::versa< DataType, af::flex_grid<> > array_t;
  typedef DataType data_t;
  array_t rawdata;  // original image
  af::versa<int, af::c_grid<3> > channels; // half-size 3-channel cont./bright. adjusted
  af::versa<int, af::c_grid<2> > export_m;
  int export_size_uncut1;
  int export_size_uncut2;
  int export_size_cut1;
  int export_size_cut2;
  int export_anchor_x;
  int export_anchor_y;
  const int nchannels;
  int color_scheme_state;
  bool show_untrusted;

  inline
  af::versa<DataType, af::c_grid<2> > raw_to_sampled(
    const af::versa<DataType, af::c_grid<2> >& raw) const {
    af::c_grid<2> rawsize = raw.accessor();
    af::c_grid<2> sample_size(rawsize[0]/binning, rawsize[1]/binning);

    af::versa<DataType, af::c_grid<2> > sampled(sample_size);

/*    if (binning==2) {
      for (std::size_t i = 0; i < sample_size[0]; ++i) {
        for (std::size_t j = 0; j < sample_size[1]; ++j) {
          sampled(i,j) = (raw(2*i,2*j) + raw(2*i+1,2*j+1))/2;
        }
      }
    } else */ if (binning==1) { return raw; }
    else {
      //binning is a larger power of two.  Take the max pixel for each chunk
      std::vector<DataType> candidate_max;
      for (std::size_t i = 0; i < sample_size[0]; ++i) {
        for (std::size_t j = 0; j < sample_size[1]; ++j) {
          for (std::size_t isample = 0; isample < binning; ++isample){
            for (std::size_t jsample = 0; jsample < binning; ++jsample){
              candidate_max.push_back(raw(binning*i+isample,binning*j+jsample));
            }
          }
          sampled(i,j) = *std::max_element(candidate_max.begin(),
                                           candidate_max.end());
          SCITBX_ASSERT(candidate_max.size()==binning*binning);
          candidate_max.clear();
          SCITBX_ASSERT(candidate_max.size()==0);
        }
      }
    }
    return sampled;
  }

  inline
  af::versa<int, af::c_grid<2> > bright_contrast(
    af::versa<DataType, af::c_grid<2> > raw) const {
      const double outscale = 256;

      af::versa<int, af::c_grid<2> > z(raw.accessor());

      bool has_pilatus_inactive_flag = false;
      ptr_area detector_location = ptr_area(new ActiveAreaDefault());
      if (vendortype=="Pilatus-6M") {
        detector_location = ptr_area(new ActiveAreaPilatus6M());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Pilatus-2M") {
        detector_location = ptr_area(new ActiveAreaPilatus2M());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Pilatus-300K") {
        detector_location = ptr_area(new ActiveAreaPilatus300K());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger-16M") {
        detector_location = ptr_area(new ActiveAreaEigerX<4>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger-9M") {
        detector_location = ptr_area(new ActiveAreaEigerX<3>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger-4M") {
        detector_location = ptr_area(new ActiveAreaEigerX<2>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger-1M") {
        detector_location = ptr_area(new ActiveAreaEigerX<1>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger2-16M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<4>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger2-9M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<3>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger2-4M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<2>());
        has_pilatus_inactive_flag = true;
      } else if (vendortype=="Eiger2-1M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<1>());
        has_pilatus_inactive_flag = true;
      }
      for (std::size_t i=0; i<raw.accessor()[0]; i++) {
        int slow = binning * i;
        int idx_slow = slow * rawdata.accessor().focus()[1];
        int idx_i = i * raw.accessor()[1];
        for (std::size_t j=0; j< raw.accessor()[1]; j++) {
          int fast = binning * j;
          int idx = idx_slow + fast;
          int idx_ij = idx_i + j;
          if (detector_location->is_active_area(slow, fast)) {
            //fractional input value:
            double corrected = raw[idx_ij]*correction;
            double outvalue  = outscale * ( 1.0 - corrected );
            if (has_pilatus_inactive_flag && raw[idx_ij]==-2){
              //an ad hoc flag to aimed at coloring red the inactive pixels (-2) on Pilatus
              z[idx_ij]=1000;  //flag value used is out of range of the 256-pt grey scale but in range int datatype.
              continue;
            } else if (raw[idx_ij]== std::numeric_limits<int>::min()){
              //flag value for explicitly-masked pixels
              z[idx_ij]=1000;
              if (has_pilatus_inactive_flag){
                //to avoid confusion, reset flagged pixels for Pilatus-like detectors to -2
                raw[idx_ij]=-2;
              }
              continue;
            } else if (raw[idx_ij]>saturation){
              z[idx_ij]=2000;  //an ad hoc flag to aimed at coloring yellow the saturated pixels
              continue;
            }
            if (outvalue<0.0)            { z[idx_ij]=0; }
            else if (outvalue>=outscale) { z[idx_ij]=int(outscale)-1; }
            else                         { z[idx_ij] = int(outvalue); }
          }
        }
      }
      return z;
  }

  inline
  double global_bright_contrast() const {

      ptr_area detector_location = ptr_area(new ActiveAreaDefault());
      if (vendortype=="Pilatus-6M") {
        detector_location = ptr_area(new ActiveAreaPilatus6M());
      } else if (vendortype=="Pilatus-2M") {
        detector_location = ptr_area(new ActiveAreaPilatus2M());
      } else if (vendortype=="Pilatus-300K") {
        detector_location = ptr_area(new ActiveAreaPilatus300K());
      } else if (vendortype=="Eiger-16M") {
        detector_location = ptr_area(new ActiveAreaEigerX<4>());
      } else if (vendortype=="Eiger-9M") {
        detector_location = ptr_area(new ActiveAreaEigerX<3>());
      } else if (vendortype=="Eiger-4M") {
        detector_location = ptr_area(new ActiveAreaEigerX<2>());
      } else if (vendortype=="Eiger-1M") {
        detector_location = ptr_area(new ActiveAreaEigerX<1>());
      } else if (vendortype=="Eiger2-16M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<4>());
      } else if (vendortype=="Eiger2-9M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<3>());
      } else if (vendortype=="Eiger2-4M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<2>());
      } else if (vendortype=="Eiger2-1M") {
        detector_location = ptr_area(new ActiveAreaEiger2X<1>());
      }

      af::shared<data_t> raw_active;

      //first pass through data calculate average
      for (std::size_t i = 0; i < rawdata.accessor().focus()[0]; i++) {
        for (std::size_t j = 0; j < rawdata.accessor().focus()[1]; j++) {
          if (detector_location->is_active_area(i, j)) {
            raw_active.push_back(rawdata(i,j));
          }
        }
      }

      /*
      data_t active_mean = af::mean(raw_active.const_ref());
      SCITBX_EXAMINE(active_mean);
      data_t active_min= af::min(raw_active.const_ref());
      SCITBX_EXAMINE(active_min);
      data_t active_max = af::max(raw_active.const_ref());
      SCITBX_EXAMINE(active_max);
      */
      std::size_t active_count = raw_active.size();
      double PERCENTILE_TARGET = 0.90; // targets the 90th percentile pixel value.
      std::size_t nth_offset = PERCENTILE_TARGET*active_count;
      std::nth_element(raw_active.begin(), raw_active.begin()+nth_offset,
                       raw_active.end());

      if (vendortype=="Pilatus-6M") {
        SCITBX_ASSERT( (active_count == 60*195*487 || active_count == 5*195*487) );
      } else if (vendortype=="Pilatus-2M") {
        SCITBX_ASSERT( (active_count == 24*195*487 || active_count == 3*195*487) );
      } else if (vendortype=="Pilatus-300K") {
        SCITBX_ASSERT( (active_count == 3*195*487) );
      } else if (vendortype=="Eiger-16M") {
        SCITBX_ASSERT( (active_count == 32*514*1030 || active_count == 4*514*1030) );
      }

      double percentile = *(raw_active.begin()+nth_offset);
      double adjlevel = 0.4;

      return (percentile>0.) ? brightness * adjlevel/percentile : brightness / 5.0;
  }

  int binning; // either 1 (unbinned) or a power of 2 (binned)
  std::string vendortype;
  double brightness, correction;
  DataType saturation;
  double zoom;
  bool supports_rotated_tiles_antialiasing_recommended; //whether the client viewer should apply extra antialiasing

  /* Relationship between binning & zoom for this class:
        zoom level  zoom    binning   magnification factor
        ----------  ----    -------   --------------------
         -3         0.125     8        1/8 x
         -2         0.25      4        1/4 x
         -1         0.5       2        1/2 x
         0          1         1        1x
         1          2         1        2x
         2          4         1        4x
         etc.
  */

public:
  inline
  FlexImage(array_t rawdata, const int& power_of_two,
            const double& brightness = 1.0,
            const DataType& saturation = 1.0,
            const bool& show_untrusted = false,
            const int& color_scheme_state = COLOR_GRAY):
    rawdata(rawdata),
    brightness(brightness),saturation(saturation), nchannels(4),
    supports_rotated_tiles_antialiasing_recommended(false),
    color_scheme_state(color_scheme_state),show_untrusted(show_untrusted),
    binning(power_of_two){}

  inline
  FlexImage(array_t rawdata, const int& power_of_two,
            const std::string& vendortype,
            const double& brightness = 1.0,
            const DataType& saturation = 65535,
            const bool& show_untrusted = false,
            const int& color_scheme_state = COLOR_GRAY):
    brightness(brightness),
    saturation(saturation),
    rawdata(rawdata),
    nchannels(4),
    supports_rotated_tiles_antialiasing_recommended(false),
    color_scheme_state(color_scheme_state),
    binning(power_of_two),
    vendortype(vendortype),
    show_untrusted(show_untrusted){
    //Assert that binning is a power of two
    SCITBX_ASSERT ( binning > 0 && (binning & (binning-1) )==0 );
    zoom = 1./ binning;
    export_size_uncut1 = size1()/binning;
    export_size_uncut2 = size2()/binning;
    channels = af::versa<int, af::c_grid<3> >(af::c_grid<3>(nchannels,
                            export_size_uncut1,export_size_uncut2),af::init_functor_null<int>());
    correction = global_bright_contrast();
  }

  //change raw data to reflect different view
  inline
  void spot_convention(const int & conv){
    if (conv==0) {return;}

    int size1 = rawdata.accessor().focus()[0]; //the slow size
    int size2 = rawdata.accessor().focus()[1]; //the fast size
    array_t z(af::flex_grid<>(size1,size2));
    DataType* zac   = z.begin();
    DataType* raw   = rawdata.begin();

    // The only conventions used so far are 0 & 2, so no need to make this
    // general yet.
    if (conv==2) {
      for (int i = 0; i<size1; ++i){
        for (int j = 0; j<size2; ++j){
          zac[size1*i+j] = raw[size1*(size1-1-i)+j];
        }
      }
    }
    //generalize using the XY convention mechanism
    // later modify the mosflm orientation writer to work with spot conventions like 4.
    else {
      iotbx::detectors::context::spot_xy_convention XY(size1,size2,1,conv);
      //third arg is dummy "1"

      for (int i = 0; i<size1; ++i){
        for (int j = 0; j<size2; ++j){
          scitbx::vec2<int> ptr(i,j);
          scitbx::vec2<int> transformed = XY.call(
            &ptr, scitbx::type_holder<int>());
          zac[size1*i+j] = raw[size1*transformed[0]+transformed[1]];
        }
      }
    }
    rawdata = z;
  }

  inline int size1() const {return rawdata.accessor().focus()[0];}
  inline int size2() const {return rawdata.accessor().focus()[1];}

  inline
  void setWindow(
    const double& wxafrac, const double& wyafrac,const double& fraction){
    int apply_zoom = (binning == 1)? zoom : 1;
    export_size_cut1 = export_size_uncut1*fraction*apply_zoom;
    export_size_cut2 = export_size_uncut2*fraction*apply_zoom;
    export_m = af::versa<int, af::c_grid<2> >(
       af::c_grid<2>(export_size_cut1,export_size_cut2));
    //Compute integer anchor index based on input fractional dimension:
    export_anchor_x = export_size_uncut1*wxafrac*apply_zoom;
    export_anchor_y = export_size_uncut2*wyafrac*apply_zoom;
  }

  inline
  void setWindowCart(
    const double& xtile, const double& ytile,const double& fraction){
      // fractional coverage of image is defined on dimension 1 (slow) only,
      // allowing square subarea display of a rectangular detector image
      // tile coverage is specified in cartesian rather than fractional image coordinates
    int apply_zoom = (binning == 1)? zoom : 1;
    export_size_cut1 = iround((double(size1())/binning)*fraction*apply_zoom);
    export_size_cut2 = iround((double(size2())/binning)*fraction*apply_zoom*(double(size1())/double(size2())));
    export_m = af::versa<int, af::c_grid<2> >(
       af::c_grid<2>(export_size_cut1,export_size_cut2));
    //Compute integer anchor index based on input fractional dimension:
    export_anchor_x = iround(export_size_uncut1*xtile*fraction*apply_zoom);
    export_anchor_y = iround(export_size_uncut2*ytile*fraction*apply_zoom*(double(size1())/double(size2())));
  }

  inline int ex_size1() const {return export_m.accessor()[0];}
  inline int ex_size2() const {return export_m.accessor()[1];}

  inline void
  setZoom( const int& zoom_level) {
    zoom = std::pow(2.,double(zoom_level));

    int potential_binning = int(std::ceil(1./double(zoom)));
    if (potential_binning != binning) {
      binning = potential_binning;
      export_size_uncut1 = size1()/binning;
      export_size_uncut2 = size2()/binning;
      channels = af::versa<int, af::c_grid<3> >(af::c_grid<3>(nchannels,
                            export_size_uncut1,export_size_uncut2));
      adjust(color_scheme_state);
    }
  }

  void adjust(int color_scheme=COLOR_GRAY) {
    color_scheme_state = color_scheme;
    using scitbx::graphics_utils::hsv2rgb;
    using scitbx::graphics_utils::get_heatmap_color;
    af::versa<DataType, af::c_grid<2> > sam(rawdata, af::c_grid<2>(rawdata.accessor()));
    af::versa<int, af::c_grid<2> > data = bright_contrast(
                                             raw_to_sampled(sam));
    for (int i=0; i< nchannels; ++i) {
      int idx_i = i * export_size_uncut1;
      for (int j=0; j< export_size_uncut1; ++j) {
        int idx_ij = (idx_i + j) * export_size_uncut2;
        int idx_j = j * export_size_uncut2;
        for (int k=0; k<export_size_uncut2; ++k) {
          int idx_jk = idx_j + k;
          int idx_ijk = idx_ij + k;

          if (data[idx_jk] == 1000){
            //ad hoc method to color inactive pixels red
            if (color_scheme == COLOR_GRAY || color_scheme == COLOR_INVERT) {
              if (i==0) {channels[idx_ijk] = 254;}
              else {channels[idx_ijk] = 1;}
            } else { // or black, if displaying a rainbow or heatmap
              channels[idx_ijk] = 0;
            }
            continue;
          } else if (data[idx_jk] == 2000){
            //ad hoc method to color saturated pixels yellow
            if (color_scheme == COLOR_GRAY || color_scheme == COLOR_INVERT) {
              if (i==0 || i==1) {channels[idx_ijk] = 254;}
              else {channels[idx_ijk] = 1;}
            } else if (color_scheme == COLOR_RAINBOW) {
              // or white, if displaying a rainbow
              channels[idx_ijk] = 255;
            } else { // heatmap
              if (i == 1) {
                channels[idx_ijk] = 255;
              } else {
                channels[idx_ijk] = 0;
              }
            }
            continue;
          }
          if (color_scheme == COLOR_GRAY) {
            channels[idx_ijk] = data[idx_jk];
          } else if (color_scheme == COLOR_INVERT) {
            channels[idx_ijk] = 255. - data[idx_jk];
          } else if (color_scheme == COLOR_RAINBOW) { // rainbow
            double h = 255 * std::pow(((double) data[idx_jk]) / 255., 0.5);
            scitbx::vec3<double> rgb = hsv2rgb(h, 1.0, 1.0);
            channels[idx_ijk] = (int) (rgb[i] * 255.);
          } else { // heatmap
            double ratio = std::pow((double) (255. - data[idx_jk]) / 255., 2);
            scitbx::vec3<double> rgb = get_heatmap_color(ratio);
            channels[idx_ijk] = (int) (rgb[i] * 255.);
          }
        }
      }
    }
  }

  inline
  af::versa<int, af::c_grid<2> > channel(const int& c){
    for (int i=export_anchor_x; i< export_anchor_x+export_size_cut1; ++i) {
        for (int j=export_anchor_y; j< export_anchor_y+export_size_cut2; ++j) {
          export_m(i,j) = channels(c,i,j);
        }
    }

    return export_m;

  }

  inline
  void point_overlay(const int& x, const int& y, const Color& c){
    if ( x >= 0 && x < size1() && y >= 0 && y < size2() ) {
    //must be a better way to do this test within af::shared framework
      int sampled_x = x/binning;
      int sampled_y = y/binning;
      for (int i=0; i<3; ++i){
        channels(i,sampled_x,sampled_y) = c.RGB[i];
      }
    }
  }

  /* there's probably considerable room to improve this if the is_valid_index()
  were not tested each time.  Instead break the export window into 4 sections,
  only one of which has valid data in it.  Populate the export_s for the
  valid section first, then put pink in the other sections.
  */
  inline
  void prep_string(){
    typedef af::c_grid<3> t_C;
    const t_C& acc = channels.accessor();
    ptr_area detector_location = ptr_area(new ActiveAreaDefault());
    if (vendortype=="Pilatus-6M") {
      detector_location = ptr_area(new ActiveAreaPilatus6M());
    } else if (vendortype=="Pilatus-2M") {
      detector_location = ptr_area(new ActiveAreaPilatus2M());
    } else if (vendortype=="Pilatus-300K") {
      detector_location = ptr_area(new ActiveAreaPilatus300K());
    } else if (vendortype=="Eiger-16M") {
      detector_location = ptr_area(new ActiveAreaEigerX<4>());
    } else if (vendortype=="Eiger-9M") {
      detector_location = ptr_area(new ActiveAreaEigerX<3>());
    } else if (vendortype=="Eiger-4M") {
      detector_location = ptr_area(new ActiveAreaEigerX<2>());
    } else if (vendortype=="Eiger-1M") {
      detector_location = ptr_area(new ActiveAreaEigerX<1>());
    } else if (vendortype=="Eiger2-16M") {
      detector_location = ptr_area(new ActiveAreaEiger2X<4>());
    } else if (vendortype=="Eiger2-9M") {
      detector_location = ptr_area(new ActiveAreaEiger2X<3>());
    } else if (vendortype=="Eiger2-4M") {
      detector_location = ptr_area(new ActiveAreaEiger2X<2>());
    } else if (vendortype=="Eiger2-1M") {
      detector_location = ptr_area(new ActiveAreaEiger2X<1>());
    }
    export_s = "";
    export_s.reserve(export_size_cut1*export_size_cut2*3);
    if (zoom > 1.){
    for (int zoom_x=export_anchor_x;
         zoom_x< export_anchor_x+export_size_cut1;
         ++zoom_x) {
        for (int zoom_y=export_anchor_y;
             zoom_y< export_anchor_y+export_size_cut2;
             ++zoom_y) {
            int i = zoom_x / zoom;
            int j = zoom_y / zoom;
            if (acc.is_valid_index(0,i,j) &&
                detector_location->is_active_area(i*binning,j*binning)){
              for (int c=0; c<3; ++c) {
                export_s.push_back((char)channels(c,i,j));
              }
            } else {
                 //pixels not in-bounds within the image get a pink tint
                 export_s.push_back((char)255);
                 export_s.push_back((char)228);
                 export_s.push_back((char)228);
            }
        }
    }
    } else {
    for (int i=export_anchor_x; i< export_anchor_x+export_size_cut1; ++i) {
        for (int j=export_anchor_y; j< export_anchor_y+export_size_cut2; ++j) {
            if (acc.is_valid_index(0,i,j) &&
                detector_location->is_active_area(i*binning,j*binning)){
              for (int c=0; c<3; ++c) {
                export_s.push_back((char)channels(c,i,j));
              }
            } else {
                 //pixels not in-bounds within the image get a pink tint
                 export_s.push_back((char)255);
                 export_s.push_back((char)228);
                 export_s.push_back((char)228);
            }
        }
    }
    }
  }

  inline
  void prep_string_monochrome(){
    typedef af::c_grid<3> t_C;
    const t_C& acc = channels.accessor();
    export_s = "";
    export_s.reserve(export_size_cut1*export_size_cut2);
    for (int i=export_anchor_x; i< export_anchor_x+export_size_cut1; ++i) {
        for (int j=export_anchor_y; j< export_anchor_y+export_size_cut2; ++j) {
            if (acc.is_valid_index(0,i,j)){
                 export_s.push_back((char)channels(0,i,j));
            } else {
                 //pixels not in-bounds within the image get a white tint
                 export_s.push_back((char)255);
            }
        }
    }
  }

  std::string export_s; //export data in string form; make public for readonly

  /**
   * Return data in bytes form
   */
  inline
  boost::python::object as_bytes(){
    // Convert to a python bytes object
    boost::python::object data_bytes(
      boost::python::handle<>(PyBytes_FromStringAndSize(export_s.c_str(), export_s.size()))
    );
    return data_bytes;
  }

  inline
  void circle_overlay(const double& pixel_size,
                      af::shared<scitbx::vec3<double> > centers,
                      const double& radius,
                      const double& thickness,
                      const Color& color){
    //centers is actually a list of circle center coordinates, given in mm
    //pixel_size is the conversion factor to integer point coordinates
    //radius, thickness are given in points
    typedef af::shared<scitbx::vec3<double> >::const_iterator afsi;
    afsi c_end=centers.end();
    for (double r = iround(radius - 0.5*thickness);
          r<iround(radius + 0.5*thickness); r+=1.) {
      double trial_angle = 0.9/r; //in radians; angle subtended by 0.9 pixel at radius r
      int increments = scitbx::constants::two_pi / trial_angle;
      double angle = scitbx::constants::two_pi / (increments - increments%4);
      for (double theta=0.0;theta<scitbx::constants::two_pi; theta+=angle){
        int delta_x = iround(r * std::cos(theta));
        int delta_y = iround(r * std::sin(theta));
        for (afsi c=centers.begin(); c!=c_end; ++c) {
          int x = 0.5 + ((*c)[0] / pixel_size); //rounding of a positive number
          int y = 0.5 + ((*c)[1] / pixel_size); //rounding of a positive number
          point_overlay(x+delta_x,y+delta_y,color);
        }
      }
    }
  }

};

class generic_flex_image: public FlexImage<double>{
 public:

  scitbx::vec3<double> axis;
  scitbx::mat3<double> rotation;
  scitbx::mat2<double> rotation2;

  scitbx::af::shared<scitbx::mat2<double> > transformations;
  scitbx::af::shared<scitbx::vec2<double> > translations;
  std::vector<int> windowed_readouts;

  int size1_readout;
  int size2_readout;

  typedef af::c_grid<3> t_C;
  t_C acc;

  // The ASIC:s stacked within rawdata must be padded in each
  // dimension to multiples of eight.  Hence, the unpadded size of a
  // readout is passed in size_readout1 and size_readout2.
  inline
  generic_flex_image(
    array_t rawdata, const int& power_of_two,
    const int& size1_readout,
    const int& size2_readout,
    const double& brightness = 1.0,
    const double& saturation = 1.0,
    const bool& show_untrusted = false,
    const int& color_scheme_state = COLOR_GRAY
    )
    : FlexImage<double>(rawdata, power_of_two, brightness, saturation,
                        show_untrusted, color_scheme_state),
      size1_readout(size1_readout),
      size2_readout(size2_readout)
  {
    supports_rotated_tiles_antialiasing_recommended=true;
    zoom = 1./ binning;
    export_size_uncut1 = size1()/binning;
    export_size_uncut2 = size2()/binning;
    channels = af::versa<int, af::c_grid<3> >(af::c_grid<3>(nchannels,
                            export_size_uncut1,export_size_uncut2),af::init_functor_null<int>());
    axis = scitbx::vec3<double>(0.,0.,1.);
    rotation = scitbx::math::r3_rotation::axis_and_angle_as_matrix<double>(axis,4.,true);
    rotation2 = scitbx::mat2<double>(rotation[0],rotation[1],rotation[3],rotation[4]);
    correction = 1.0; //requires followup brightness scaling; use separate function
  }

  inline
  void setWindow(
    const double& wxafrac, const double& wyafrac,const double& fraction){
    int apply_zoom = (binning == 1)? zoom : 1;
    //Compute integer anchor index based on input fractional dimension:
    export_anchor_x = export_size_uncut1*wxafrac*apply_zoom;
    export_anchor_y = export_size_uncut2*wyafrac*apply_zoom;

    export_size_cut1 = export_size_uncut1;
    export_size_cut2 = export_size_uncut2;

    // Calculate the limits of the output window by iterating over all tiles
    for (size_t k = 0; k < transformations.size(); k++) {
      for (size_t islow=0; islow <= size1_readout; islow+=size1_readout) {
      for (size_t ifast=0; ifast <= size2_readout; ifast+=size2_readout) {
        scitbx::vec2<double> point_p = tile_readout_to_picture(k,islow,ifast);
        export_size_cut1 = std::max(
          export_size_cut1, int(std::ceil(point_p[0])));
        export_size_cut2 = std::max(
          export_size_cut2, int(std::ceil(point_p[1])));
      }
      }
    }

    //Find out which readouts have any intersection with this Window
    windowed_readouts.clear();
    //Define the 4 corners of the Window rectangle A B C D, counterclockwise from top left,
    af::shared<scitbx::vec2<double> > window_polygon;
    /*A*/window_polygon.push_back(
      scitbx::vec2<double>(export_anchor_x - 1.  ,export_anchor_y - 1.));
    /*B*/window_polygon.push_back(
      scitbx::vec2<double>(export_size_cut1 + 1. ,export_anchor_y - 1.));
    /*C*/window_polygon.push_back(
      scitbx::vec2<double>(export_size_cut1 + 1. ,export_size_cut2 + 1. ));
    /*D*/window_polygon.push_back(
      scitbx::vec2<double>(export_anchor_x - 1.  ,export_size_cut2 + 1. ));

    for (size_t k = 0; k < transformations.size(); k++) { //loop through all readout tiles
      af::shared<scitbx::vec2<double> > readout_polygon;
      for (size_t islow=0; islow <= size1_readout; islow+=size1_readout) {
        for (size_t ifast=0; ifast <= size2_readout; ifast+=size2_readout) {
          scitbx::vec2<double> point_p = tile_readout_to_picture(k,islow,ifast);
          readout_polygon.push_back(point_p);
        }
      }
      std::swap<scitbx::vec2<double> >(readout_polygon[2],readout_polygon[3]);
      if (scitbx::math::convex_polygons_intersect_2D(window_polygon, readout_polygon)) {
        windowed_readouts.push_back(k);
      }
    }

    export_size_cut1 = iround((export_size_cut1/binning)*fraction*apply_zoom);
    export_size_cut2 = iround((export_size_cut2/binning)*fraction*apply_zoom);
    export_m = af::versa<int, af::c_grid<2> >(
       af::c_grid<2>(export_size_cut1,export_size_cut2));

  }

  // Identical to base class, except for computation of export_anchor &
  //   just-in-time calculation of window/readout intersections
  inline
  void setWindowCart(
    const double& xtile, const double& ytile,const double& fraction){
      // fractional coverage of image is defined on dimension 1 (slow) only,
      // allowing square subarea display of a rectangular detector image
      // tile coverage is specified in cartesian rather than fractional image coordinates
    int apply_zoom = (binning == 1)? zoom : 1;
    export_size_cut1 = iround((double(size1())/binning)*fraction*apply_zoom);
    export_size_cut2 = iround((double(size2())/binning)*fraction*apply_zoom*(double(size1())/double(size2())));

    export_m = af::versa<int, af::c_grid<2> >(
       af::c_grid<2>(export_size_cut1,export_size_cut2));
    //Compute integer anchor index based on input fractional dimension:
    export_anchor_x = xtile * export_size_cut2;
    export_anchor_y = ytile * export_size_cut1;

    //Find out which readouts have any intersection with this WindowCart
    windowed_readouts.clear();
    //Define the 4 corners of the WindowCart rectangle A B C D, counterclockwise from top left,
    af::shared<scitbx::vec2<double> > window_polygon;
    /*A*/window_polygon.push_back( scitbx::vec2<double>((export_anchor_x/zoom) - 1.                ,(export_anchor_y/zoom) - 1.                ));
    /*B*/window_polygon.push_back( scitbx::vec2<double>((((xtile+1) * export_size_cut2)/zoom) + 1. ,(export_anchor_y/zoom) - 1.                ));
    /*C*/window_polygon.push_back( scitbx::vec2<double>((((xtile+1) * export_size_cut2)/zoom) + 1. ,(((ytile+1) * export_size_cut1)/zoom) + 1. ));
    /*D*/window_polygon.push_back( scitbx::vec2<double>((export_anchor_x/zoom) - 1.                ,(((ytile+1) * export_size_cut1)/zoom) + 1. ));

    for (size_t k = 0; k < transformations.size(); k++) { //loop through all readout tiles
      af::shared<scitbx::vec2<double> > readout_polygon;

      for (size_t islow=0; islow <= size1_readout; islow+=size1_readout) {
      for (size_t ifast=0; ifast <= size2_readout; ifast+=size2_readout) {
        scitbx::vec2<double> point_p = tile_readout_to_picture(k,islow,ifast);
        readout_polygon.push_back(point_p);
      }
      }
      std::swap<scitbx::vec2<double> >(readout_polygon[2],readout_polygon[3]);
      if (scitbx::math::convex_polygons_intersect_2D(window_polygon, readout_polygon)) {
        windowed_readouts.push_back(k);
      }
    }

  }

  inline scitbx::vec2<double> tile_readout_to_picture_f(
    int const& tile, double const& islow, double const& ifast) const {

    scitbx::vec2<double> fpicture =
        transformations[tile].inverse() * (
          scitbx::vec2<double>(islow, ifast) - translations[tile]);

    return fpicture;
  }

  inline scitbx::vec2<double> tile_readout_to_picture(
    int const& tile, int const& islow, int const& ifast) const {

    return tile_readout_to_picture_f(
      tile, double(islow), double(ifast));
  }

  inline af::shared<double> picture_to_readout_f(double const& i,double const& j)
   const {
    af::shared<double> z;

    if (transformations.size() == 0) {
      scitbx::vec2<double> rdout = rotation2 * scitbx::vec2<double>(i,j);
      z.push_back(rdout[0]); z.push_back(rdout[1]);
      return z;
    }

    // The length of a padded readout in the slow dimension, assuming
    // all readouts are the same size, and that the number of
    // transformation matrices is equal to the number of readouts.
    const std::size_t dim_slow = size1() / transformations.size();

    for (size_t k = 0; k < transformations.size(); k++) {
      scitbx::vec2<double> rdout =
        transformations[k] * scitbx::vec2<double>(i, j) + translations[k];

      scitbx::vec2<int> irdout(iround(rdout[0]), iround(rdout[1]));
      if (irdout[0] >= 0 && irdout[0] < size1_readout &&
          irdout[1] >= 0 && irdout[1] < size2_readout) {

        // Since acc may be binned, irdout must take it into account.
        if (acc.is_valid_index(
              0, (k * dim_slow + irdout[0]) / binning, irdout[1] / binning)) {
          z.push_back(rdout[0]); z.push_back(rdout[1]); z.push_back(k);
          return z;
        }
      }
    }
    z.push_back(0); z.push_back(0); z.push_back(-1);
    return z;
  }
  inline scitbx::vec2<int> picture_to_readout(double const& i,double const& j)
    const {
    if (transformations.size() == 0) {
      //return scitbx::vec2<int>(iround(i),iround(j));
      scitbx::vec2<double> rdout = rotation2 * scitbx::vec2<double>(i,j);
      return scitbx::vec2<int>(iround(rdout[0]),iround(rdout[1]));
    }

    // The length of a binned and padded readout in the slow
    // dimension.
    const std::size_t dim_slow = size1() / transformations.size() / binning;

    for (size_t k = 0; k < windowed_readouts.size(); k++) {
      scitbx::vec2<double> rdout =
        transformations[windowed_readouts[k]] * scitbx::vec2<double>(i, j) +
        translations[windowed_readouts[k]] / binning;

      scitbx::vec2<int> irdout(iround(rdout[0]), iround(rdout[1]));
      if (irdout[0] >= 0 && irdout[0] < size1_readout / binning &&
          irdout[1] >= 0 && irdout[1] < size2_readout / binning) {

        irdout[0] += windowed_readouts[k] * dim_slow;
        if (acc.is_valid_index(0, irdout[0], irdout[1])){
          return irdout;
        }
      }
    }
    return scitbx::vec2<int>(-1, -1);
  }
  inline void prep_string(){
    //const t_C& acc = channels.accessor();
    acc = channels.accessor();
    export_s = "";
    export_s.reserve(export_size_cut1*export_size_cut2*3);
    //scitbx::vec3<int> datafocus = channels.accessor().focus();
    if (zoom <= 1.){
      for (int i=export_anchor_x; i< export_anchor_x+export_size_cut1; ++i) {
          for (int j=export_anchor_y; j< export_anchor_y+export_size_cut2; ++j) {

              scitbx::vec2<int> irdjrd = picture_to_readout(i,j);
              if (acc.is_valid_index(0,irdjrd[0],irdjrd[1])){
                for (int c=0; c<3; ++c) {
                  export_s.push_back((char)channels(c,irdjrd[0],irdjrd[1]));
                }
              } else {
                  //pixels not in-bounds within the image get a pink tint
                   export_s.push_back((char)255);
                   export_s.push_back((char)228);
                   export_s.push_back((char)228);
              }
          }
      }
    } else {
      for (int zoom_x=export_anchor_x;
           zoom_x< export_anchor_x+export_size_cut1;
           ++zoom_x) {
          double i = zoom_x / zoom;
          for (int zoom_y=export_anchor_y;
               zoom_y< export_anchor_y+export_size_cut2;
               ++zoom_y) {
              double j = zoom_y / zoom;
              scitbx::vec2<int> irdjrd = picture_to_readout(i,j);
              if (acc.is_valid_index(0,irdjrd[0],irdjrd[1])){
                for (int c=0; c<3; ++c) {
                  export_s.push_back((char)channels(c,irdjrd[0],irdjrd[1]));
                }
              } else {
                  //pixels not in-bounds within the image get a pink tint
                   export_s.push_back((char)255);
                   export_s.push_back((char)228);
                   export_s.push_back((char)228);
              }
          }
      }
    }

    application_specific_anti_aliasing();
  }
  inline void application_specific_anti_aliasing(){
    // no implementation at present
  }

  inline void add_transformation_and_translation(
    const scitbx::mat2<double>& T, const scitbx::vec2<double>& t)
  {
    transformations.push_back(T);
    translations.push_back(t);
  }

  inline
  void followup_brightness_scale(){

      //first pass through data calculate average
      double qave = af::mean(rawdata.const_ref());
      //std::cout<<"ave shown pixel value is "<<qave<<std::endl;

      //second pass calculate histogram
      data_t* data_ptr = rawdata.begin();
      int hsize=100;
      array_t histogram(hsize);
      std::size_t data_sz = rawdata.size();
      double bins_per_pixel_unit = (hsize/2)/qave;
      int temp;
      for (std::size_t i = 0; i < data_sz; i++) {
          temp = int(bins_per_pixel_unit*(*data_ptr++));
          if (temp<0){histogram[0]+=1;}
          else if (temp>=hsize){histogram[hsize-1]+=1;}
          else {histogram[temp]+=1;}
      }

      //third pass calculate 90%
      double percentile=0;
      double accum=0;
      for (std::size_t i = 0; i<hsize; i++) {
        accum+=histogram[i];
        if (accum > 0.9*rawdata.size()) { percentile=i*qave/(hsize/2); break; }
      }
      //std::cout<<"the 90-percentile pixel value is "<<percentile<<std::endl;

      double adjlevel = 0.4;
      correction = (percentile>0.) ? brightness * adjlevel/percentile : brightness / 5.0;
  }
};

}}} //namespace

#endif //DETECTORS_IMAGE_DISPV_H
