/* -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*- */

#ifndef DETECTORS_IMAGE_DISPV_H
#define DETECTORS_IMAGE_DISPV_H

#include <vector>
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

inline int iround(double const& x)
    {if (x < 0) return static_cast<int>(x-0.5);
      return static_cast<int>(x+.5);}

#define COLOR_GRAY 0
#define COLOR_RAINBOW 1
#define COLOR_HEAT 2
template <typename DataType = int>
class FlexImage {

public:
  typedef af::versa< DataType, af::flex_grid<> > array_t;
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

  inline
  af::versa<DataType, af::c_grid<2> > raw_to_sampled(
    const af::versa<DataType, af::c_grid<2> >& raw) const {
    af::c_grid<2> rawsize = raw.accessor();
    af::c_grid<2> sample_size(rawsize[0]/binning, rawsize[1]/binning);

    af::versa<DataType, af::c_grid<2> > sampled(sample_size);

    if (binning==2) {
      for (std::size_t i = 0; i < sample_size[0]; ++i) {
        for (std::size_t j = 0; j < sample_size[1]; ++j) {
          sampled(i,j) = (raw(2*i,2*j) + raw(2*i+1,2*j+1))/2;
        }
      }
    } else if (binning==1) { return raw; }
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

      ptr_area detector_location = ptr_area(new ActiveAreaDefault());
      if (vendortype=="Pilatus-6M") {
        detector_location = ptr_area(new ActiveAreaPilatus6M());
      } else if (vendortype=="Pilatus-2M") {
        detector_location = ptr_area(new ActiveAreaPilatus2M());
      } else if (vendortype=="Pilatus-300K") {
        detector_location = ptr_area(new ActiveAreaPilatus300K());
      }
      for (std::size_t i = 0; i < raw.size(); i++) {
        int fast = binning*(i%raw.accessor()[1]);
        int slow = binning*(i/raw.accessor()[1]);
        int idx = slow*rawdata.accessor().focus()[1]+ fast;
        if (detector_location->is_active_area_by_linear_index(idx)) {
          //fractional input value:
          double corrected = raw[i]*correction;
          double outvalue  = outscale * ( 1.0 - corrected );
          if (raw[i]==-2){
            //an ad hoc flag to aimed at coloring red the inactive pixels (-2) on Pilatus
            z[i]=1000;  //flag value used is out of range of the 256-pt grey scale but in range int datatype.
            continue;
          } else if (raw[i]>saturation){
            z[i]=2000;  //an ad hoc flag to aimed at coloring yellow the saturated pixels
            continue;
          }
          if (outvalue<0.0)            { z[i]=0; }
          else if (outvalue>=outscale) { z[i]=int(outscale)-1; }
          else                         { z[i] = int(outvalue); }
        }
      }
      return z;
  }

  inline
  double global_bright_contrast() const {
      /* yes, this could probably be made more efficient using the std::algorithm library */

      ptr_area detector_location = ptr_area(new ActiveAreaDefault());
      if (vendortype=="Pilatus-6M") {
        detector_location = ptr_area(new ActiveAreaPilatus6M());
      } else if (vendortype=="Pilatus-2M") {
        detector_location = ptr_area(new ActiveAreaPilatus2M());
      } else if (vendortype=="Pilatus-300K") {
        detector_location = ptr_area(new ActiveAreaPilatus300K());
      }
      std::size_t active_count = 0;

      //first pass through data calculate average
      double qave = 0;
      for (std::size_t i = 0; i < rawdata.size(); i++) {
        if (detector_location->is_active_area_by_linear_index(i)) {
          qave+=rawdata[i];
          active_count++;
        }
      }
      if (vendortype=="Pilatus-6M") {
        SCITBX_ASSERT( (active_count == 60*195*487 || active_count == 5*195*487) );
      } else if (vendortype=="Pilatus-2M") {
        SCITBX_ASSERT( (active_count == 24*195*487 || active_count == 3*195*487) );
      } else if (vendortype=="Pilatus-300K") {
        SCITBX_ASSERT( (active_count == 3*195*487) );
      }

      qave/=active_count;
      //std::cout<<"ave shown pixel value is "<<qave<<std::endl;

      //second pass calculate histogram
      int hsize=100;
      array_t histogram(hsize);
      for (std::size_t i = 0; i < rawdata.size(); i++) {
        if (detector_location->is_active_area_by_linear_index(i)) {
          int temp = int((hsize/2)*rawdata[i]/qave);
          if (temp<0){histogram[0]+=1;}
          else if (temp>=hsize){histogram[hsize-1]+=1;}
          else {histogram[temp]+=1;}
        }
      }

      //third pass calculate 90%
      double percentile=0;
      double accum=0;
      for (std::size_t i = 0; i<hsize; i++) {
        accum+=histogram[i];
        if (accum > 0.9*active_count) { percentile=i*qave/(hsize/2); break; }
      }
      //std::cout<<"the 90-percentile pixel value is "<<percentile<<std::endl;

      double adjlevel = 0.4;
      return (percentile>0.) ? brightness * adjlevel/percentile : 1.0;
  }

  int binning; // either 1 (unbinned) or a power of 2 (binned)
  std::string vendortype;
  double brightness, correction;
  int saturation;
  double zoom;
  bool use_antialiasing; //whether the client viewer should apply extra antialiasing

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
  FlexImage(array_t rawdata,const double& brightness = 1.0, const double& saturation = 1.0):
    rawdata(rawdata),
    brightness(brightness),saturation(saturation), nchannels(4),use_antialiasing(false),
    color_scheme_state(COLOR_GRAY){}

  inline
  FlexImage(array_t rawdata, const int& power_of_two,
            const std::string& vendortype, const double& brightness = 1.0,
            int const& saturation = 65535):
    brightness(brightness),
    saturation(saturation),
    rawdata(rawdata),
    nchannels(4),
    use_antialiasing(false),
    color_scheme_state(COLOR_GRAY),
    binning(power_of_two),
    vendortype(vendortype){
    //Assert that binning is a power of two
    SCITBX_ASSERT ( binning > 0 && (binning & (binning-1) )==0 );
    zoom = 1./ binning;
    export_size_uncut1 = size1()/binning;
    export_size_uncut2 = size2()/binning;
    channels = af::versa<int, af::c_grid<3> >(af::c_grid<3>(nchannels,
                            export_size_uncut1,export_size_uncut2));
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
      for (int j=0; j< export_size_uncut1; ++j) {
        for (int k=0; k<export_size_uncut2; ++k) {
          if (data(j,k) == 1000){
            //ad hoc method to color inactive pixels red
            if (color_scheme == COLOR_GRAY) {
              if (i==0) {channels(i,j,k) = 254;}
              else {channels(i,j,k) = 1;}
            } else { // or black, if displaying a rainbow or heatmap
              channels(i,j,k) = 0;
            }
            continue;
          } else if (data(j,k) == 2000){
            //ad hoc method to color saturated pixels yellow
            if (color_scheme == COLOR_GRAY) {
              if (i==0 || i==1) {channels(i,j,k) = 254;}
              else {channels(i,j,k) = 1;}
            } else if (color_scheme == COLOR_RAINBOW) {
              // or white, if displaying a rainbow
              channels(i,j,k) = 255;
            } else { // heatmap
              if (i == 1) {
                channels(i,j,k) = 255;
              } else {
                channels(i,j,k) = 0;
              }
            }
            continue;
          }
          if (color_scheme == COLOR_GRAY) {
            channels(i,j,k) = data(j,k);
          } else if (color_scheme == COLOR_RAINBOW) { // rainbow
            double h = 255 * std::pow(((double) data(j,k)) / 255., 0.5);
            scitbx::vec3<double> rgb = hsv2rgb(h, 1.0, 1.0);
            channels(i,j,k) = (int) (rgb[i] * 255.);
          } else { // heatmap
            double ratio = std::pow((double) (255. - data(j,k)) / 255., 2);
            scitbx::vec3<double> rgb = get_heatmap_color(ratio);
            channels(i,j,k) = (int) (rgb[i] * 255.);
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

  int size1_readout;
  int size2_readout;

  typedef af::c_grid<3> t_C;
  t_C acc;

  // The ASIC:s stacked within rawdata must be padded in each
  // dimension to multiples of eight.  Hence, the unpadded size of a
  // readout is passed in size_readout1 and size_readout2.
  inline
  generic_flex_image(
    array_t rawdata,
    const int& size1_readout,
    const int& size2_readout,
    const double& brightness = 1.0,
    const double& saturation = 1.0)
    : FlexImage<double>(rawdata,brightness, saturation),
      size1_readout(size1_readout),
      size2_readout(size2_readout)
  {
    use_antialiasing=true;
    binning=1;
    zoom = 1./ binning;
    export_size_uncut1 = size1()/binning;
    export_size_uncut2 = size2()/binning;
    channels = af::versa<int, af::c_grid<3> >(af::c_grid<3>(nchannels,
                            export_size_uncut1,export_size_uncut2));
    axis = scitbx::vec3<double>(0.,0.,1.);
    rotation = scitbx::math::r3_rotation::axis_and_angle_as_matrix<double>(axis,4.,true);
    rotation2 = scitbx::mat2<double>(rotation[0],rotation[1],rotation[3],rotation[4]);
    correction = 1.0;
  }

  // Identical to base class, except for computation of export_anchor.
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
  }

  inline scitbx::vec2<double> tile_readout_to_picture(
    int const& tile, int const& islow,int const& ifast) const {

    scitbx::vec2<double> fpicture =
        transformations[tile].inverse() * (
          scitbx::vec2<double>(islow, ifast) - translations[tile]);

    return fpicture;
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

    for (size_t k = 0; k < transformations.size(); k++) {
      scitbx::vec2<double> rdout =
        transformations[k] * scitbx::vec2<double>(i, j) +
        translations[k] / binning;

      scitbx::vec2<int> irdout(iround(rdout[0]), iround(rdout[1]));
      if (irdout[0] >= 0 && irdout[0] < size1_readout / binning &&
          irdout[1] >= 0 && irdout[1] < size2_readout / binning) {

        irdout[0] += k * dim_slow;
        if (acc.is_valid_index(0, irdout[0], irdout[1]))
          return irdout;
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
};

}}} //namespace

#endif //DETECTORS_IMAGE_DISPV_H
