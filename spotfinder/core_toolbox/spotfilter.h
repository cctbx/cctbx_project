#ifndef SPOTFILTER_H
#define SPOTFILTER_H

#include <string>
#include <cmath>

#include <scitbx/error.h>//SCITBX_CHECK_POINT
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/constants.h>

#include <spotfinder/core_toolbox/libdistl.h>

namespace af = scitbx::af;
namespace spotfinder {
namespace distltbx {

typedef af::shared<Distl::spot> spot_list_t;
typedef Distl::spot w_spot;
typedef scitbx::vec3<double> point;
typedef scitbx::mat3<double> matrix;

struct SpotFilterAgent {
  double const pixel_size;
  double const xbeam;
  double const ybeam;
  double const distance;
  double const wavelength;
  af::shared<Distl::icering> icerings;
  af::shared<double> arguments;

  inline SpotFilterAgent (double pixel_size,
                          double xbeam,
                          double ybeam,
                          double distance,
                          double wavelength,
                          af::shared<Distl::icering> icerings):
      pixel_size(pixel_size),xbeam(xbeam),ybeam(ybeam),distance(distance),
      wavelength(wavelength),icerings(icerings){}
  void precompute_resolution(spot_list_t,
                             double const&,
                             af::shared<double>,
                             af::shared<double>);
  af::shared<int>
  filter(spot_list_t, af::shared<int>, std::string);
  af::shared<int>
  resolution_sort(spot_list_t, af::shared<int>);
  af::shared<int>
  resolution_sort_nztt(spot_list_t, af::shared<int>);
  af::shared<double>
  get_resolution(spot_list_t, af::shared<int>);
  af::shared<double>
  get_property(spot_list_t, af::shared<int>,std::string);
  void set_arguments(af::ref<double>const&);
  af::shared<int>
  order_by(spot_list_t, af::shared<int>, std::string);
};

class resolution_on_image {

  typedef scitbx::vec3<double> point;
  typedef scitbx::vec3<double> axis;
  typedef scitbx::mat3<double> matrix;
  typedef scitbx::af::tiny_plain<double,4> quaternion;
  matrix forward_W;
  matrix inverse_W;
  double xbeam_mm;
  double ybeam_mm;
  double distance_mm;
  double wavelength_ang;
  double twotheta_rad;
  matrix image_to_detector;
  point B_XTD_at_zero;
  axis direct_beam;
 public:
  inline resolution_on_image(
    double const& xbeam_mm,
    double const& ybeam_mm,
    double const& distance_mm,
    double const& wavelength_ang,
    double const& twotheta_rad):
    xbeam_mm(xbeam_mm),
    ybeam_mm(ybeam_mm),
    distance_mm(distance_mm),
    wavelength_ang(wavelength_ang),
    twotheta_rad(twotheta_rad),
    image_to_detector( -1., 0, 0, 0,1,0,0,0,1), //Not sure if this will be general
    direct_beam( 0.,0.,1.)  //beam travels in forward z direction
    {

      axis rotation_axis(0.,1.,0.); //assume for the moment that detector rotation (twotheta) is on y-hat.
      //express rotation through angle two theta as a unit quaternion;
      double h = twotheta_rad * 0.5;
      double c = std::cos(h); double s = std::sin(h);
      quaternion q(c, rotation_axis[0] * s, rotation_axis[1] * s, rotation_axis[2] * s);
      //re-express this rotation as a 3x3 matrix operator
      forward_W = matrix ( 2*(q[0]*q[0]+q[1]*q[1])-1, 2*(q[1]*q[2]-q[0]*q[3]), 2*(q[1]*q[3]+q[0]*q[2]),
                           2*(q[1]*q[2]+q[0]*q[3]), 2*(q[0]*q[0]+q[2]*q[2])-1, 2*(q[2]*q[3]-q[0]*q[1]),
                           2*(q[1]*q[3]-q[0]*q[2]), 2*(q[2]*q[3]+q[0]*q[1]), 2*(q[0]*q[0]+q[3]*q[3])-1);
      inverse_W = forward_W.inverse();

      //get detector
      point D = inverse_W * point(0,0,(distance_mm/std::cos(twotheta_rad)));
      point solve1 = image_to_detector * point (xbeam_mm,ybeam_mm,distance_mm);
      double B0x = D[0] - solve1[0];
      double B0y = D[1] - solve1[1];
      B_XTD_at_zero = point(B0x,B0y,distance_mm);

    }
  inline double resolution_at_point(double const& xpoint_mm, double const& ypoint_mm) {
    //laboratory coordinates of the detected point:
    point lab( forward_W * (image_to_detector * point(xpoint_mm,ypoint_mm,distance_mm) + B_XTD_at_zero) );
    //use dot product with direct beam to get the point's two theta; shortcut knowing direct beam is on z-axis
    double theta = std::acos(lab[2]/lab.length()) * 0.5;
    return wavelength_ang / (2* std::sin(theta));
  }



};


} //namespace

} //namespace

#endif //SPOTFILTER_H
