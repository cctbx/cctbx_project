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

  double resolution_at_point (
    double xpoint,
    double ypoint,
    double xbeam,
    double ybeam,
    double distance,
    double wavelength,
    double twotheta,
    double pixel_size);

} //namespace

} //namespace

#endif //SPOTFILTER_H
