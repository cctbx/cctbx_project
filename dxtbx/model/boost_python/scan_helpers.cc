/*
 * scan.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/format.hpp>
#include <string>
#include <scitbx/constants.h>
#include <dxtbx/model/scan_helpers.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  static
  vec2<double> deg_as_rad(vec2<double> angles) {
    angles[0] = deg_as_rad(angles[0]);
    angles[1] = deg_as_rad(angles[1]);
    return angles;
  }

  static
  vec2<double> rad_as_deg(vec2<double> angles) {
    angles[0] = rad_as_deg(angles[0]);
    angles[1] = rad_as_deg(angles[1]);
    return angles;
  }
  
  static 
  bool is_angle_in_range_wrapper(vec2<double> range, double angle, 
      bool deg) {
    return is_angle_in_range(
      deg ? deg_as_rad(range) : range,
      deg ? deg_as_rad(angle) : angle);
  }

  static vec2 <double> 
  get_range_of_mod2pi_angles_wrapper(vec2<double> range, double angle, 
      bool deg) {
    return rad_as_deg(get_range_of_mod2pi_angles(
      deg ? deg_as_rad(range) : range,
      deg ? deg_as_rad(angle) : angle));
  }

  static flex_double 
  get_mod2pi_angles_in_range_wrapper(vec2 <double> range, double angle, 
      bool deg) {
    flex_double result = get_mod2pi_angles_in_range(
      deg ? deg_as_rad(range) : range,
      deg ? deg_as_rad(angle) : angle);    
    if (deg) {
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = rad_as_deg(result[i]);
      }
    }
    return result;
  }
  
  void export_scan_helpers()
  {
    def("is_angle_in_range", &is_angle_in_range_wrapper, (
        arg("range"), arg("angle"), arg("deg") = false));

    def("get_range_of_mod2pi_angles", &get_range_of_mod2pi_angles_wrapper, (
        arg("range"), arg("angle"), arg("deg") = false));

    def("get_mod2pi_angles_in_range", &get_mod2pi_angles_in_range_wrapper, (
        arg("range"), arg("angle"), arg("deg") = false));
  }

}}} // namespace = dxtbx::model::boost_python
