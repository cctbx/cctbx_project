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
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <dxtbx/model/scan.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  std::string scan_to_string(const ScanData &scan) {
    std::stringstream ss;
    ss << scan;
    return ss.str();
  }

  static ScanData* make_scan(vec2 <int> image_range, vec2 <double> oscillation,
      double exposure_time, bool deg) {
    ScanData *scan = NULL;
    if (deg) {
      scan = new ScanData(image_range, 
        vec2 <double> (
          deg_as_rad(oscillation[0]), 
          deg_as_rad(oscillation[1])), 
        exposure_time);
    } else {
      scan = new ScanData(image_range, oscillation, exposure_time);
    }
    return scan;
  }

  static ScanData* make_scan_w_epoch(vec2 <int> image_range, 
      vec2 <double> oscillation, double exposure_time, 
      const flex_double &epochs, bool deg) {
    ScanData *scan = NULL;
    if (deg) {
      scan = new ScanData(image_range, 
        vec2 <double> (
          deg_as_rad(oscillation[0]), 
          deg_as_rad(oscillation[1])), 
        exposure_time, epochs);
    } else {
      scan = new ScanData(image_range, oscillation, exposure_time, epochs);
    }
    return scan;
  }
  
  static
  vec2<double> rad_as_deg(vec2<double> angles) {
    angles[0] = rad_as_deg(angles[0]);
    angles[1] = rad_as_deg(angles[1]);
    return angles;
  }
  
  static
  vec2<double> get_oscillation_range(const ScanData &scan, bool deg) {
    vec2<double> range = scan.get_oscillation_range();
    return deg ? rad_as_deg(range) : range;
  }
  
  
  static 
  bool is_angle_valid(const ScanData &scan, double angle, bool deg) {
    return scan.is_angle_valid(deg ? deg_as_rad(angle) : angle);
  }

  static 
  double get_angle_from_frame(const ScanData &scan, double frame, bool deg) {
    double angle = scan.get_angle_from_frame(frame);
    return deg ? rad_as_deg(angle) : angle;
  }

  static 
  double get_frame_from_angle(const ScanData &scan, double angle, bool deg) {
    return scan.get_frame_from_angle(deg ? deg_as_rad(angle) : angle);
  }

  static 
  flex_double get_frames_with_angle(const ScanData &scan, double angle, 
      bool deg) {
    return scan.get_frames_with_angle(deg ? deg_as_rad(angle) : angle);
  }
  
  void export_scan()
  {
    // Export ScanBase
    class_ <ScanBase> ("ScanBase");

    // Export Scan : ScanBase
    class_ <ScanData, bases <ScanBase> > ("ScanData")
      .def(init <vec2 <int>, vec2 <double>, double> ((
          arg("image_range"), 
          arg("oscillation"),
          arg("exposure_time"))))
      .def(init <vec2 <int>, vec2 <double>, double, const flex_double &> ((
          arg("image_range"), 
          arg("oscillation"),
          arg("exposure_time"),
          arg("epochs"))))
      .def("__init__",
          make_constructor(
          &make_scan, 
          default_call_policies(), (
          arg("image_range"),
          arg("oscillation"),
          arg("exposure_time"),
          arg("deg"))))
      .def("__init__",
          make_constructor(
          &make_scan_w_epoch, 
          default_call_policies(), (
          arg("image_range"),
          arg("oscillation"),
          arg("exposure_time"),
          arg("epochs"),          
          arg("deg"))))
      .def("get_image_range",  
        &ScanData::get_image_range)
      .def("set_image_range",
        &ScanData::set_image_range)
      .def("get_oscillation",  
        &ScanData::get_oscillation)
      .def("set_oscillation",
        &ScanData::set_oscillation)
      .def("get_exposure_time",
        &ScanData::get_exposure_time)
      .def("set_exposure_time",
        &ScanData::set_exposure_time)
      .def("get_epochs",
        &ScanData::get_epochs)
      .def("set_epochs",
        &ScanData::set_epochs)
      .def("get_num_images",
        &ScanData::get_num_images)
      .def("get_image_oscillation",
        &ScanData::get_image_oscillation, (
          arg("index")))
      .def("get_image_epoch",
        &ScanData::get_image_epoch, (
          arg("index")))
      .def("get_oscillation_range",
        &get_oscillation_range, (
          arg("deg") = false))          
      .def("is_angle_valid",
        &is_angle_valid, (
          arg("angle"),
          arg("deg") = false))
      .def("is_frame_valid",
        &ScanData::is_frame_valid, (
          arg("frame")))
      .def("get_angle_from_frame",
        &get_angle_from_frame, (
          arg("frame"),
          arg("deg") = false))
      .def("get_frame_from_angle",
        &get_frame_from_angle, (
          arg("angle"),
          arg("deg") = false))
      .def("get_frames_with_angle",
        &get_frames_with_angle, (
          arg("angle"),
          arg("deg") = false))
      .def("__eq__", &ScanData::operator==)
      .def("__nq__", &ScanData::operator!=)
      .def("__len__", &ScanData::get_num_images)
      .def("__str__", &scan_to_string);
  }

}}} // namespace = dxtbx::model::boost_python
