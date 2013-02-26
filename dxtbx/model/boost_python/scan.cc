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
#include <dxtbx/model/scan.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  std::string scan_to_string(const ScanData &scan) {
    boost::format fmt(
      "Scan:\n"
      "    image range:       (%1%, %2%)\n"
      "    oscillation:       (%3%, %4%)\n"
      "    exposure time:     %5%");
        
    fmt % scan.get_image_range()[0];
    fmt % scan.get_image_range()[1];
    fmt % scan.get_oscillation()[0];
    fmt % scan.get_oscillation()[1];
    fmt % scan.get_exposure_time();
    return fmt.str();
  }

  static ScanData* make_scan(vec2 <int> image_range, vec2 <double> oscillation,
      double exposure_time, const flex_double &epochs, bool deg) {
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
  
  void export_scan()
  {
    // Export ScanBase
    class_ <ScanBase> ("ScanBase");

    // Export Scan : ScanBase
    class_ <ScanData, bases <ScanBase> > ("ScanData")
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
          arg("epochs"),          
          arg("deg"))))
      .add_property("image_range",  
        &ScanData::get_image_range,
        &ScanData::set_image_range)
      .add_property("oscillation",  
        &ScanData::get_oscillation,
        &ScanData::set_oscillation)
      .add_property("exposure_time",
        &ScanData::get_exposure_time,
        &ScanData::set_exposure_time)
      .add_property("epochs",
        &ScanData::get_epochs,
        &ScanData::set_epochs)
      .add_property("oscillation_range",
        &ScanData::get_oscillation_range)
      .add_property("num_images",
        &ScanData::get_num_images)
      .def("get_image_oscillation",
        &ScanData::get_image_oscillation, (
          arg("index")))
      .def("get_image_epoch",
        &ScanData::get_image_epoch, (
          arg("index")))
      .def("__eq__", &ScanData::operator==)
      .def("__nq__", &ScanData::operator!=)
      .def("__len__", &ScanData::get_num_images)
      .def("__str__", &scan_to_string);
  }

}}} // namespace = dxtbx::model::boost_python
