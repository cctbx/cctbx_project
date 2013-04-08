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
#include <boost/python/slice.hpp>
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <dxtbx/model/scan.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  static
  vec2<double> rad_as_deg(vec2<double> angles) {
    angles[0] = rad_as_deg(angles[0]);
    angles[1] = rad_as_deg(angles[1]);
    return angles;
  }

  std::string scan_to_string(const Scan &scan) {
    std::stringstream ss;
    ss << scan;
    return ss.str();
  }

  struct ScanPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Scan &obj) {
      return boost::python::make_tuple(
        obj.get_image_range(),
        rad_as_deg(obj.get_oscillation()),
        obj.get_exposure_time(),
        obj.get_epochs());
    }
  };

  static Scan* make_scan(vec2 <int> image_range, vec2 <double> oscillation,
      double exposure_time, bool deg) {
    Scan *scan = NULL;
    if (deg) {
      scan = new Scan(image_range, 
        vec2 <double> (
          deg_as_rad(oscillation[0]), 
          deg_as_rad(oscillation[1])), 
        exposure_time);
    } else {
      scan = new Scan(image_range, oscillation, exposure_time);
    }
    return scan;
  }

  static Scan* make_scan_w_epoch(vec2 <int> image_range, 
      vec2 <double> oscillation, double exposure_time, 
      const flex_double &epochs, bool deg) {
    Scan *scan = NULL;
    if (deg) {
      scan = new Scan(image_range, 
        vec2 <double> (
          deg_as_rad(oscillation[0]), 
          deg_as_rad(oscillation[1])), 
        exposure_time, epochs);
    } else {
      scan = new Scan(image_range, oscillation, exposure_time, epochs);
    }
    return scan;
  }
  
 
  static
  vec2<double> get_oscillation_range(const Scan &scan, bool deg) {
    vec2<double> range = scan.get_oscillation_range();
    return deg ? rad_as_deg(range) : range;
  }

  static  
  vec2<double> get_oscillation(const Scan &scan, bool deg) {
    vec2<double> oscillation = scan.get_oscillation();
    return deg ? rad_as_deg(oscillation) : oscillation;
  }

  static
  void set_oscillation(Scan &scan, vec2<double> oscillation,
      bool deg) {
    if (deg) {
      oscillation = rad_as_deg(oscillation);
    }
    scan.set_oscillation(oscillation);
  }

   static  
  vec2<double> get_image_oscillation(const Scan &scan, int image, 
      bool deg) {
    vec2<double> oscillation = scan.get_image_oscillation(image);
    return deg ? rad_as_deg(oscillation) : oscillation;
  }

  
  static 
  bool is_angle_valid(const Scan &scan, double angle, bool deg) {
    return scan.is_angle_valid(deg ? deg_as_rad(angle) : angle);
  }

  static 
  double get_angle_from_image_index(const Scan &scan, double index, 
      bool deg) {
    double angle = scan.get_angle_from_image_index(index);
    return deg ? rad_as_deg(angle) : angle;
  }

  static 
  double get_angle_from_array_index(const Scan &scan, double index, 
      bool deg) {
    double angle = scan.get_angle_from_array_index(index);
    return deg ? rad_as_deg(angle) : angle;
  }

  static 
  double get_image_index_from_angle(const Scan &scan, double angle, 
      bool deg) {
    return scan.get_image_index_from_angle(deg ? deg_as_rad(angle) : angle);
  }

  static 
  double get_array_index_from_angle(const Scan &scan, double angle,
      bool deg) {
    return scan.get_array_index_from_angle(
      deg ? deg_as_rad(angle) : angle);
  }

  static 
  flex_double get_image_indices_with_angle(const Scan &scan, double angle, 
      bool deg) {
    return scan.get_image_indices_with_angle(deg ? deg_as_rad(angle) : angle);
  }
  
  static 
  flex_double get_array_indices_with_angle(const Scan &scan, 
      double angle, bool deg) {
    return scan.get_array_indices_with_angle(
      deg ? deg_as_rad(angle) : angle);
  }  
  
  static 
  Scan getitem_single(const Scan &scan, int index) 
  {
    // Check index
    DXTBX_ASSERT(scan.get_image_range()[0] <= index && 
      index <= scan.get_image_range()[1]);

    // Create the new epoch array
    flex_double new_epochs(1);
    new_epochs[0] = scan.get_image_epoch(index);

    // Return scan
    return Scan(vec2<int>(index, index), 
      scan.get_image_oscillation(index),
      scan.get_exposure_time(),
      new_epochs);
  }
  
  static
  Scan getitem_slice(const Scan &scan, const slice index)
  {
    vec2<int> image_range = scan.get_image_range();

    // Ensure no step
    DXTBX_ASSERT(index.step() == object());
    
    // Get start index
    int start = 0, stop = 0;
    if (index.start() == object()) {
      start = image_range[0];
    } else {
      start = extract<int>(index.start());
    }
    
    // Get stop index
    if (index.stop() == object()) {
      stop = image_range[1];
    } else {
      stop = extract<int>(index.stop());
    }

    // Check ranges
    DXTBX_ASSERT(start >= image_range[0]);
    DXTBX_ASSERT(stop <= image_range[1]);
    DXTBX_ASSERT(start <= stop);

    // Create the new epoch array
    flex_double new_epochs(stop - start + 1);
    for (std::size_t i = 0; i < new_epochs.size(); ++i) {
      new_epochs[i] = scan.get_image_epoch(i + start);
    }

    // Create the new scan object
    return Scan(vec2<int>(start, stop), 
      scan.get_image_oscillation(start), 
      scan.get_exposure_time(), new_epochs);
  }  
  
  void export_scan()
  {
    // Export ScanBase
    class_ <ScanBase> ("ScanBase");

    // Export Scan : ScanBase
    class_ <Scan, bases <ScanBase> > ("Scan")
      .def("__init__",
          make_constructor(
          &make_scan, 
          default_call_policies(), (
          arg("image_range"),
          arg("oscillation"),
          arg("exposure_time"),
          arg("deg") = true)))
      .def("__init__",
          make_constructor(
          &make_scan_w_epoch, 
          default_call_policies(), (
          arg("image_range"),
          arg("oscillation"),
          arg("exposure_time"),
          arg("epochs"),          
          arg("deg") = true)))
      .def("get_image_range",  
        &Scan::get_image_range)
      .def("set_image_range",
        &Scan::set_image_range)
      .def("get_array_range",  
        &Scan::get_array_range)
      .def("get_oscillation",  
        &get_oscillation, (
          arg("deg") = true))
      .def("set_oscillation",
        &set_oscillation, (
          arg("deg") = true))
      .def("get_exposure_time",
        &Scan::get_exposure_time)
      .def("set_exposure_time",
        &Scan::set_exposure_time)
      .def("get_epochs",
        &Scan::get_epochs)
      .def("set_epochs",
        &Scan::set_epochs)
      .def("get_num_images",
        &Scan::get_num_images)
      .def("get_image_oscillation",
        &get_image_oscillation, (
          arg("index"),
          arg("deg") = true))
      .def("get_image_epoch",
        &Scan::get_image_epoch, (
          arg("index")))
      .def("get_oscillation_range",
        &get_oscillation_range, (
          arg("deg") = true))          
      .def("is_angle_valid",
        &is_angle_valid, (
          arg("angle"),
          arg("deg") = true))
      .def("is_image_index_valid",
        &Scan::is_image_index_valid, (
          arg("index")))
      .def("is_array_index_valid",
        &Scan::is_array_index_valid, (
          arg("index")))
      .def("get_angle_from_image_index",
        &get_angle_from_image_index, (
          arg("index"),
          arg("deg") = true))
      .def("get_angle_from_array_index",
        &get_angle_from_array_index, (
          arg("index"),
          arg("deg") = true))          
      .def("get_image_index_from_angle",
        &get_image_index_from_angle, (
          arg("angle"),
          arg("deg") = true))
      .def("get_array_index_from_angle",
        &get_array_index_from_angle, (
          arg("angle"),
          arg("deg") = true))
      .def("get_image_indices_with_angle",
        &get_image_indices_with_angle, (
          arg("angle"),
          arg("deg") = true))
      .def("get_array_indices_with_angle",
        &get_array_indices_with_angle, (
          arg("angle"),
          arg("deg") = true))
      .def("__getitem__", &getitem_single)
      .def("__getitem__", &getitem_slice)
      .def("__eq__", &Scan::operator==)
      .def("__nq__", &Scan::operator!=)
      .def("__lt__", &Scan::operator<)
      .def("__le__", &Scan::operator<=)
      .def("__gt__", &Scan::operator>)
      .def("__ge__", &Scan::operator>=)
      .def("__add__", &Scan::operator+)
      .def("__len__", &Scan::get_num_images)
      .def("__str__", &scan_to_string)
      .def_pickle(ScanPickleSuite());
  }

}}} // namespace = dxtbx::model::boost_python
