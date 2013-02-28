/*
 * detector.cc
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
#include <string>
#include <iostream>
#include <sstream>
#include <boost_adaptbx/std_pair_conversion.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/multi_panel_detector.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  
  template <typename T>
  std::ostream& operator<<(std::ostream &os, vec2 <T> a) {
    return os << a.const_ref();
  }

  template <typename T>
  std::ostream& operator<<(std::ostream &os, vec3 <T> a) {
    return os << a.const_ref();
  }

  std::string detector_to_string(const Detector &detector)
  {
    std::stringstream os;
    os << "Detector:\n";
    os << "    type:          " << detector.get_type() << "\n";
    os << "    fast axis:     " << detector.get_fast_axis() << "\n";
    os << "    slow axis:     " << detector.get_slow_axis() << "\n";
    os << "    normal:        " << detector.get_normal() << "\n";
    os << "    origin:        " << detector.get_origin() << "\n";
    os << "    pixel size:    " << detector.get_pixel_size() << "\n";
    os << "    image size:    " << detector.get_image_size() << "\n";
    os << "    trusted range: " << detector.get_trusted_range() << "\n";
    return os.str();
  }

  std::string multi_panel_detector_to_string(
      const MultiPanelDetector &detector)
  {
    std::string str;
    str += "MultiPanelDetector:\n";
    for (std::size_t i = 0; i < detector.num_panels(); ++i) {
      str += "\n";
      str += detector_to_string(detector[i]);
    }
    return str;
  }

  void multi_panel_detector_set_item(MultiPanelDetector &d, 
      std::size_t i, const Detector &v) {
    d[i] = v;
  }

  void multi_panel_detector_del_item(MultiPanelDetector &d, 
      std::size_t i) {
    d.remove_panel(i);
  }

  Detector& multi_panel_detector_get_item(
      MultiPanelDetector &d, std::size_t i) {
    return d[i];
  }

  void export_detector() 
  {
    // Export a flex array - should probably move somewhere else
    scitbx::af::boost_python::flex_wrapper <int4>::plain("flex_int4");

    // Export the DetectorBase class
    class_ <DetectorBase> ("DetectorBase");

    // Export the Detector class
    class_ <Detector, bases <DetectorBase> > ("Detector")
      .def(init <std::string,
                 vec3 <double>,
                 vec3 <double>,
                 vec3 <double>,
                 vec2 <double>,
                 vec2 <std::size_t>,
                vec2 <double> > ((                 
          arg("type"),
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin"),
          arg("pixel_size"),
          arg("image_size"),
           arg("trusted_range"))))
      .add_property("type",
        &Detector::get_type,
        &Detector::set_type)    
      .add_property("fast_axis",
        &Detector::get_fast_axis,
        &Detector::set_fast_axis)
      .add_property("slow_axis",
        &Detector::get_slow_axis,
        &Detector::set_slow_axis)
      .add_property("normal",
        &Detector::get_normal)
      .add_property("origin",
        &Detector::get_origin,
        &Detector::set_origin)
      .add_property("pixel_size",
        &Detector::get_pixel_size,
        &Detector::set_pixel_size)
      .add_property("image_size",
        &Detector::get_image_size,
        &Detector::set_image_size)
      .add_property("trusted_range",
        &Detector::get_trusted_range,
        &Detector::set_trusted_range)
      .add_property("mask",
        &Detector::get_mask,
        &Detector::set_mask)
      .def("get_d_matrix",
        &Detector::get_d_matrix)
      .def("get_D_matrix",
        &Detector::get_D_matrix)
      .def("set_d_matrix",
        &Detector::set_d_matrix)
      .def("set_D_matrix",
        &Detector::set_D_matrix)
      .def("add_mask",
        &Detector::add_mask)
      .def("get_pixel_lab_coord",
        &Detector::get_pixel_lab_coord<vec2<double> >)
      .def("get_image_size_mm",
        &Detector::get_image_size_mm)
      .def("is_value_in_trusted_range",
        &Detector::is_value_in_trusted_range)
      .def("is_coord_valid",
        &Detector::is_coord_valid)
      .def("is_coord_valid_mm",
        &Detector::is_coord_valid_mm)
      .def("get_ray_intersection",
        &Detector::get_ray_intersection)
      .def("get_distance",
        &Detector::get_distance)
      .def("get_beam_centre",
        &Detector::get_beam_centre)
      .def("get_beam_centre_lab",
        &Detector::get_beam_centre_lab)
      .def("get_resolution_at_pixel",
        &Detector::get_resolution_at_pixel)
      .def("get_max_resolution_at_corners",
        &Detector::get_max_resolution_at_corners)
      .def("get_max_resolution_elipse",
        &Detector::get_max_resolution_elipse)
      .def("millimeter_to_pixel",
        &Detector::millimeter_to_pixel<vec2<double> >)
      .def("pixel_to_millimeter",
        &Detector::pixel_to_millimeter<vec2<double> >)
      .def("__eq__", &Detector::operator==)
      .def("__ne__", &Detector::operator!=)
      .def("__str__", &detector_to_string);

    // Register std::pair conversion for MultiPanelDetector coordinate type 
    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();

    // Export a MultiPanelDetector class
    class_ <MultiPanelDetector, 
            bases <DetectorBase> > ("MultiPanelDetector")
      .def(init <std::string> ((
          arg("type"))))
      .def("add_panel",
        &MultiPanelDetector::add_panel, (
          arg("panel")))
      .def("num_panels",
        &MultiPanelDetector::num_panels)
      .def("get_d_matrices",
        &MultiPanelDetector::get_d_matrices)
      .def("get_D_matrices",
        &MultiPanelDetector::get_D_matrices)
      .def("is_value_in_trusted_range",
        &MultiPanelDetector::is_value_in_trusted_range)
      .def("is_coord_valid",
        &MultiPanelDetector::is_coord_valid)
      .def("millimeter_to_pixel",
        &MultiPanelDetector::millimeter_to_pixel)
      .def("pixel_to_millimeter",
        &MultiPanelDetector::pixel_to_millimeter)   
//      .def("do_panels_intersect",
//        &MultiPanelDetector::do_panels_intersect)
      .def("__len__", 
        &MultiPanelDetector::num_panels)
      .def("__setitem__", 
        &multi_panel_detector_set_item)
      .def("__delitem__", 
        &multi_panel_detector_del_item)
      .def("__getitem__", 
        &multi_panel_detector_get_item, 
        return_internal_reference <> ())
      .def("__iter__", 
        iterator <
          MultiPanelDetector, 
          return_internal_reference<> >())
      .def("__eq__", &MultiPanelDetector::operator==)
      .def("__ne__", &MultiPanelDetector::operator!=)
      .def("__str__", &multi_panel_detector_to_string);
  }

}}} // namespace dials::model::boost_python
