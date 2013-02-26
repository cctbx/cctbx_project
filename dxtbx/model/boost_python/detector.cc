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

  std::string flat_panel_detector_to_string(const FlatPanelDetector &detector)
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

  std::string multi_flat_panel_detector_to_string(
      const MultiFlatPanelDetector &detector)
  {
    std::string str;
    str += "MultiFlatPanelDetector:\n";
    for (std::size_t i = 0; i < detector.num_panels(); ++i) {
      str += "\n";
      str += flat_panel_detector_to_string(detector[i]);
    }
    return str;
  }

  void multi_flat_panel_detector_set_item(MultiFlatPanelDetector &d, 
      std::size_t i, const FlatPanelDetector &v) {
    d[i] = v;
  }

  void multi_flat_panel_detector_del_item(MultiFlatPanelDetector &d, 
      std::size_t i) {
    d.remove_panel(i);
  }

  FlatPanelDetector& multi_flat_panel_detector_get_item(
      MultiFlatPanelDetector &d, std::size_t i) {
    return d[i];
  }

  void export_detector() 
  {
    // Export a flex array - should probably move somewhere else
    scitbx::af::boost_python::flex_wrapper <int4>::plain("flex_int4");

    // Export the DetectorBase class
    class_ <DetectorBase> ("DetectorBase");

    // Export the FlatPanelDetector class
    class_ <FlatPanelDetector, bases <DetectorBase> > ("FlatPanelDetector")
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
        &FlatPanelDetector::get_type,
        &FlatPanelDetector::set_type)    
      .add_property("fast_axis",
        &FlatPanelDetector::get_fast_axis,
        &FlatPanelDetector::set_fast_axis)
      .add_property("slow_axis",
        &FlatPanelDetector::get_slow_axis,
        &FlatPanelDetector::set_slow_axis)
      .add_property("normal",
        &FlatPanelDetector::get_normal)
      .add_property("origin",
        &FlatPanelDetector::get_origin,
        &FlatPanelDetector::set_origin)
      .add_property("pixel_size",
        &FlatPanelDetector::get_pixel_size,
        &FlatPanelDetector::set_pixel_size)
      .add_property("image_size",
        &FlatPanelDetector::get_image_size,
        &FlatPanelDetector::set_image_size)
      .add_property("trusted_range",
        &FlatPanelDetector::get_trusted_range,
        &FlatPanelDetector::set_trusted_range)
      .add_property("d_matrix",
        &FlatPanelDetector::get_d_matrix,
        &FlatPanelDetector::set_d_matrix)
      .add_property("inverse_d_matrix",
        &FlatPanelDetector::get_D_matrix,
        &FlatPanelDetector::set_D_matrix)
      .add_property("mask",
        &FlatPanelDetector::get_mask,
        &FlatPanelDetector::set_mask)
      .def("add_mask",
        &FlatPanelDetector::add_mask)
      .def("get_pixel_lab_coord",
        &FlatPanelDetector::get_pixel_lab_coord<vec2<double> >)
      .def("get_image_rectangle",
        &FlatPanelDetector::get_image_rectangle)
      .def("get_image_size_mm",
        &FlatPanelDetector::get_image_size_mm)
      .def("is_value_in_trusted_range",
        &FlatPanelDetector::is_value_in_trusted_range)
      .def("is_coord_valid",
        &FlatPanelDetector::is_coord_valid)
      .def("millimeter_to_pixel",
        &FlatPanelDetector::millimeter_to_pixel<vec2<double> >)
      .def("pixel_to_millimeter",
        &FlatPanelDetector::pixel_to_millimeter<vec2<double> >)
      .def("__eq__", &FlatPanelDetector::operator==)
      .def("__ne__", &FlatPanelDetector::operator!=)
      .def("__str__", &flat_panel_detector_to_string);

    // Register std::pair conversion for MultiFlatPanelDetector coordinate type 
    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();

    // Export a MultiFlatPanelDetector class
    class_ <MultiFlatPanelDetector, 
            bases <DetectorBase> > ("MultiFlatPanelDetector")
      .def(init <std::string> ((
          arg("type"))))
      .def("add_panel",
        &MultiFlatPanelDetector::add_panel, (
          arg("panel")))
      .def("num_panels",
        &MultiFlatPanelDetector::num_panels)
      .def("is_value_in_trusted_range",
        &MultiFlatPanelDetector::is_value_in_trusted_range)
      .def("is_coord_valid",
        &MultiFlatPanelDetector::is_coord_valid)
      .def("millimeter_to_pixel",
        &MultiFlatPanelDetector::millimeter_to_pixel)
      .def("pixel_to_millimeter",
        &MultiFlatPanelDetector::pixel_to_millimeter)   
      .def("do_panels_intersect",
        &MultiFlatPanelDetector::do_panels_intersect)
      .def("__len__", 
        &MultiFlatPanelDetector::num_panels)
      .def("__setitem__", 
        &multi_flat_panel_detector_set_item)
      .def("__delitem__", 
        &multi_flat_panel_detector_del_item)
      .def("__getitem__", 
        &multi_flat_panel_detector_get_item, 
        return_internal_reference <> ())
      .def("__iter__", 
        iterator <
          MultiFlatPanelDetector, 
          return_internal_reference<> >())
      .def("__eq__", &MultiFlatPanelDetector::operator==)
      .def("__ne__", &MultiFlatPanelDetector::operator!=)
      .def("__str__", &multi_flat_panel_detector_to_string);
  }

}}} // namespace dials::model::boost_python
