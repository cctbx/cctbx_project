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

  
  std::string detector_to_string(const Detector &detector) {
    std::stringstream os;
    os << detector;
    return os.str();
  }

  void detector_set_item(Detector &d, 
      std::size_t i, const Panel &v) {
    d[i] = v;
  }

  void detector_del_item(Detector &d, 
      std::size_t i) {
    d.remove_panel(i);
  }

  Panel& detector_get_item(
      Detector &d, std::size_t i) {
    return d[i];
  }

  scitbx::af::shared<vec3<double> >
  get_lab_coord(const Detector &detector,
                scitbx::af::flex<vec2<double> >::type const& xy)
  {
    scitbx::af::shared<vec3<double> > result((scitbx::af::reserve(xy.size())));
    for(std::size_t i=0;i<xy.size();i++) {
      result.push_back(detector.get_lab_coord(xy[i]));
    }
    return result;
  }

  void export_detector() 
  {
    using namespace boost::python;
    // Register std::pair conversion for Detector coordinate type 
    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();

    // Export a panel list type
    scitbx::af::boost_python::flex_wrapper 
      <Panel>::plain("PanelList")
        .enable_pickling(); 

    // Export a Detector base class
    class_ <DetectorBase> ("DetectorBase")
      .def(init<const Panel&>((
        arg("panel"))))
      .def(init<const Detector::panel_list_type&>((
        arg("panel_list"))))
      .def("add_panel",
        &Detector::add_panel, (
          arg("panel")))
      .def("num_panels",
        &Detector::num_panels)
      .def("get_names",
        &Detector::get_names)
      .def("get_d_matrices",
        &Detector::get_d_matrices)
      .def("get_D_matrices",
        &Detector::get_D_matrices)
      .def("get_max_resolution",
        &Detector::get_max_resolution, (
          arg("s0"),
          arg("wavelength")))
      .def("get_ray_intersection",
        &Detector::get_ray_intersection, (
          arg("s1")))
//      .def("do_panels_intersect",
//        &Detector::do_panels_intersect)
      .def("__len__", 
        &Detector::num_panels)
      .def("__setitem__", 
        &detector_set_item)
      .def("__delitem__", 
        &detector_del_item)
      .def("__getitem__", 
        &detector_get_item, 
        return_internal_reference <> ())
      .def("__iter__", 
        iterator <
          Detector, 
          return_internal_reference<> >())
      .def("__eq__", &Detector::operator==)
      .def("__ne__", &Detector::operator!=)
      .def("__str__", &detector_to_string);

    // Export the detector class
    class_ <Detector, bases<DetectorBase> >("Detector")
      .def(init<const Detector&>())
      .def(init<const Panel&>((
        arg("panel"))))
      .def(init<const Detector::panel_list_type&>((
        arg("panel_list"))))
      .def("get_type",
        &Detector::get_type)
      .def("set_type",
        &Detector::set_type)   
      .def("get_name",
        &Detector::get_name)
      .def("set_name",
        &Detector::set_name)   
      .def("get_fast_axis",
        &Detector::get_fast_axis)
      .def("get_slow_axis",
        &Detector::get_slow_axis)
      .def("get_origin",
        &Detector::get_origin)
      .def("set_frame",
        &Detector::set_frame, (
          arg("fast_axis"), 
          arg("slow_axis"),
          arg("origin")))
      .def("get_normal",
        &Detector::get_normal)
      .def("get_pixel_size",
        &Detector::get_pixel_size)
      .def("set_pixel_size",
        &Detector::set_pixel_size)
      .def("get_image_size",
        &Detector::get_image_size)
      .def("set_image_size",
        &Detector::set_image_size)
      .def("get_trusted_range",
        &Detector::get_trusted_range)
      .def("set_trusted_range",
        &Detector::set_trusted_range)
      .def("get_mask",
        &Detector::get_mask)
      .def("set_mask",
        &Detector::set_mask)
      .def("get_d_matrix",
        &Detector::get_d_matrix)
      .def("get_D_matrix",
        &Detector::get_D_matrix)
      .def("add_mask",
        &Detector::add_mask)
      .def("get_lab_coord",
        &Detector::get_lab_coord)
      .def("get_lab_coord",
        get_lab_coord)
      .def("get_pixel_lab_coord",
        &Detector::get_pixel_lab_coord)
      .def("get_image_size_mm",
        &Detector::get_image_size_mm)
      .def("is_value_in_trusted_range",
        &Detector::is_value_in_trusted_range)
      .def("is_coord_valid",
        &Detector::is_coord_valid)
      .def("is_coord_valid_mm",
        &Detector::is_coord_valid_mm)
      .def("get_distance",
        &Detector::get_distance)
      .def("get_beam_centre",
        &Detector::get_beam_centre, (
          arg("s0")))
      .def("get_beam_centre_lab",
        &Detector::get_beam_centre_lab, (
          arg("s0")))
      .def("get_resolution_at_pixel",
        &Detector::get_resolution_at_pixel, (
          arg("s0"),
          arg("wavelength"),
          arg("xy")))
      .def("get_max_resolution_at_corners",
        &Detector::get_max_resolution_at_corners, (
          arg("s0"),
          arg("wavelength")))
      .def("get_max_resolution_elipse",
        &Detector::get_max_resolution_elipse, (
          arg("s0"),
          arg("wavelength")))
      .def("millimeter_to_pixel",
        &Detector::millimeter_to_pixel)
      .def("pixel_to_millimeter",
        &Detector::pixel_to_millimeter)
      .enable_pickling();
  }

}}} // namespace dxtbx::model::boost_python
