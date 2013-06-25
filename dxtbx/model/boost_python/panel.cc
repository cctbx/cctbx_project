/*
 * panel.cc
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
#include <dxtbx/model/panel.h>

namespace dxtbx { namespace model { namespace boost_python {

  
  std::string panel_to_string(const Panel &panel) {
    std::stringstream os;
    os << panel;
    return os.str();
  }

  struct PanelPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Panel &obj) {
      return boost::python::make_tuple(
        obj.get_type(),
        obj.get_fast_axis(),
        obj.get_slow_axis(),
        obj.get_origin(),
        obj.get_pixel_size(),
        obj.get_image_size(),
        obj.get_trusted_range());
    }
  };

  void export_panel() 
  {
    using namespace boost::python;
    // Export a flex array - should probably move somewhere else
    scitbx::af::boost_python::flex_wrapper <int4>::plain("flex_int4");

    // Export the Panel class
    class_ <Panel> ("Panel")
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
      .def(init <std::string,
                 vec3 <double>,
                 vec3 <double>,
                 vec3 <double>,
                 vec2 <double>,
                 vec2 <std::size_t>,
                 vec2 <double>,
                 shared_ptr<PxMmStrategy> > ((                 
          arg("type"),
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin"),
          arg("pixel_size"),
          arg("image_size"),
          arg("trusted_range"),
          arg("px_mm"))))          
      .def("get_type",
        &Panel::get_type)
      .def("set_type",
        &Panel::set_type)    
      .def("get_fast_axis",
        &Panel::get_fast_axis)
      .def("get_slow_axis",
        &Panel::get_slow_axis)
      .def("get_origin",
        &Panel::get_origin)
      .def("set_frame",
        &Panel::set_frame, (
          arg("fast_axis"), 
          arg("slow_axis"),
          arg("origin")))
      .def("get_normal",
        &Panel::get_normal)
      .def("get_pixel_size",
        &Panel::get_pixel_size)
      .def("set_pixel_size",
        &Panel::set_pixel_size)
      .def("get_image_size",
        &Panel::get_image_size)
      .def("set_image_size",
        &Panel::set_image_size)
      .def("get_trusted_range",
        &Panel::get_trusted_range)
      .def("set_trusted_range",
        &Panel::set_trusted_range)
      .def("get_mask",
        &Panel::get_mask)
      .def("set_mask",
        &Panel::set_mask)
      .def("get_d_matrix",
        &Panel::get_d_matrix)
      .def("get_D_matrix",
        &Panel::get_D_matrix)
      .def("add_mask",
        &Panel::add_mask)
      .def("get_lab_coord",
        &Panel::get_lab_coord)
      .def("get_pixel_lab_coord",
        &Panel::get_pixel_lab_coord)
      .def("get_image_size_mm",
        &Panel::get_image_size_mm)
      .def("is_value_in_trusted_range",
        &Panel::is_value_in_trusted_range)
      .def("is_coord_valid",
        &Panel::is_coord_valid)
      .def("is_coord_valid_mm",
        &Panel::is_coord_valid_mm)
      .def("get_ray_intersection",
        &Panel::get_ray_intersection)
      .def("get_ray_intersection_px",
        &Panel::get_ray_intersection_px)
      .def("get_distance",
        &Panel::get_distance)
      .def("get_beam_centre",
        &Panel::get_beam_centre, (
          arg("s0")))
      .def("get_beam_centre_lab",
        &Panel::get_beam_centre_lab, (
          arg("s0")))
      .def("get_resolution_at_pixel",
        &Panel::get_resolution_at_pixel, (
          arg("s0"),
          arg("wavelength"),
          arg("xy")))
      .def("get_max_resolution_at_corners",
        &Panel::get_max_resolution_at_corners, (
          arg("s0"),
          arg("wavelength")))
      .def("get_max_resolution_elipse",
        &Panel::get_max_resolution_elipse, (
          arg("s0"),
          arg("wavelength")))
      .def("millimeter_to_pixel",
        &Panel::millimeter_to_pixel)
      .def("pixel_to_millimeter",
        &Panel::pixel_to_millimeter)
      .def("__eq__", &Panel::operator==)
      .def("__ne__", &Panel::operator!=)
      .def("__str__", &panel_to_string)
      .def_pickle(PanelPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
