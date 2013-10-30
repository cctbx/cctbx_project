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

  scitbx::af::shared<vec3<double> >
  get_lab_coord_multiple(const Panel &panel,
                scitbx::af::flex<vec2<double> >::type const& xy) {
    scitbx::af::shared<vec3<double> > result((scitbx::af::reserve(xy.size())));
    for(std::size_t i = 0; i < xy.size(); i++) {
      result.push_back(panel.get_lab_coord(xy[i]));
    }
    return result;
  }

  void export_panel() 
  {
    using namespace boost::python;

    class_<PanelFrame>("PanelFrame")
      .def("set_local_frame", 
        &PanelFrame::set_local_frame, (
          arg("fast_axis"), 
          arg("slow_axis"), 
          arg("origin")))
      .def("set_parent_frame",
        &PanelFrame::set_parent_frame, (
          arg("fast_axis"), 
          arg("slow_axis"), 
          arg("origin")))
      .def("get_local_d_matrix", 
        &PanelFrame::get_local_d_matrix)
      .def("get_parent_d_matrix", 
        &PanelFrame::get_parent_d_matrix)
      .def("get_d_matrix", 
        &PanelFrame::get_d_matrix)
      .def("get_D_matrix", 
        &PanelFrame::get_D_matrix)
      .def("get_origin", 
        &PanelFrame::get_origin)
      .def("get_fast_axis", 
        &PanelFrame::get_fast_axis)
      .def("get_slow_axis", 
        &PanelFrame::get_slow_axis)
      .def("get_normal", 
        &PanelFrame::get_normal)
      .def("get_normal_origin", 
        &PanelFrame::get_normal_origin)
      .def("get_distance",
        &PanelFrame::get_distance)
      .def("get_beam_centre", 
        &PanelFrame::get_beam_centre, (
          arg("s0")))
      .def("get_beam_centre_lab", 
        &PanelFrame::get_beam_centre_lab, (
          arg("s0")))
      .def("get_lab_coord", 
        &PanelFrame::get_lab_coord, (
          arg("xy")))
      .def("get_lab_coord", 
        &get_lab_coord_multiple, (
          arg("xy")))
      .def("get_ray_intersection", 
        &PanelFrame::get_ray_intersection, (
          arg("s1")))
      .def("get_bidirectional_ray_intersection", 
        &PanelFrame::get_bidirectional_ray_intersection, (
          arg("s1")))
      .def("__eq__", &PanelFrame::operator==)
      .def("__ne__", &PanelFrame::operator!=);
          
    class_<PanelBase, bases<PanelFrame> >("PanelBase")
      .def("get_name", &PanelBase::get_name)
      .def("set_name", &PanelBase::set_name)
      .def("get_type", &PanelBase::get_type)
      .def("set_type", &PanelBase::set_type)
      .def("__eq__", &PanelBase::operator==)
      .def("__ne__", &PanelBase::operator!=);
      
    class_<Panel, bases<PanelBase> >("Panel")
      .def(init<std::string,
                std::string,
                vec3 <double>,
                vec3 <double>,
                vec3 <double>,
                vec2 <double>,
                vec2 <std::size_t>,
                vec2 <double> >((
        arg("type"),
        arg("name"),
        arg("fast_axis"),
        arg("slow_axis"),
        arg("origin"),
        arg("pixel_size"),
        arg("image_size"),
        arg("trusted_range"))))
      .def(init<std::string,
                std::string,
                vec3 <double>,
                vec3 <double>,
                vec3 <double>,
                vec2 <double>,
                vec2 <std::size_t>,
                vec2 <double>,
                shared_ptr<PxMmStrategy> >((
        arg("type"),
        arg("name"),
        arg("fast_axis"),
        arg("slow_axis"),
        arg("origin"),
        arg("pixel_size"),
        arg("image_size"),
        arg("trusted_range"),
        arg("px_mm"))))
      .def("get_pixel_size", &Panel::get_pixel_size)
      .def("set_pixel_size", &Panel::set_pixel_size)
      .def("get_image_size", &Panel::get_image_size)
      .def("set_image_size", &Panel::set_image_size)
      .def("get_trusted_range", &Panel::get_trusted_range)
      .def("set_trusted_range", &Panel::set_trusted_range)
      .def("get_image_size_mm", &Panel::get_image_size_mm)
      .def("get_px_mm_strategy", &Panel::get_px_mm_strategy)
      .def("set_px_mm_strategy", &Panel::set_px_mm_strategy)
      .def("set_mask", &Panel::set_mask)
      .def("get_mask", &Panel::get_mask)
      .def("add_mask", &Panel::add_mask)
      .def("is_value_in_trusted_range", &Panel::is_value_in_trusted_range)
      .def("is_coord_valid", &Panel::is_coord_valid)
      .def("is_coord_valid_mm", &Panel::is_coord_valid_mm)
      .def("get_normal_origin_px", &Panel::get_normal_origin_px)
      .def("get_beam_centre_px", &Panel::get_beam_centre_px)
      .def("get_pixel_lab_coord", &Panel::get_pixel_lab_coord)
      .def("get_ray_intersection_px", &Panel::get_ray_intersection_px)
      .def("get_bidirectional_ray_intersection_px", 
        &Panel::get_bidirectional_ray_intersection_px)
      .def("millimeter_to_pixel", &Panel::millimeter_to_pixel)
      .def("pixel_to_millimeter", &Panel::pixel_to_millimeter)
      .def("get_resolution_at_pixel", &Panel::get_resolution_at_pixel)
      .def("get_max_resolution_at_corners", 
        &Panel::get_max_resolution_at_corners)
      .def("get_max_resolution_ellipse", &Panel::get_max_resolution_ellipse)
      .def("__eq__", &Panel::operator==)
      .def("__ne__", &Panel::operator!=);
  }

}}} // namespace dxtbx::model::boost_python
