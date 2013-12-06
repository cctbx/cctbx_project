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
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  std::string panel_to_string(const Panel &panel) {
    std::stringstream ss;
    ss << panel;
    return ss.str();
  }

  static
  scitbx::af::shared<vec3<double> >
  get_lab_coord_multiple(const Panel &panel,
                scitbx::af::flex<vec2<double> >::type const& xy) {
    scitbx::af::shared<vec3<double> > result((scitbx::af::reserve(xy.size())));
    for(std::size_t i = 0; i < xy.size(); i++) {
      result.push_back(panel.get_lab_coord(xy[i]));
    }
    return result;
  }

  static
  Panel panel_deepcopy(const Panel &panel, boost::python::object dict) {
    return Panel(panel);
  }

  static
  bool panel_is(const Panel *lhs, const Panel *rhs) {
    return lhs == rhs;
  }

  struct PanelBasePickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(boost::python::object obj) {
      using namespace boost::python;

      unsigned int version = 1;
      const PanelBase &p = extract<const PanelBase&>(obj)();
      return boost::python::make_tuple(
        version,
        obj.attr("__dict__"),
        p.get_name(),
        p.get_type(),
        p.get_local_fast_axis(),
        p.get_local_slow_axis(),
        p.get_local_origin(),
        p.get_parent_fast_axis(),
        p.get_parent_slow_axis(),
        p.get_parent_origin());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      using namespace boost::python;

      DXTBX_ASSERT(len(state) == 10);
      unsigned int version = extract<unsigned int>(state[0]);
      DXTBX_ASSERT(version == 1);
      extract<dict>(obj.attr("__dict__"))().update(state[1]);
      PanelBase &p = extract<PanelBase&>(obj)();
      p.set_name(extract<std::string>(state[2]));
      p.set_type(extract<std::string>(state[3]));
      p.set_local_frame(
        extract< vec3<double> >(state[4]),
        extract< vec3<double> >(state[5]),
        extract< vec3<double> >(state[6]));
      p.set_parent_frame(
        extract< vec3<double> >(state[7]),
        extract< vec3<double> >(state[8]),
        extract< vec3<double> >(state[9]));
    }

    static bool getstate_manages_dict() { return true; }
  };

  struct PanelPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(boost::python::object obj) {
      using namespace boost::python;

      unsigned int version = 1;
      const Panel &p = extract<const Panel&>(obj)();
      return boost::python::make_tuple(
        version,
        obj.attr("__dict__"),
        p.get_name(),
        p.get_type(),
        p.get_local_fast_axis(),
        p.get_local_slow_axis(),
        p.get_local_origin(),
        p.get_parent_fast_axis(),
        p.get_parent_slow_axis(),
        p.get_parent_origin(),
        p.get_pixel_size(),
        p.get_image_size(),
        p.get_trusted_range(),
        p.get_px_mm_strategy());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      using namespace boost::python;

      DXTBX_ASSERT(len(state) == 14);
      unsigned int version = extract<unsigned int>(state[0]);
      DXTBX_ASSERT(version == 1);
      extract<dict>(obj.attr("__dict__"))().update(state[1]);
      Panel &p = extract<Panel&>(obj)();
      p.set_name(extract<std::string>(state[2]));
      p.set_type(extract<std::string>(state[3]));
      p.set_local_frame(
        extract< vec3<double> >(state[4]),
        extract< vec3<double> >(state[5]),
        extract< vec3<double> >(state[6]));
      p.set_parent_frame(
        extract< vec3<double> >(state[7]),
        extract< vec3<double> >(state[8]),
        extract< vec3<double> >(state[9]));
      p.set_pixel_size(extract< vec2<double> >(state[10]));
      p.set_image_size(extract< vec2<std::size_t> >(state[11]));
      p.set_trusted_range(extract< vec2<double> >(state[12]));
      p.set_px_mm_strategy(extract< shared_ptr<PxMmStrategy> >(state[13]));
    }

    static bool getstate_manages_dict() { return true; }
  };

  template <>
  boost::python::dict to_dict<PanelBase>(const PanelBase &obj) {
    boost::python::dict result;
    result["name"] = obj.get_name();
    result["type"] = obj.get_type();
    result["fast_axis"] = obj.get_local_fast_axis();
    result["slow_axis"] = obj.get_local_slow_axis();
    result["origin"] = obj.get_local_origin();
    return result;
  }

  template <>
  PanelBase* from_dict<PanelBase>(boost::python::dict obj) {
    PanelBase *result = new PanelBase();
    if (obj.has_key("name")) {
      result->set_name(boost::python::extract<std::string>(obj["name"]));
    }
    if (obj.has_key("type")) {
      result->set_type(boost::python::extract<std::string>(obj["type"]));
    }
    if (obj.has_key("fast_axis") &&
        obj.has_key("slow_axis") &&
        obj.has_key("origin")) {
      result->set_local_frame(
        boost::python::extract< vec3<double> >(obj["fast_axis"]),
        boost::python::extract< vec3<double> >(obj["slow_axis"]),
        boost::python::extract< vec3<double> >(obj["origin"]));
    }
    return result;
  }

  template <>
  boost::python::dict to_dict<Panel>(const Panel &obj) {
    boost::python::dict result;
    result["name"] = obj.get_name();
    result["type"] = obj.get_type();
    result["fast_axis"] = obj.get_local_fast_axis();
    result["slow_axis"] = obj.get_local_slow_axis();
    result["origin"] = obj.get_local_origin();
    result["image_size"] = obj.get_image_size();
    result["pixel_size"] = obj.get_pixel_size();
    result["trusted_range"] = obj.get_trusted_range();
    result["mask"] = boost::python::list(obj.get_mask());
    return result;
  }

  template <>
  Panel* from_dict<Panel>(boost::python::dict obj) {
    Panel *result = new Panel();
    if (obj.has_key("name")) {
      result->set_name(boost::python::extract<std::string>(obj["name"]));
    }
    if (obj.has_key("type")) {
      result->set_type(boost::python::extract<std::string>(obj["type"]));
    }
    if (obj.has_key("fast_axis") &&
        obj.has_key("slow_axis") &&
        obj.has_key("origin")) {
      result->set_local_frame(
        boost::python::extract< vec3<double> >(obj["fast_axis"]),
        boost::python::extract< vec3<double> >(obj["slow_axis"]),
        boost::python::extract< vec3<double> >(obj["origin"]));
    }
    if (obj.has_key("mask")) {
      scitbx::af::shared<int4> mask =
        boost::python::extract< scitbx::af::shared<int4> >(
          boost::python::extract<boost::python::list>(obj["mask"]));
      result->set_mask(mask.const_ref());
    }
    result->set_image_size(
      boost::python::extract< vec2<std::size_t> >(obj["image_size"]));
    result->set_pixel_size(
      boost::python::extract< vec2<double> >(obj["pixel_size"]));
    result->set_trusted_range(
      boost::python::extract< vec2<double> >(obj["trusted_range"]));

    return result;
  }

  void export_panel()
  {
    using namespace boost::python;

    class_<PanelFrame>("PanelFrame")
      .def("set_frame",
        &PanelFrame::set_frame, (
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin")))
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
      .def("get_local_fast_axis",
        &PanelFrame::get_local_fast_axis)
      .def("get_local_slow_axis",
        &PanelFrame::get_local_slow_axis)
      .def("get_local_origin",
        &PanelFrame::get_local_origin)
      .def("get_parent_fast_axis",
        &PanelFrame::get_parent_fast_axis)
      .def("get_parent_slow_axis",
        &PanelFrame::get_parent_slow_axis)
      .def("get_parent_origin",
        &PanelFrame::get_parent_origin)
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
      .def("__ne__", &PanelBase::operator!=)
      .def("to_dict", &to_dict<PanelBase>)
      .def("from_dict", &from_dict<PanelBase>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(PanelBasePickleSuite());

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
      .def("__ne__", &Panel::operator!=)
      .def("__deepcopy__", &panel_deepcopy)
      .def("__copy__", &panel_deepcopy)
      .def("__str__", &panel_to_string)
      .def("to_dict", &to_dict<Panel>)
      .def("from_dict", &from_dict<Panel>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def("is_", &panel_is)
      .def("is_similar_to", &Panel::is_similar_to)
      .def_pickle(PanelPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
