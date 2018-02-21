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
#include <scitbx/constants.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/boost_python/to_from_dict.h>
#include <dxtbx/model/boost_python/pickle_suite.h>

namespace dxtbx { namespace model { namespace boost_python {

  using scitbx::deg_as_rad;

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
  scitbx::af::shared<vec2<double> >
  pixel_to_millimeter_multiple(const Panel &panel,
                scitbx::af::flex<vec2<double> >::type const& xy) {
    scitbx::af::shared<vec2<double> > result((scitbx::af::reserve(xy.size())));
    for(std::size_t i = 0; i < xy.size(); i++) {
      result.push_back(panel.pixel_to_millimeter(xy[i]));
    }
    return result;
  }

  static
  scitbx::af::shared<vec2<double> >
  millimeter_to_pixel_multiple(const Panel &panel,
                scitbx::af::flex<vec2<double> >::type const& xy) {
    scitbx::af::shared<vec2<double> > result((scitbx::af::reserve(xy.size())));
    for(std::size_t i = 0; i < xy.size(); i++) {
      result.push_back(panel.millimeter_to_pixel(xy[i]));
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

  static
  void rotate_around_origin(Panel &panel, vec3<double> axis, double angle, bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    panel.rotate_around_origin(axis, angle_rad);
  }

  template <>
  boost::python::dict to_dict<VirtualPanel>(const VirtualPanel &obj) {
    boost::python::dict result;
    result["name"] = obj.get_name();
    result["type"] = obj.get_type();
    result["fast_axis"] = obj.get_local_fast_axis();
    result["slow_axis"] = obj.get_local_slow_axis();
    result["origin"] = obj.get_local_origin();
    return result;
  }

  template <>
  VirtualPanel* from_dict<VirtualPanel>(boost::python::dict obj) {
    VirtualPanel *result = new VirtualPanel();
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
  boost::python::dict to_dict< shared_ptr<PxMmStrategy> >(
      const shared_ptr<PxMmStrategy> &obj) {
    boost::python::dict result;
    std::string name = obj->name();
    result["type"] = name;
    /* if (name == "SimplePxMmStrategy") { */
    /* } else if (name == "ParallaxCorrectedPxMmStrategy") { */
    /*   boost::shared_ptr<ParallaxCorrectedPxMmStrategy> d = */
    /*     boost::dynamic_pointer_cast<ParallaxCorrectedPxMmStrategy>(obj); */
    /*   result["mu"] = d->mu(); */
    /*   result["t0"] = d->t0(); */
    /* } else if (name == "OffsetParallaxCorrectedPxMmStrategy") { */
    /*   boost::shared_ptr<OffsetParallaxCorrectedPxMmStrategy> d = */
    /*     boost::dynamic_pointer_cast<OffsetParallaxCorrectedPxMmStrategy>(obj); */
    /*   result["mu"] = d->mu(); */
    /*   result["t0"] = d->t0(); */
    /*   result["dx"] = boost::python::list(d->dx()); */
    /*   result["dy"] = boost::python::list(d->dy()); */
    /* } else { */
    /*   throw DXTBX_ERROR("Unknown PxMmStrategy"); */
    /* } */
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
    result["raw_image_offset"] = obj.get_raw_image_offset();
    result["image_size"] = obj.get_image_size();
    result["pixel_size"] = obj.get_pixel_size();
    result["trusted_range"] = obj.get_trusted_range();
    result["thickness"] = obj.get_thickness();
    result["material"] = obj.get_material();
    result["mu"] = obj.get_mu();
    result["identifier"] = obj.get_identifier();
    result["mask"] = boost::python::list(obj.get_mask());
    result["gain"] = obj.get_gain();
    result["px_mm_strategy"] = to_dict(obj.get_px_mm_strategy());
    return result;
  }

  Panel* basic_panel_from_dict(boost::python::dict obj) {
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
    if (obj.has_key("thickness")) {
      result->set_thickness(boost::python::extract<double>(obj["thickness"]));
    }
    if (obj.has_key("material")) {
      result->set_material(boost::python::extract<std::string>(obj["material"]));
    }
    if (obj.has_key("mu")) {
      result->set_mu(boost::python::extract<double>(obj["mu"]));
    }
    if (obj.has_key("identifier")) {
      result->set_identifier(boost::python::extract<std::string>(obj["identifier"]));
    }
    if (obj.has_key("mask")) {
      scitbx::af::shared<int4> mask =
        boost::python::extract< scitbx::af::shared<int4> >(
          boost::python::extract<boost::python::list>(obj["mask"]));
      result->set_mask(mask.const_ref());
    }
    if (obj.has_key("gain")) {
      result->set_gain(boost::python::extract<double>(obj["gain"]));
    }
    if (obj.has_key("raw_image_offset")) {
      result->set_raw_image_offset(
        boost::python::extract<int2>(obj["raw_image_offset"]));
    }
    if (obj.has_key("image_size")) {
      result->set_image_size(
        boost::python::extract< tiny<std::size_t,2> >(obj["image_size"]));
    }
    if (obj.has_key("pixel_size")) {
      result->set_pixel_size(
        boost::python::extract< tiny<double,2> >(obj["pixel_size"]));
    }
    if (obj.has_key("trusted_range")) {
      result->set_trusted_range(
        boost::python::extract< tiny<double,2> >(obj["trusted_range"]));
    }
    return result;
  }

  template <>
  Panel* from_dict<Panel>(boost::python::dict obj) {
    Panel *result = basic_panel_from_dict(obj);
    if (obj.has_key("px_mm_strategy")) {
      boost::python::dict st = boost::python::extract<boost::python::dict>(obj["px_mm_strategy"]);
      std::string name = boost::python::extract<std::string>(st["type"]);
      if (name == "SimplePxMmStrategy") {
        shared_ptr<PxMmStrategy> strategy(new SimplePxMmStrategy());
        result->set_px_mm_strategy(strategy);
      } else if (name == "ParallaxCorrectedPxMmStrategy") {
        double mu = result->get_mu();
        double t0 = result->get_thickness();
        if (st.has_key("mu") && st.has_key("t0")) {
          mu = boost::python::extract<double>(st["mu"]);
          t0 = boost::python::extract<double>(st["t0"]);
          result->set_mu(mu);
          result->set_thickness(t0);
        }
        if (mu > 0 && t0 > 0) {
          shared_ptr<PxMmStrategy> strategy(
              new ParallaxCorrectedPxMmStrategy(mu, t0));
          result->set_px_mm_strategy(strategy);
        }
      } else if (name == "OffsetPxMmStrategy") {
        shared_ptr<PxMmStrategy> strategy(
            new SimplePxMmStrategy());
        result->set_px_mm_strategy(strategy);
        //throw DXTBX_ERROR("The offset arrays need to be set for OffsetParallaxCorrectedPxMmStrategy");
      } else if (name == "OffsetParallaxCorrectedPxMmStrategy") {
        double mu = result->get_mu();
        double t0 = result->get_thickness();
        if (st.has_key("mu") && st.has_key("t0")) {
          mu = boost::python::extract<double>(st["mu"]);
          t0 = boost::python::extract<double>(st["t0"]);
          result->set_mu(mu);
          result->set_thickness(t0);
        }
        if (mu > 0 && t0 > 0) {
          shared_ptr<PxMmStrategy> strategy(
              new ParallaxCorrectedPxMmStrategy(mu, t0));
          result->set_px_mm_strategy(strategy);
          //throw DXTBX_ERROR("The offset arrays need to be set for OffsetParallaxCorrectedPxMmStrategy");
        }
      } else {
        throw DXTBX_ERROR("Unknown px_mm_strategy");
      }
    }
    return result;
  }

  Panel* panel_from_dict_with_offset(
        boost::python::dict obj,
        scitbx::af::versa<double, scitbx::af::c_grid<2> > dx,
        scitbx::af::versa<double, scitbx::af::c_grid<2> > dy) {
    Panel *result = basic_panel_from_dict(obj);
    DXTBX_ASSERT(dx.accessor()[0] == result->get_image_size()[1]);
    DXTBX_ASSERT(dx.accessor()[1] == result->get_image_size()[0]);
    DXTBX_ASSERT(dx.accessor().all_eq(dy.accessor()));
    if (obj.has_key("px_mm_strategy")) {
      boost::python::dict st = boost::python::extract<boost::python::dict>(obj["px_mm_strategy"]);
      std::string name = boost::python::extract<std::string>(st["type"]);
      if (name == "SimplePxMmStrategy") {
        shared_ptr<PxMmStrategy> strategy(new SimplePxMmStrategy());
        result->set_px_mm_strategy(strategy);
      } else if (name == "ParallaxCorrectedPxMmStrategy") {
        double mu = result->get_mu();
        double t0 = result->get_thickness();
        if (st.has_key("mu") && st.has_key("t0")) {
          mu = boost::python::extract<double>(st["mu"]);
          t0 = boost::python::extract<double>(st["t0"]);
          result->set_mu(mu);
          result->set_thickness(t0);
        }
        if (mu > 0 && t0 > 0) {
          shared_ptr<PxMmStrategy> strategy(
              new ParallaxCorrectedPxMmStrategy(mu, t0));
         result->set_px_mm_strategy(strategy);
        }
      } else if (name == "OffsetPxMmStrategy") {
        shared_ptr<PxMmStrategy> strategy(
            new OffsetPxMmStrategy(
              dx,
              dy));
        result->set_px_mm_strategy(strategy);
      } else if (name == "OffsetParallaxCorrectedPxMmStrategy") {
        double mu = result->get_mu();
        double t0 = result->get_thickness();
        if (st.has_key("mu") && st.has_key("t0")) {
          mu = boost::python::extract<double>(st["mu"]);
          t0 = boost::python::extract<double>(st["t0"]);
          result->set_mu(mu);
          result->set_thickness(t0);
        }
        if (mu > 0 && t0 > 0) {
          shared_ptr<PxMmStrategy> strategy(
              new OffsetParallaxCorrectedPxMmStrategy(
                mu,
                t0,
                dx,
                dy));
          result->set_px_mm_strategy(strategy);
        }
      } else {
        DXTBX_ASSERT(false);
      }
    }
    return result;
  }

  Panel* panel_from_dict_with_offset_wrapper(
        boost::python::dict obj,
        scitbx::af::versa<double, scitbx::af::flex_grid<> > dx,
        scitbx::af::versa<double, scitbx::af::flex_grid<> > dy) {
    DXTBX_ASSERT(dx.accessor().all().size() == 2);
    DXTBX_ASSERT(dy.accessor().all().size() == 2);
    DXTBX_ASSERT(dx.accessor().all().all_eq(dy.accessor().all()));

    std::size_t ysize = dx.accessor().all()[0];
    std::size_t xsize = dx.accessor().all()[1];

    scitbx::af::c_grid<2> grid(ysize, xsize);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dx2(dx.handle(), grid);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dy2(dy.handle(), grid);

    return panel_from_dict_with_offset(obj, dx2, dy2);
  }

  void export_panel()
  {
    using namespace boost::python;

    class_<VirtualPanelFrame>("VirtualPanelFrame")
      .def("set_frame",
        &VirtualPanelFrame::set_frame, (
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin")))
      .def("set_local_frame",
        &VirtualPanelFrame::set_local_frame, (
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin")))
      .def("set_parent_frame",
        &VirtualPanelFrame::set_parent_frame, (
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin")))
      .def("get_local_fast_axis",
        &VirtualPanelFrame::get_local_fast_axis)
      .def("get_local_slow_axis",
        &VirtualPanelFrame::get_local_slow_axis)
      .def("get_local_origin",
        &VirtualPanelFrame::get_local_origin)
      .def("get_parent_fast_axis",
        &VirtualPanelFrame::get_parent_fast_axis)
      .def("get_parent_slow_axis",
        &VirtualPanelFrame::get_parent_slow_axis)
      .def("get_parent_origin",
        &VirtualPanelFrame::get_parent_origin)
      .def("get_local_d_matrix",
        &VirtualPanelFrame::get_local_d_matrix)
      .def("get_parent_d_matrix",
        &VirtualPanelFrame::get_parent_d_matrix)
      .def("get_d_matrix",
        &VirtualPanelFrame::get_d_matrix)
      .def("get_D_matrix",
        &VirtualPanelFrame::get_D_matrix)
      .def("get_origin",
        &VirtualPanelFrame::get_origin)
      .def("get_fast_axis",
        &VirtualPanelFrame::get_fast_axis)
      .def("get_slow_axis",
        &VirtualPanelFrame::get_slow_axis)
      .def("get_normal",
        &VirtualPanelFrame::get_normal)
      .def("get_normal_origin",
        &VirtualPanelFrame::get_normal_origin)
      .def("get_distance",
        &VirtualPanelFrame::get_distance)
      .def("get_directed_distance",
        &VirtualPanelFrame::get_directed_distance)
      .def("get_beam_centre",
        &VirtualPanelFrame::get_beam_centre, (
          arg("s0")))
      .def("get_beam_centre_lab",
        &VirtualPanelFrame::get_beam_centre_lab, (
          arg("s0")))
      .def("get_lab_coord",
        &VirtualPanelFrame::get_lab_coord, (
          arg("xy")))
      .def("get_lab_coord",
        &get_lab_coord_multiple, (
          arg("xy")))
      .def("get_ray_intersection",
        &VirtualPanelFrame::get_ray_intersection, (
          arg("s1")))
      .def("get_bidirectional_ray_intersection",
        &VirtualPanelFrame::get_bidirectional_ray_intersection, (
          arg("s1")))
      .def("__eq__", &VirtualPanelFrame::operator==)
      .def("__ne__", &VirtualPanelFrame::operator!=);

    class_<VirtualPanel, bases<VirtualPanelFrame> >("VirtualPanel")
      .def("get_name", &VirtualPanel::get_name)
      .def("set_name", &VirtualPanel::set_name)
      .def("get_type", &VirtualPanel::get_type)
      .def("set_type", &VirtualPanel::set_type)
      .def("__eq__", &VirtualPanel::operator==)
      .def("__ne__", &VirtualPanel::operator!=)
      /* .def("to_dict", &to_dict<VirtualPanel>) */
      /* .def("from_dict", &from_dict<VirtualPanel>, */
      /*   return_value_policy<manage_new_object>()) */
      /* .staticmethod("from_dict") */
      /* .def_pickle(VirtualPanelPickleSuite()) */
      ;

    class_<PanelData, bases<VirtualPanel> >("Panel")
      .def(init<std::string,
                std::string,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,2>,
                tiny<std::size_t,2>,
                tiny<double,2>,
                double,
                std::string,
                double>((
        arg("type"),
        arg("name"),
        arg("fast_axis"),
        arg("slow_axis"),
        arg("origin"),
        arg("pixel_size"),
        arg("image_size"),
        arg("trusted_range"),
        arg("thickness"),
        arg("material"),
        arg("mu")=0.0)))
      .def("get_pixel_size", &PanelData::get_pixel_size)
      .def("set_pixel_size", &PanelData::set_pixel_size)
      .def("get_image_size", &PanelData::get_image_size)
      .def("set_image_size", &PanelData::set_image_size)
      .def("get_trusted_range", &PanelData::get_trusted_range)
      .def("set_trusted_range", &PanelData::set_trusted_range)
      .def("get_thickness", &PanelData::get_thickness)
      .def("set_thickness", &PanelData::set_thickness)
      .def("get_material", &PanelData::get_material)
      .def("set_material", &PanelData::set_material)
      .def("get_mu", &PanelData::get_mu)
      .def("set_mu", &PanelData::set_mu)
      .def("get_raw_image_offset", &PanelData::get_raw_image_offset)
      .def("set_raw_image_offset", &PanelData::set_raw_image_offset)
      .def("set_mask", &PanelData::set_mask)
      .def("get_mask", &PanelData::get_mask)
      .def("add_mask", &PanelData::add_mask)
      .def("__eq__", &PanelData::operator==)
      .def("__ne__", &PanelData::operator!=)
      .def("is_similar_to", &PanelData::is_similar_to, (
            arg("other"),
            arg("fast_axis_tolerance")=1e-6,
            arg("slow_axis_tolerance")=1e-6,
            arg("origin_tolerance")=1e-6,
            arg("static_only")=false,
            arg("ignore_trusted_range")=false))
      ;

    class_<Panel, bases<PanelData> >("Panel")
      .def(init<std::string,
                std::string,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,2>,
                tiny<std::size_t,2>,
                tiny<double,2>,
                double,
                std::string,
                double,
                std::string>((
        arg("type"),
        arg("name"),
        arg("fast_axis"),
        arg("slow_axis"),
        arg("origin"),
        arg("pixel_size"),
        arg("image_size"),
        arg("trusted_range"),
        arg("thickness"),
        arg("material"),
        arg("mu")=0.0,
        arg("identifier")="")))
      .def(init<std::string,
                std::string,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,3>,
                tiny<double,2>,
                tiny<std::size_t,2>,
                tiny<double,2>,
                double,
                std::string,
                shared_ptr<PxMmStrategy>,
                double,
                std::string>((
        arg("type"),
        arg("name"),
        arg("fast_axis"),
        arg("slow_axis"),
        arg("origin"),
        arg("pixel_size"),
        arg("image_size"),
        arg("trusted_range"),
        arg("thickness"),
        arg("material"),
        arg("px_mm"),
        arg("mu")=0.0,
        arg("identifier")="")))
      .def("get_identifier", &Panel::get_identifier)
      .def("set_identifier", &Panel::set_identifier)
      .def("get_gain", &Panel::get_gain)
      .def("set_gain", &Panel::set_gain)
      .def("get_image_size_mm", &Panel::get_image_size_mm)
      .def("get_px_mm_strategy", &Panel::get_px_mm_strategy)
      .def("set_px_mm_strategy", &Panel::set_px_mm_strategy)
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
      .def("millimeter_to_pixel", millimeter_to_pixel_multiple)
      .def("pixel_to_millimeter", &Panel::pixel_to_millimeter)
      .def("pixel_to_millimeter", pixel_to_millimeter_multiple)
      .def("get_two_theta_at_pixel", &Panel::get_two_theta_at_pixel)
      .def("get_two_theta_array", &Panel::get_two_theta_array)
      .def("get_resolution_at_pixel", &Panel::get_resolution_at_pixel)
      .def("get_max_resolution_at_corners",
        &Panel::get_max_resolution_at_corners)
      .def("get_max_resolution_ellipse", &Panel::get_max_resolution_ellipse)
      .def("__deepcopy__", &panel_deepcopy)
      .def("__copy__", &panel_deepcopy)
      .def("__str__", &panel_to_string)
      .def("rotate_around_origin",
          &rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("get_trusted_range_mask", &Panel::get_trusted_range_mask<int>)
      .def("get_trusted_range_mask", &Panel::get_trusted_range_mask<double>)
      .def("get_untrusted_rectangle_mask", &Panel::get_untrusted_rectangle_mask)
      .def("to_dict", &to_dict<Panel>)
      .def("from_dict", &from_dict<Panel>,
        return_value_policy<manage_new_object>())
      .def("from_dict", &panel_from_dict_with_offset_wrapper,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def("is_", &panel_is)
      .def_pickle(PanelPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
