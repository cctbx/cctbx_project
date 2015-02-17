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
#include <scitbx/constants.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using scitbx::deg_as_rad;

  std::string detector_to_string(const Detector &detector) {
    std::stringstream ss;
    ss << detector;
    return ss.str();
  }

  static
  void detector_set_item(Detector &d, std::size_t i, const Panel &v) {
    d[i] = v;
  }

  static
  Panel* detector_get_item(Detector &d, std::size_t i) {
    return d.at(i);
  }

  static
  scitbx::af::shared<std::string> get_names(const Detector &d) {
    scitbx::af::shared<std::string> result(d.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = d[i].get_name();
    }
    return result;
  }

  static
  void rotate_around_origin(Detector &detector, vec3<double> axis, double angle, bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    detector.rotate_around_origin(axis, angle_rad);
  }

  template <>
  boost::python::dict to_dict<Detector>(const Detector &obj) {
    boost::python::dict result;
    boost::python::list panels;
    for (std::size_t i = 0; i < obj.size(); ++i) {
      panels.append(to_dict(obj[i]));
    }
    result["panels"] = panels;
    return result;
  }

  template <>
  Detector* from_dict<Detector>(boost::python::dict obj) {
    Detector *result = new Detector();
    boost::python::list panels =
      boost::python::extract<boost::python::list>(obj["panels"]);
    for (std::size_t i = 0; i < boost::python::len(panels); ++i) {
      result->add_panel(from_dict<Panel>(
        boost::python::extract<boost::python::dict>(panels[i])));
    }
    return result;
  }

  void export_detector()
  {
    using namespace boost::python;

    // Export a Detector base class
    class_ <Detector> ("Detector")
      .def(init<const Panel&>())
      .def("add_panel",
        (Panel&(Detector::*)())&Detector::add_panel,
        return_internal_reference<>())
      .def("add_panel",
        (Panel&(Detector::*)(const Panel&))&Detector::add_panel,
        return_internal_reference<>())
      .def("add_panel_by_pointer",
        (void(Detector::*)(Panel*))&Detector::add_panel,
        with_custodian_and_ward<1,2>())
      .def("__len__",
        &Detector::size)
      .def("__setitem__",
        &detector_set_item)
      .def("__getitem__",
        &detector_get_item,
        return_internal_reference<>())
      .def("__iter__",
        iterator<Detector, return_internal_reference<> >())
      .def("__eq__", &Detector::operator==)
      .def("__ne__", &Detector::operator!=)
      .def("is_similar_to", &Detector::is_similar_to)
      .def("get_max_resolution",
        &Detector::get_max_resolution, (arg("s0")))
      .def("get_ray_intersection",
        &Detector::get_ray_intersection, (arg("s1")))
      .def("get_panel_intersection",
        &Detector::get_panel_intersection, (arg("s1")))
      //.def("do_panels_intersect",
      //  &Detector::do_panels_intersect)
      .def("get_names", &get_names)
      .def("rotate_around_origin",
          &rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("__str__", &detector_to_string)
      .def("to_dict", &to_dict<Detector>)
      .def("from_dict", &from_dict<Detector>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .enable_pickling();

    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();
  }

}}} // namespace dxtbx::model::boost_python
