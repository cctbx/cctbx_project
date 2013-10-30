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

  static
  void detector_set_item(Detector &d, std::size_t i, const Panel &v) {
    d[i] = v;
  }

  static
  Panel& detector_get_item(Detector &d, std::size_t i) {
    return d[i];
  }
  
  static
  scitbx::af::shared<std::string> get_names(const Detector &d) {
    scitbx::af::shared<std::string> result(d.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = d[i].get_name();
    }
    return result;
  }

  void export_detector() 
  {
    using namespace boost::python;
      
    // Export a Detector base class
    class_ <Detector> ("DetectorBase")
      .def(init<const Panel&>())
      .def(init<const scitbx::af::const_ref<Panel>&>())
      .def("add_panel",
        (Panel&(Detector::*)())&Detector::add_panel,
        return_internal_reference<>())
      .def("add_panel",
        (Panel&(Detector::*)(const Panel&))&Detector::add_panel,
        return_internal_reference<>())
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
      .def("get_max_resolution",
        &Detector::get_max_resolution, (arg("s0")))
      .def("get_ray_intersection",
        &Detector::get_ray_intersection, (arg("s1")))
      //.def("do_panels_intersect",
      //  &Detector::do_panels_intersect)
      .def("get_names", &get_names);
      
    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();
  }

}}} // namespace dxtbx::model::boost_python
