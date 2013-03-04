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

  void export_detector() 
  {
    // Register std::pair conversion for Detector coordinate type 
    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();

    // Export a Detector class
    class_ <Detector> ("Detector")
      .def(init<const Panel&>((
        arg("panel"))))
      .def(init<const Detector::panel_list_type&>((
        arg("panel_list"))))
      .def("add_panel",
        &Detector::add_panel, (
          arg("panel")))
      .def("num_panels",
        &Detector::num_panels)
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
  }

}}} // namespace dials::model::boost_python
