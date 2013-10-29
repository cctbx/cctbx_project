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
#include <dxtbx/model/detector2.h>

namespace dxtbx { namespace model { namespace boost_python {

  void detector2_set_item(Detector2 &d, std::size_t i, const Panel2 &v) {
    d[i] = v;
  }

  Panel2& detector2_get_item(Detector2 &d, std::size_t i) {
    return d[i];
  }

  void export_detector2() 
  {
    using namespace boost::python;
      
    // Export a Detector base class
    class_ <Detector2> ("DetectorBase2")
      .def("add_panel",
        &Detector2::add_panel,
        return_internal_reference<>())
      .def("__len__",
        &Detector2::size)
      .def("__setitem__",
        &detector2_set_item)
      .def("__getitem__", 
        &detector2_get_item, 
        return_internal_reference<>())
      .def("__iter__", 
        iterator<Detector2, return_internal_reference<> >())
      .def("__eq__", &Detector2::operator==)
      .def("__ne__", &Detector2::operator!=);
  }

}}} // namespace dxtbx::model::boost_python
