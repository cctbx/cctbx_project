/*
 * goniometer.cc
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
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <sstream>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/kappa_goniometer.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string goniometer_to_string(const Goniometer &goniometer) {
    std::stringstream ss;
    ss << goniometer;
    return ss.str();
  }

  void export_goniometer() 
  {
    class_ <GoniometerBase> ("GoniometerBase");

    class_ <Goniometer, bases <GoniometerBase> > ("Goniometer")
      .def(init <vec3 <double> > ((
          arg("rotation_axis"))))
      .def(init <vec3 <double>,
                 mat3 <double> > ((
          arg("rotation_axis"), 
          arg("fixed_rotation_matrix"))))
      .def("get_rotation_axis",  
        &Goniometer::get_rotation_axis)
      .def("set_rotation_axis",
        &Goniometer::set_rotation_axis)
      .def("get_fixed_rotation",  
        &Goniometer::get_fixed_rotation)
      .def("set_fixed_rotation",
        &Goniometer::set_fixed_rotation)
      .def("__eq__", &Goniometer::operator==)
      .def("__ne__", &Goniometer::operator!=)
      .def("__str__", &goniometer_to_string);
  }

}}} // namespace dxtbx::model::boost_python
