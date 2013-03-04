/*
 * beam.cc
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
#include <sstream>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/polarized_beam.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string beam_to_string(const Beam &beam) {
    std::stringstream ss;
    ss << beam;
    return ss.str();
  }

  void export_beam()
  {
    // Export BeamBase
    class_ <BeamBase> ("BeamBase");

    // Export Beam : BeamBase
    class_ <Beam, bases <BeamBase> > ("Beam")
      .def(init <vec3 <double>,
                 double> ((
          arg("direction"), 
          arg("wavelength"))))
      .def(init <vec3 <double> > ((
          arg("s0"))))
      .def("get_direction", 
        &Beam::get_direction)
      .def("set_direction",
        &Beam::set_direction)
      .def("get_wavelength", 
        &Beam::get_wavelength)
      .def("set_wavelength",
        &Beam::set_wavelength)
      .def("get_s0",
        &Beam::get_s0)
      .def("__eq__", &Beam::operator==)
      .def("__ne__", &Beam::operator!=)
      .def("__str__", &beam_to_string);
  }

}}} // namespace = dxtbx::model::boost_python
