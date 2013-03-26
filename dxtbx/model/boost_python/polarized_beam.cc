/*
 * polarized_beam.cc
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
#include <dxtbx/model/polarized_beam.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string polarized_beam_to_string(const PolarizedBeam &beam) {
    std::stringstream ss;
    ss << beam;
    return ss.str();
  }

  struct PolarizedBeamPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const PolarizedBeam &obj) {
      return boost::python::make_tuple(
        obj.get_direction(),
        obj.get_wavelength(),
        obj.get_polarization(),
        obj.get_polarization_fraction());
    }
  };

  void export_polarized_beam() 
  {
    // Export PolarizedBeam : Beam
    class_ <PolarizedBeam, bases <Beam> > ("PolarizedBeam")
      .def(init <vec3 <double>, 
                 double,
                 vec3 <double>, 
                 double> ((
          arg("direction"),
          arg("wavelength"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .def(init <vec3 <double>, 
                 vec3 <double>, 
                 double> ((
          arg("s0"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .def("get_polarization",
        &PolarizedBeam::get_polarization)
      .def("set_polarization",
        &PolarizedBeam::set_polarization)
      .def("get_polarization_fraction",
        &PolarizedBeam::get_polarization_fraction)
      .def("set_polarization_fraction",
        &PolarizedBeam::set_polarization_fraction)
      .def("__str__", &polarized_beam_to_string)
      .def_pickle(PolarizedBeamPickleSuite());
  }

}}} // namespace = dxtbx::model::boost_python
