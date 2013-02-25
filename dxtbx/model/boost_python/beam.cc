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
#include <boost/format.hpp>
#include <string>
#include <dxtbx/model/beam.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string beam_to_string(const Beam &beam) {
    boost::format fmt(
      "Beam:\n"
      "    direction :            (%1%, %2%, %3%)\n"
      "    wavelength:            %4%");
        
    fmt % beam.get_direction()[0];
    fmt % beam.get_direction()[1];
    fmt % beam.get_direction()[2];
    fmt % beam.get_wavelength();
    return fmt.str();
  }

  std::string polarized_beam_to_string(const PolarizedBeam &beam) {
    boost::format fmt(
      "%1%\n"
      "    polarization:          (%2%, %3%, %4%)\n"
      "    polarization fraction: %5%");
        
    fmt % beam_to_string(beam);
    fmt % beam.get_polarization()[0];
    fmt % beam.get_polarization()[1];
    fmt % beam.get_polarization()[2];
    fmt % beam.get_polarization_fraction();
    return fmt.str();
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
          arg("direction"),
          arg("wavelength"))))
      .add_property("direction",  
        &Beam::get_direction,
        &Beam::set_direction)
      .add_property("wavelength", 
        &Beam::get_wavelength)
      .def("__eq__", &Beam::operator==)
      .def("__ne__", &Beam::operator!=)
      .def("__str__", &beam_to_string);

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
          arg("direction"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .add_property("polarization",
        &PolarizedBeam::get_polarization,
        &PolarizedBeam::set_polarization)
      .add_property("polarization_fraction",
        &PolarizedBeam::get_polarization_fraction,
        &PolarizedBeam::set_polarization_fraction)
      .def("__str__", &polarized_beam_to_string);
  }

}}} // namespace = dxtbx::model::boost_python
