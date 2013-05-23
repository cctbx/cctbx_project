/*
 * attenuation_coefficient.cc
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
#include <dxtbx/model/attenuation_coefficient.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_attenuation_coefficient()
  {
    class_<XrayMassCoeffTable>("XrayMassCoeffTable", no_init)
      .def(init<std::size_t>((arg("z"))))
      .def("size", &XrayMassCoeffTable::size)
      .def("energy", &XrayMassCoeffTable::energy)
      .def("mu_rho", &XrayMassCoeffTable::mu_rho)
      .def("mu_en_rho", &XrayMassCoeffTable::mu_en_rho)
      .def("mu_rho_at_ev", &XrayMassCoeffTable::mu_rho_at_ev)
      .def("mu_rho_at_kev", &XrayMassCoeffTable::mu_rho_at_kev)
      .def("mu_rho_at_angstrom", &XrayMassCoeffTable::mu_rho_at_angstrom)
      .def("mu_en_rho_at_ev", &XrayMassCoeffTable::mu_rho_at_ev)
      .def("mu_en_rho_at_kev", &XrayMassCoeffTable::mu_rho_at_kev)
      .def("mu_en_rho_at_angstrom", &XrayMassCoeffTable::mu_rho_at_angstrom);
  }

}}} // namespace = dxtbx::model::boost_python
