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
#include <cctbx/eltbx/attenuation_coefficient.h>

namespace cctbx { namespace eltbx { namespace attenuation_coefficient {
    namespace boost_python {

  using namespace boost::python;

  void init_module()
  {
    double (table::*energy_single)(std::size_t) const = &table::energy;
    double (table::*mu_rho_single)(std::size_t) const = &table::mu_rho;
    double (table::*mu_en_rho_single)(std::size_t) const = &table::mu_en_rho;
    flex_double (table::*energy_list)() const = &table::energy;
    flex_double (table::*mu_rho_list)() const = &table::mu_rho;
    flex_double (table::*mu_en_rho_list)() const = &table::mu_en_rho;

    class_<table>("table", no_init)
      .def(init<std::size_t>((arg("z"))))
      .def("size", &table::size)
      .def("density", &table::density)
      .def("min_energy", &table::min_energy)
      .def("max_energy", &table::max_energy)
      .def("energy", energy_single)
      .def("energy", energy_list)
      .def("mu_rho", mu_rho_single)
      .def("mu_rho", mu_rho_list)
      .def("mu_en_rho", mu_en_rho_single)
      .def("mu_en_rho", mu_en_rho_list)
      .def("mu_rho_at_ev", &table::mu_rho_at_ev)
      .def("mu_rho_at_kev", &table::mu_rho_at_kev)
      .def("mu_rho_at_angstrom", &table::mu_rho_at_angstrom)
      .def("mu_en_rho_at_ev", &table::mu_en_rho_at_ev)
      .def("mu_en_rho_at_kev", &table::mu_en_rho_at_kev)
      .def("mu_en_rho_at_angstrom", &table::mu_en_rho_at_angstrom)
      .def("mu_at_ev", &table::mu_at_ev)
      .def("mu_at_kev", &table::mu_at_kev)
      .def("mu_at_angstrom", &table::mu_at_angstrom)
      .def("lambda_at_ev", &table::lambda_at_ev)
      .def("lambda_at_kev", &table::lambda_at_kev)
      .def("lambda_at_angstrom", &table::lambda_at_angstrom);
  }

}}}} // namespace = cctbx::eltbx::attenuation_coefficient::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_attenuation_coefficient_ext)
{
  cctbx::eltbx::attenuation_coefficient::boost_python::init_module();
}
