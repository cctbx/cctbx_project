/*
 * kappa_goniometer.cc
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

  std::string kappa_goniometer_to_string(const KappaGoniometer &goniometer) {
    std::stringstream ss;
    ss << goniometer;
    return ss.str();
  }

  struct KappaGoniometerPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const KappaGoniometer &obj) {
      KappaGoniometer::Direction direction(obj.get_direction());
      std::string direction_str;
      if (direction == KappaGoniometer::PlusY) {
        direction_str = "+y";
      } else if (direction == KappaGoniometer::PlusZ) {
        direction_str = "+z";
      } else if (direction == KappaGoniometer::MinusY) {
        direction_str = "-y";
      } else if (direction == KappaGoniometer::MinusZ) {
        direction_str = "-y";
      } else {
        direction_str = "";
      }

      KappaGoniometer::ScanAxis scan_axis(obj.get_scan_axis());
      std::string scan_axis_str;

      if (scan_axis == KappaGoniometer::Omega) {
        scan_axis_str = "omega";
      } else if (scan_axis == KappaGoniometer::Phi) {
        scan_axis_str = "phi";
      } else {
        scan_axis_str = "";
      }

      return boost::python::make_tuple(
        obj.get_alpha_angle(),
        obj.get_omega_angle(),
        obj.get_kappa_angle(),
        obj.get_phi_angle(),
        direction_str,
        scan_axis_str);
    }
  };

  static boost::shared_ptr<KappaGoniometer> make_kappa_goniometer(
      double alpha, double omega, double kappa, double phi,
      std::string direction_str, std::string scan_axis_str)
  {
    KappaGoniometer::Direction direction = KappaGoniometer::NoDirection;
    boost::algorithm::to_lower(direction_str);
    if (direction_str == "+y") {
      direction = KappaGoniometer::PlusY;
    } else if (direction_str == "+z") {
      direction = KappaGoniometer::PlusZ;
    } else if (direction_str == "-y")  {
      direction = KappaGoniometer::MinusY;
    } else if (direction_str == "-z") {
      direction = KappaGoniometer::MinusZ;
    } else {
      direction = KappaGoniometer::NoDirection;
    }

    KappaGoniometer::ScanAxis scan_axis = KappaGoniometer::NoAxis;
    boost::algorithm::to_lower(scan_axis_str);
    if (scan_axis_str == "omega") {
      scan_axis = KappaGoniometer::Omega;
    } else if (scan_axis_str == "phi") {
      scan_axis = KappaGoniometer::Phi;
    } else {
      scan_axis = KappaGoniometer::NoAxis;
    }

    return boost::shared_ptr<KappaGoniometer>(new KappaGoniometer(
      alpha, omega, kappa, phi, direction, scan_axis));
  }

  void export_kappa_goniometer()
  {
    enum_ <KappaGoniometer::Direction> ("KappaDirection")
      .value("PlusY", KappaGoniometer::PlusY)
      .value("PlusZ", KappaGoniometer::PlusZ)
      .value("MinusY", KappaGoniometer::MinusY)
      .value("MinusZ", KappaGoniometer::MinusZ)
      .export_values();

    enum_ <KappaGoniometer::ScanAxis> ("KappaScanAxis")
      .value("Omega", KappaGoniometer::Omega)
      .value("Phi", KappaGoniometer::Phi)
      .export_values();

    class_ <KappaGoniometer, bases <Goniometer> > ("KappaGoniometer")
      .def(init <double,
                 double,
                 double,
                 double,
                 KappaGoniometer::Direction,
                 KappaGoniometer::ScanAxis> ((
          arg("alpha"),
          arg("omega"),
          arg("kappa"),
          arg("phi"),
          arg("direction"),
          arg("scan_axis"))))
      .def("__init__", make_constructor(&make_kappa_goniometer))
      .def("get_alpha_angle",
        &KappaGoniometer::get_alpha_angle)
      .def("get_omega_angle",
        &KappaGoniometer::get_omega_angle)
      .def("get_kappa_angle",
        &KappaGoniometer::get_kappa_angle)
      .def("get_phi_angle",
        &KappaGoniometer::get_phi_angle)
      .def("get_direction",
        &KappaGoniometer::get_direction)
      .def("get_scan_axis",
        &KappaGoniometer::get_scan_axis)
      .def("get_omega_axis",
        &KappaGoniometer::get_omega_axis)
      .def("get_kappa_axis",
        &KappaGoniometer::get_kappa_axis)
      .def("get_phi_axis",
        &KappaGoniometer::get_phi_axis)
      .def("__str__", &kappa_goniometer_to_string)
      .def_pickle(KappaGoniometerPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
