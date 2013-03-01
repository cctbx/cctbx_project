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
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/kappa_goniometer.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string goniometer_to_string(const Goniometer &goniometer) {
    boost::format fmt(
      "Goniometer:\n"
      "    rotation axis:  (%1%, %2%, %3%)\n"
      "    fixed rotation: [%4%, %5%, %6%,\n"
      "                     %7%, %8%, %9%,\n"
      "                     %10%, %11%, %12%]");        
    fmt % goniometer.get_rotation_axis()[0];
    fmt % goniometer.get_rotation_axis()[1];
    fmt % goniometer.get_rotation_axis()[2];
    fmt % goniometer.get_fixed_rotation()[0];
    fmt % goniometer.get_fixed_rotation()[1];
    fmt % goniometer.get_fixed_rotation()[2];
    fmt % goniometer.get_fixed_rotation()[3];
    fmt % goniometer.get_fixed_rotation()[4];
    fmt % goniometer.get_fixed_rotation()[5];
    fmt % goniometer.get_fixed_rotation()[6];
    fmt % goniometer.get_fixed_rotation()[7];
    fmt % goniometer.get_fixed_rotation()[8];
    return fmt.str();
  }

  std::string kappa_goniometer_direction_to_string(
      KappaGoniometer::Direction d) {
    if (d == KappaGoniometer::PlusY) {
      return "+y";
    } else if (d == KappaGoniometer::PlusZ) {
      return "+z";
    } else if (d == KappaGoniometer::MinusY) {
      return "-y";
    } else if (d == KappaGoniometer::MinusZ) {
      return "-z";
    }
    return "none";
  }

  std::string kappa_goniometer_scan_axis_to_string(
      KappaGoniometer::ScanAxis a) {
    if (a == KappaGoniometer::Omega) {
      return "omega";
    } else if (a == KappaGoniometer::Phi) {
      return "phi";
    }
    return "none";
  }

  std::string kappa_goniometer_to_string(const KappaGoniometer &goniometer) {
    boost::format fmt(
      "%1%\n"
      "    alpha angle:    %2%\n"
      "    omega angle:    %3%\n"
      "    kappa angle:    %4%\n"
      "    phi angle:      %5%\n"
      "    direction:      %6%\n"
      "    scan axis:      %7%\n"
      "    omega axis:     (%8%, %9%, %10%)\n"
      "    kappa axis:     (%11%, %12%, %13%)\n"
      "    phi axis:       (%14%, %15%, %16%)");
    fmt % goniometer_to_string(goniometer);
    fmt % goniometer.get_alpha_angle();
    fmt % goniometer.get_omega_angle();
    fmt % goniometer.get_kappa_angle();
    fmt % goniometer.get_phi_angle();
    fmt % kappa_goniometer_direction_to_string(goniometer.get_direction());
    fmt % kappa_goniometer_scan_axis_to_string(goniometer.get_scan_axis());
    fmt % goniometer.get_omega_axis()[0];
    fmt % goniometer.get_omega_axis()[1];
    fmt % goniometer.get_omega_axis()[2];
    fmt % goniometer.get_kappa_axis()[0];
    fmt % goniometer.get_kappa_axis()[1];
    fmt % goniometer.get_kappa_axis()[2];
    fmt % goniometer.get_phi_axis()[0];
    fmt % goniometer.get_phi_axis()[1];
    fmt % goniometer.get_phi_axis()[2];
    return fmt.str();
  }

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
    }

    KappaGoniometer::ScanAxis scan_axis = KappaGoniometer::NoAxis;
    boost::algorithm::to_lower(scan_axis_str);
    if (scan_axis_str == "omega") {
      scan_axis = KappaGoniometer::Omega;
    } else if (scan_axis_str == "phi") {
      scan_axis = KappaGoniometer::Phi;
    }

    return boost::shared_ptr<KappaGoniometer>(new KappaGoniometer(
      alpha, omega, kappa, phi, direction, scan_axis));
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
     
    enum_ <KappaGoniometer::Direction> ("KappaDirection")
      .value("PlusY", KappaGoniometer::PlusY)
      .value("PlusZ", KappaGoniometer::PlusZ)
      .value("MinusY", KappaGoniometer::MinusY)
      .value("MinusZ", KappaGoniometer::MinusZ);

    enum_ <KappaGoniometer::ScanAxis> ("KappaScanAxis")
      .value("Omega", KappaGoniometer::Omega)
      .value("Phi", KappaGoniometer::Phi);

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
      .def("__str__", &kappa_goniometer_to_string);
  }

}}} // namespace dxtbx::model::boost_python
