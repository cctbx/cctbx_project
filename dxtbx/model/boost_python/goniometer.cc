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
#include <scitbx/constants.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/kappa_goniometer.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;

  static
  void rotate_around_origin(Goniometer &goniometer, vec3<double> axis, double angle, bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    goniometer.rotate_around_origin(axis, angle_rad);
  }

  std::string goniometer_to_string(const Goniometer &goniometer) {
    std::stringstream ss;
    ss << goniometer;
    return ss.str();
  }

  struct GoniometerPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Goniometer &obj) {
      return boost::python::make_tuple(
        obj.get_rotation_axis(),
        obj.get_fixed_rotation(),
        obj.get_setting_rotation());
    }
  };

  template <>
  boost::python::dict to_dict<Goniometer>(const Goniometer &obj) {
    boost::python::dict result;
    result["rotation_axis"] = obj.get_rotation_axis();
    result["fixed_rotation"] = obj.get_fixed_rotation();
    result["setting_rotation"] = obj.get_setting_rotation();
    return result;
  }

  template <>
  Goniometer* from_dict<Goniometer>(boost::python::dict obj) {
    return new Goniometer(
      boost::python::extract< vec3<double> >(obj["rotation_axis"]),
      boost::python::extract< mat3<double> >(obj.get("fixed_rotation",
        mat3<double>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))),
      boost::python::extract< mat3<double> >(obj.get("setting_rotation",
        mat3<double>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))));
  }

  void export_goniometer()
  {
    class_ <GoniometerBase> ("GoniometerBase");

    class_ <Goniometer, bases <GoniometerBase> > ("Goniometer")
      .def(init<const Goniometer&>())
      .def(init <vec3 <double> > ((
          arg("rotation_axis"))))
      .def(init <vec3 <double>,
                 mat3 <double> > ((
          arg("rotation_axis"),
          arg("fixed_rotation_matrix"))))
      .def(init <vec3 <double>,
                 mat3 <double>,
                 mat3 <double> > ((
          arg("rotation_axis"),
          arg("fixed_rotation_matrix"),
          arg("setting_rotation_matrix"))))
      .def("get_rotation_axis",
        &Goniometer::get_rotation_axis)
      .def("set_rotation_axis",
        &Goniometer::set_rotation_axis)
      .def("get_fixed_rotation",
        &Goniometer::get_fixed_rotation)
      .def("set_fixed_rotation",
        &Goniometer::set_fixed_rotation)
      .def("get_setting_rotation",
        &Goniometer::get_setting_rotation)
      .def("set_setting_rotation",
        &Goniometer::set_setting_rotation)
      .def("rotate_around_origin",
          &rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("__eq__", &Goniometer::operator==)
      .def("__ne__", &Goniometer::operator!=)
      .def("is_similar_to", &Goniometer::is_similar_to, (
            arg("other"),
            arg("rotation_axis_tolerance")=1e-6,
            arg("fixed_rotation_tolerance")=1e-6,
            arg("setting_rotation_tolerance")=1e-6))
      .def("__str__", &goniometer_to_string)
      .def("to_dict", &to_dict<Goniometer>)
      .def("from_dict", &from_dict<Goniometer>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(GoniometerPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
