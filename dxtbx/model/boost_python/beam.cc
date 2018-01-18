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
#include <scitbx/constants.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/boost_python/to_from_dict.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  std::string beam_to_string(const Beam &beam) {
    std::stringstream ss;
    ss << beam;
    return ss.str();
  }

  struct BeamPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Beam &obj) {
      return boost::python::make_tuple(
        obj.get_direction(),
        obj.get_wavelength(),
        obj.get_divergence(),
        obj.get_sigma_divergence(),
        obj.get_polarization_normal(),
        obj.get_polarization_fraction(),
        obj.get_flux(),
        obj.get_transmission());
    }
  };

  template <>
  boost::python::dict to_dict<Beam>(const Beam &obj) {
    boost::python::dict result;
    result["direction"] = obj.get_direction();
    result["wavelength"] = obj.get_wavelength();
    result["divergence"] = rad_as_deg(obj.get_divergence());
    result["sigma_divergence"] = rad_as_deg(obj.get_sigma_divergence());
    result["polarization_normal"] = obj.get_polarization_normal();
    result["polarization_fraction"] = obj.get_polarization_fraction();
    result["flux"] = obj.get_flux();
    result["transmission"] = obj.get_transmission();
    return result;
  }

  template <>
  Beam* from_dict<Beam>(boost::python::dict obj) {
    return new Beam(
      boost::python::extract< vec3<double> >(obj["direction"]),
      boost::python::extract< double >(obj["wavelength"]),
      deg_as_rad(
        boost::python::extract< double >(obj.get("divergence", 0.0))),
      deg_as_rad(
        boost::python::extract< double >(obj.get("sigma_divergence", 0.0))),
      boost::python::extract< vec3<double> >(
        obj.get("polarization_normal", vec3<double>(0.0, 1.0, 0.0))),
      boost::python::extract< double >(obj.get("polarization_fraction", 0.999)),
      boost::python::extract< double >(obj.get("flux", 0)),
      boost::python::extract< double >(obj.get("transmission", 1)));
  }

  static Beam* make_beam(vec3<double> sample_to_source, double wavelength,
                         double divergence, double sigma_divergence, bool deg) {
    Beam *beam = NULL;
    if (deg) {
      beam = new Beam(sample_to_source, wavelength,
                      deg_as_rad(divergence),
                      deg_as_rad(sigma_divergence));
    } else {
      beam = new Beam(sample_to_source, wavelength,
                      divergence, sigma_divergence);
    }
    return beam;
  }

  static Beam* make_beam_w_s0(vec3<double> s0, double divergence,
                              double sigma_divergence, bool deg) {
    Beam *beam = NULL;
    if (deg) {
      beam = new Beam(s0, deg_as_rad(divergence),
                      deg_as_rad(sigma_divergence));
    } else {
      beam = new Beam(s0, divergence, sigma_divergence);
    }
    return beam;
  }

  static Beam* make_beam_w_all(vec3<double> sample_to_source,
                               double wavelength,
                               double divergence, double sigma_divergence,
                               vec3<double> polarization_normal,
                               double polarization_fraction,
                               double flux,
                               double transmission,
                               bool deg) {
    Beam *beam = NULL;
    if (deg) {
      beam = new Beam(sample_to_source, wavelength,
                      deg_as_rad(divergence),
                      deg_as_rad(sigma_divergence),
                      polarization_normal,
                      polarization_fraction,
                      flux,
                      transmission);
    } else {
      beam = new Beam(sample_to_source, wavelength,
                      divergence, sigma_divergence,
                      polarization_normal,
                      polarization_fraction,
                      flux,
                      transmission);
    }
    return beam;
  }

  static
  double get_divergence(const Beam &beam, bool deg) {
    double divergence = beam.get_divergence();
    return deg ? rad_as_deg(divergence) : divergence;
  }

  static
  double get_sigma_divergence(const Beam &beam, bool deg) {
    double sigma_divergence = beam.get_sigma_divergence();
    return deg ? rad_as_deg(sigma_divergence) : sigma_divergence;
  }

  static
  void set_divergence(Beam &beam, double divergence, bool deg) {
    beam.set_divergence(deg ? deg_as_rad(divergence) : divergence);
  }

  static
  void set_sigma_divergence(Beam &beam, double sigma_divergence,
                            bool deg) {
    beam.set_sigma_divergence(
      deg ? deg_as_rad(sigma_divergence) : sigma_divergence);
  }

  static
  void rotate_around_origin(Beam &beam, vec3<double> axis, double angle, bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    beam.rotate_around_origin(axis, angle_rad);
  }

  static
  void Beam_set_s0_at_scan_points_from_tuple(Beam &beam, boost::python::tuple l) {
    scitbx::af::shared< vec3<double> > s0_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      vec3<double> s0 = boost::python::extract< vec3<double> >(l[i]);
      s0_list.push_back(s0);
    }
    beam.set_s0_at_scan_points(s0_list.const_ref());
  }

  static
  void Beam_set_s0_at_scan_points_from_list(Beam &beam, boost::python::list l) {
    scitbx::af::shared< vec3<double> > s0_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      vec3<double> s0 = boost::python::extract< vec3<double> >(l[i]);
      s0_list.push_back(s0);
    }
    beam.set_s0_at_scan_points(s0_list.const_ref());
  }

  void export_beam()
  {
    // Export BeamBase
    class_ <BeamBase, boost::noncopyable> ("BeamBase", no_init)
      .def("get_direction",
        &BeamBase::get_direction)
      .def("set_direction",
        &BeamBase::set_direction)
      .def("get_wavelength",
        &BeamBase::get_wavelength)
      .def("set_wavelength",
        &BeamBase::set_wavelength)
      .def("get_s0",
        &BeamBase::get_s0)
      .def("set_s0",
        &BeamBase::set_s0)
      .def("get_unit_s0",
        &BeamBase::get_unit_s0)
      .def("set_unit_s0",
        &BeamBase::set_unit_s0)
      .def("get_divergence",
        &get_divergence, (
          arg("deg") = true))
      .def("set_divergence",
        &set_divergence, (
          arg("divergence"),
          arg("deg") = true))
      .def("get_sigma_divergence",
        &get_sigma_divergence, (
          arg("deg") = true))
      .def("set_sigma_divergence",
        &set_sigma_divergence, (
          arg("sigma_divergence"),
          arg("deg") = true))
      .def("get_polarization_normal",
        &BeamBase::get_polarization_normal)
      .def("set_polarization_normal",
        &BeamBase::set_polarization_normal)
      .def("get_polarization_fraction",
        &BeamBase::get_polarization_fraction)
      .def("set_polarization_fraction",
        &BeamBase::set_polarization_fraction)
      .def("get_flux",
        &BeamBase::get_flux)
      .def("set_flux",
        &BeamBase::set_flux)
      .def("get_transmission",
        &BeamBase::get_transmission)
      .def("set_transmission",
        &BeamBase::set_transmission)
      .add_property("num_scan_points", &BeamBase::get_num_scan_points)
      .def("get_num_scan_points",
        &BeamBase::get_num_scan_points)
      .def("set_s0_at_scan_points",
        &BeamBase::set_s0_at_scan_points)
      .def("set_s0_at_scan_points",
        &Beam_set_s0_at_scan_points_from_tuple)
      .def("set_s0_at_scan_points",
        &Beam_set_s0_at_scan_points_from_list)
      .def("get_s0_at_scan_points",
        &BeamBase::get_s0_at_scan_points)
      .def("get_s0_at_scan_point",
        &BeamBase::get_s0_at_scan_point)
      .def("reset_scan_points",
        &BeamBase::reset_scan_points)
      .def("rotate_around_origin",
          &rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("__eq__", &BeamBase::operator==)
      .def("__ne__", &BeamBase::operator!=)
      .def("is_similar_to", &BeamBase::is_similar_to, (
            arg("other"),
            arg("wavelength_tolerance")=1e-6,
            arg("direction_tolerance")=1e-6,
            arg("polarization_normal_tolerance")=1e-6,
            arg("polarization_fraction_tolerance")=1e-6));

    // Export Beam : BeamBase
    class_ <Beam, boost::shared_ptr<Beam>, bases <BeamBase> > ("Beam")
      .def(init<const Beam&>())
      .def(init <vec3 <double>,
                 double> ((
          arg("direction"),
          arg("wavelength"))))
      .def(init <vec3 <double> > ((
          arg("s0"))))
      .def("__init__",
          make_constructor(
          &make_beam,
          default_call_policies(), (
            arg("direction"),
            arg("wavelength"),
            arg("divergence"),
            arg("sigma_divergence"),
            arg("deg") = true)))
      .def("__init__",
          make_constructor(
          &make_beam_w_s0,
          default_call_policies(), (
            arg("s0"),
            arg("divergence"),
            arg("sigma_divergence"),
            arg("deg") = true)))
      .def("__init__",
          make_constructor(
          &make_beam_w_all,
          default_call_policies(), (
            arg("direction"),
            arg("wavelength"),
            arg("divergence"),
            arg("sigma_divergence"),
            arg("polarization_normal"),
            arg("polarization_fraction"),
            arg("flux"),
            arg("transmission"),
            arg("deg") = true)))
      .def("__str__", &beam_to_string)
      .def("to_dict", &to_dict<Beam>)
      .def("from_dict", &from_dict<Beam>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(BeamPickleSuite());

    scitbx::af::boost_python::flex_wrapper<Beam>::plain("flex_Beam");
  }

}}} // namespace = dxtbx::model::boost_python
