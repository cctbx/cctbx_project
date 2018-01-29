/*
 * experiment.cc
 *
 *  Copyright (C) 2017 Diamond Light Source
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
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;


  struct ExperimentPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const Experiment &obj) {
      return boost::python::make_tuple(
        obj.get_beam(),
        obj.get_detector(),
        obj.get_goniometer(),
        obj.get_scan(),
        obj.get_crystal(),
        obj.get_profile(),
        obj.get_imageset(),
        obj.get_scaling_model());
    }
  };

  template <>
  boost::python::dict to_dict<Experiment>(const Experiment &obj) {
    boost::python::dict result;
    return result;
  }

  template <>
  Experiment* from_dict<Experiment>(boost::python::dict obj) {
    return new Experiment();
  }

  /**
   * Return function pointers to overrides for different types
   */
  struct experiment_contains_pointers {

    typedef bool (Experiment::*beam_type)(const boost::shared_ptr<BeamBase>&)const;
    typedef bool (Experiment::*detector_type)(const boost::shared_ptr<Detector>&)const;
    typedef bool (Experiment::*goniometer_type)(const boost::shared_ptr<Goniometer>&)const;
    typedef bool (Experiment::*scan_type)(const boost::shared_ptr<Scan>&)const;
    typedef bool (Experiment::*crystal_type)(const boost::shared_ptr<CrystalBase>&)const;
    typedef bool (Experiment::*object_type)(boost::python::object)const;

    static beam_type beam() {
      return &Experiment::contains;
    }

    static detector_type detector() {
      return &Experiment::contains;
    }

    static goniometer_type goniometer() {
      return &Experiment::contains;
    }

    static scan_type scan() {
      return &Experiment::contains;
    }

    static crystal_type crystal() {
      return &Experiment::contains;
    }

    static object_type object() {
      return &Experiment::contains;
    }

  };

  void export_experiment()
  {
    class_ <Experiment> ("Experiment")
      .def(init<
          boost::shared_ptr<BeamBase>,
          boost::shared_ptr<Detector>,
          boost::shared_ptr<Goniometer>,
          boost::shared_ptr<Scan>,
          boost::shared_ptr<CrystalBase>,
          boost::python::object,
          boost::python::object,
          boost::python::object>((
              arg("beam") = boost::shared_ptr<BeamBase>(),
              arg("detector") = boost::shared_ptr<Detector>(),
              arg("goniometer") = boost::shared_ptr<Goniometer>(),
              arg("scan") = boost::shared_ptr<Scan>(),
              arg("crystal") = boost::shared_ptr<CrystalBase>(),
              arg("profile") = boost::python::object(),
              arg("imageset") = boost::python::object(),
              arg("scaling_model") = boost::python::object())))
      .add_property(
          "beam",
          &Experiment::get_beam,
          &Experiment::set_beam)
      .add_property(
          "detector",
          &Experiment::get_detector,
          &Experiment::set_detector)
      .add_property(
          "goniometer",
          &Experiment::get_goniometer,
          &Experiment::set_goniometer)
      .add_property(
          "scan",
          &Experiment::get_scan,
          &Experiment::set_scan)
      .add_property(
          "crystal",
          &Experiment::get_crystal,
          &Experiment::set_crystal)
      .add_property(
          "profile",
          &Experiment::get_profile,
          &Experiment::set_profile)
      .add_property(
          "imageset",
          &Experiment::get_imageset,
          &Experiment::set_imageset)
      .add_property(
          "scaling_model",
          &Experiment::get_scaling_model,
          &Experiment::set_scaling_model)
      .def("__contains__", experiment_contains_pointers::beam())
      .def("__contains__", experiment_contains_pointers::detector())
      .def("__contains__", experiment_contains_pointers::goniometer())
      .def("__contains__", experiment_contains_pointers::scan())
      .def("__contains__", experiment_contains_pointers::crystal())
      .def("__contains__", experiment_contains_pointers::object())
      .def("__eq__", &Experiment::operator==)
      .def("is_consistent", &Experiment::is_consistent)
      .def("to_dict", &to_dict<Experiment>)
      .def("from_dict", &from_dict<Experiment>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(ExperimentPickleSuite());
  }

}}} // namespace = dxtbx::model::boost_python
