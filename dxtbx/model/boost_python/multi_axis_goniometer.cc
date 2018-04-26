/*
 * multi_axis_goniometer.cc
 *
 *  Copyright (C) 2016 Diamond Light Source
 *
 *  Author: Richard Gildea
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
#include <dxtbx/model/multi_axis_goniometer.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string multi_axis_goniometer_to_string(const MultiAxisGoniometer &goniometer) {
    std::stringstream ss;
    ss << goniometer;
    return ss.str();
  }

  struct MultiAxisGoniometerPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const MultiAxisGoniometer &obj) {
      return boost::python::make_tuple(
        obj.get_axes(),
        obj.get_angles(),
        obj.get_names(),
        obj.get_scan_axis());
    }

    static
    boost::python::tuple getstate(boost::python::object obj)
    {
      const MultiAxisGoniometer &goniometer = boost::python::extract<const MultiAxisGoniometer &>(obj)();
      return boost::python::make_tuple(
          obj.attr("__dict__"),
          goniometer.get_setting_rotation_at_scan_points());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state)
    {
      MultiAxisGoniometer &goniometer = boost::python::extract<MultiAxisGoniometer&>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 2);

      // restore the object's __dict__
      boost::python::dict d = boost::python::extract<boost::python::dict>(
          obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref< mat3<double> > S_list = boost::python::extract<
        scitbx::af::const_ref< mat3<double> > >(state[1]);
      goniometer.set_setting_rotation_at_scan_points(S_list);
    }

    static bool getstate_manages_dict() { return true; }
  };

  boost::python::dict to_dict(const MultiAxisGoniometer &obj) {
    boost::python::dict result;
    result["axes"] = boost::python::list(obj.get_axes());
    result["angles"] = boost::python::list(obj.get_angles());
    result["names"] = boost::python::list(obj.get_names());
    result["scan_axis"] = obj.get_scan_axis();
    if(obj.get_num_scan_points() > 0){
      boost::python::list l;
      scitbx::af::shared< mat3<double> > setting_rotation_at_scan_points =
          obj.get_setting_rotation_at_scan_points();
      for (scitbx::af::shared< mat3<double> >::iterator it = setting_rotation_at_scan_points.begin();
           it != setting_rotation_at_scan_points.end();
           ++it) {
        l.append(boost::python::tuple(*it));
      }
      result["setting_rotation_at_scan_points"] = l;
    }
    return result;
    return result;
  };

  MultiAxisGoniometer* from_dict(boost::python::dict obj) {
    scitbx::af::shared<vec3<double> > axes =
      boost::python::extract< scitbx::af::shared<vec3<double> > >(obj["axes"]);
    scitbx::af::shared<double> angles =
      boost::python::extract< scitbx::af::shared<double> >(obj["angles"]);
    scitbx::af::shared<std::string> names =
      boost::python::extract< scitbx::af::shared<std::string> >(obj["names"]);
    MultiAxisGoniometer* g = new MultiAxisGoniometer(
      axes.const_ref(), angles.const_ref(), names.const_ref(),
      boost::python::extract< std::size_t >(obj["scan_axis"]));
    if(obj.has_key("setting_rotation_at_scan_points")){
      boost::python::list l = boost::python::extract<boost::python::list>(obj["setting_rotation_at_scan_points"]);
      scitbx::af::shared< mat3<double> > S_array;
      for (std::size_t i = 0; i < boost::python::len(l); ++i) {
        mat3<double> S = boost::python::extract< mat3<double> >(l[i]);
        S_array.push_back(S);
      }
      g->set_setting_rotation_at_scan_points(S_array.const_ref());
    }
    return g;
  };

  static boost::shared_ptr<MultiAxisGoniometer> make_multi_axis_goniometer(
    const scitbx::af::const_ref<vec3<double> > &axes,
    const scitbx::af::const_ref<double> &angles,
    const scitbx::af::const_ref<std::string> &names,
    std::size_t scan_axis)
  {
    return boost::shared_ptr<MultiAxisGoniometer>(new MultiAxisGoniometer(
      axes, angles, names, scan_axis));
  }

  void export_multi_axis_goniometer()
  {

    class_ <MultiAxisGoniometer, bases <Goniometer> > ("MultiAxisGoniometer")
      .def(init <const scitbx::af::const_ref<vec3<double> > &,
                 const scitbx::af::const_ref<double> &,
                 const scitbx::af::const_ref<std::string> &,
                 std::size_t> ((
          arg("axes"),
          arg("angles"),
          arg("names"),
          arg("scan_axis"))))
      .def("__init__", make_constructor(&make_multi_axis_goniometer))
      .def("get_axes",
        &MultiAxisGoniometer::get_axes)
      .def("set_axes",
        &MultiAxisGoniometer::set_axes)
      .def("get_angles",
        &MultiAxisGoniometer::get_angles)
      .def("set_angles",
        &MultiAxisGoniometer::set_angles)
      .def("get_names",
        &MultiAxisGoniometer::get_names)
      .def("get_scan_axis",
        &MultiAxisGoniometer::get_scan_axis)
      .def("__str__", &multi_axis_goniometer_to_string)
      .def("to_dict", &to_dict)
      .def("from_dict", &from_dict,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(MultiAxisGoniometerPickleSuite());
  }

}}} // namespace dxtbx::model::boost_python
