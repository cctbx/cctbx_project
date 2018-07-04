/*
 * pickle_suite.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_BOOST_PYTHON_PICKLE_SUITE_H
#define DXTBX_MODEL_BOOST_PYTHON_PICKLE_SUITE_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/scan.h>

namespace dxtbx { namespace model { namespace boost_python {

  struct VirtualPanelPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(boost::python::object obj) {
      using namespace boost::python;

      unsigned int version = 1;
      const VirtualPanel &p = extract<const VirtualPanel&>(obj)();
      return boost::python::make_tuple(
        version,
        obj.attr("__dict__"),
        p.get_name(),
        p.get_type(),
        p.get_local_fast_axis(),
        p.get_local_slow_axis(),
        p.get_local_origin(),
        p.get_parent_fast_axis(),
        p.get_parent_slow_axis(),
        p.get_parent_origin());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      using namespace boost::python;

      DXTBX_ASSERT(len(state) == 10);
      unsigned int version = extract<unsigned int>(state[0]);
      DXTBX_ASSERT(version == 1);
      extract<dict>(obj.attr("__dict__"))().update(state[1]);
      VirtualPanel &p = extract<VirtualPanel&>(obj)();
      p.set_name(extract<std::string>(state[2]));
      p.set_type(extract<std::string>(state[3]));
      p.set_local_frame(
        extract< vec3<double> >(state[4]),
        extract< vec3<double> >(state[5]),
        extract< vec3<double> >(state[6]));
      p.set_parent_frame(
        extract< vec3<double> >(state[7]),
        extract< vec3<double> >(state[8]),
        extract< vec3<double> >(state[9]));
    }

    static bool getstate_manages_dict() { return true; }
  };

  struct PanelPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(boost::python::object obj) {
      using namespace boost::python;

      unsigned int version = 2;
      const Panel &p = extract<const Panel&>(obj)();
      boost::python::dict data;
      data["name"] = p.get_name();
      data["type"] = p.get_type();
      data["fast_axis"] = p.get_local_fast_axis();
      data["slow_axis"] = p.get_local_slow_axis();
      data["origin"] = p.get_local_origin();
      data["parent_fast_axis"] = p.get_parent_fast_axis();
      data["parent_slow_axis"] = p.get_parent_slow_axis();
      data["parent_origin"] = p.get_parent_origin();
      data["raw_image_offset"] = p.get_raw_image_offset();
      data["image_size"] = p.get_image_size();
      data["pixel_size"] = p.get_pixel_size();
      data["trusted_range"] = p.get_trusted_range();
      data["thickness"] = p.get_thickness();
      data["gain"] = p.get_gain();
      data["material"] = p.get_material();
      data["identifier"] = p.get_identifier();
      data["mu"] = p.get_mu();
      data["mask"] = boost::python::list(p.get_mask());
      data["px_mm_strategy"] = p.get_px_mm_strategy();
      return boost::python::make_tuple(
        version,
        obj.attr("__dict__"),
        data);
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      using namespace boost::python;

      DXTBX_ASSERT(len(state) == 3);
      unsigned int version = extract<unsigned int>(state[0]);
      DXTBX_ASSERT(version == 2);
      extract<dict>(obj.attr("__dict__"))().update(state[1]);
      boost::python::dict data = extract<dict>(state[2]);
      Panel &p = extract<Panel&>(obj)();
      p.set_name(extract<std::string>(data["name"]));
      p.set_type(extract<std::string>(data["type"]));
      p.set_local_frame(
        extract< vec3<double> >(data["fast_axis"]),
        extract< vec3<double> >(data["slow_axis"]),
        extract< vec3<double> >(data["origin"]));
      p.set_parent_frame(
        extract< vec3<double> >(data["parent_fast_axis"]),
        extract< vec3<double> >(data["parent_slow_axis"]),
        extract< vec3<double> >(data["parent_origin"]));
      if (data.has_key("pixel_size")) {
        p.set_pixel_size(extract< tiny<double,2> >(data["pixel_size"]));
      }
      if (data.has_key("raw_image_offset")) {
        p.set_raw_image_offset(extract<int2>(data["raw_image_offset"]));
      }
      if (data.has_key("image_size")) {
        p.set_image_size(extract< tiny<std::size_t,2> >(data["image_size"]));
      }
      if (data.has_key("trusted_range")) {
        p.set_trusted_range(extract< tiny<double,2> >(data["trusted_range"]));
      }
      if (data.has_key("thickness")) {
        p.set_thickness(extract<double>(data["thickness"]));
      }
      if (data.has_key("gain")) {
        p.set_gain(extract<double>(data["gain"]));
      }
      if (data.has_key("material")) {
        p.set_material(extract<std::string>(data["material"]));
      }
      if (data.has_key("identifier")) {
        p.set_identifier(extract<std::string>(data["identifier"]));
      }
      if (data.has_key("mu")) {
        p.set_mu(extract<double>(data["mu"]));
      }
      if (data.has_key("px_mm_strategy")) {
        p.set_px_mm_strategy(extract< shared_ptr<PxMmStrategy> >(data["px_mm_strategy"]));
      }
      if (data.has_key("mask")) {
        scitbx::af::shared<int4> mask =
          boost::python::extract< scitbx::af::shared<int4> >(
            boost::python::extract<boost::python::list>(data["mask"]));
        p.set_mask(mask.const_ref());
      }
    }

    static bool getstate_manages_dict() { return true; }
  };


}}} // namespace dxtbx::model::boost_python

#endif /* DXTBX_MODEL_BOOST_PYTHON_PICKLE_SUITE_H */
