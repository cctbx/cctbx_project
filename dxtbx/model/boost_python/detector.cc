/*
 * detector.cc
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
#include <iostream>
#include <sstream>
#include <boost_adaptbx/std_pair_conversion.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/constants.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/boost_python/to_from_dict.h>
#include <dxtbx/model/boost_python/pickle_suite.h>

namespace dxtbx { namespace model { namespace boost_python {

  using scitbx::deg_as_rad;

  std::string detector_to_string(const Detector &detector) {
    std::stringstream ss;
    ss << detector;
    return ss.str();
  }

  static
  void detector_set_item(Detector &d, std::size_t i, const Panel &v) {
    d[i] = v;
  }

  static
  Panel* detector_get_item(Detector &d, std::size_t i) {
    return d.at(i);
  }

  static
  scitbx::af::shared<std::string> get_names(const Detector &d) {
    scitbx::af::shared<std::string> result(d.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = d[i].get_name();
    }
    return result;
  }

  static
  void rotate_around_origin(Detector &detector, vec3<double> axis, double angle, bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    detector.rotate_around_origin(axis, angle_rad);
  }

  template <>
  boost::python::dict to_dict<Detector::Node>(const Detector::Node &obj) {
    boost::python::dict result = to_dict<Panel>(obj);
    boost::python::list children;
    for (std::size_t i = 0; i < obj.size(); ++i) {
      Detector::const_node_pointer node = obj[i];
      if (node->is_panel()) {
        std::size_t idx = node->index();
        boost::python::dict panel;
        panel["panel"] = idx;
        children.append(panel);
      } else {
        children.append(to_dict(*node));
      }
    }
    result["children"] = children;
    return result;
  }


  void node_from_dict(
      boost::python::dict obj,
      Detector::Node *node,
      boost::python::list panels,
      scitbx::af::ref<bool> used) {
    if (obj.contains("panel")) {
      std::size_t idx = boost::python::extract<std::size_t>(obj["panel"]);
      DXTBX_ASSERT(used.size() == boost::python::len(panels));
      DXTBX_ASSERT(idx < boost::python::len(panels));
      DXTBX_ASSERT(used[idx] == false);
      used[idx] = true;
      boost::python::dict panel_dict = boost::python::extract<boost::python::dict>(panels[idx]);
      Panel *panel = from_dict<Panel>(panel_dict);
      node->add_panel(*panel, idx);
      delete panel;
    } else {
      Panel *panel = from_dict<Panel>(obj);
      node = node->add_group(*panel);
      for (std::size_t i = 0; i < boost::python::len(obj["children"]); ++i) {
        boost::python::dict child = boost::python::extract<boost::python::dict>(obj["children"][i]);
        node_from_dict(child, node, panels, used);
      }
      delete panel;
    }
  }


  static
  Detector detector_deepcopy(const Detector &detector, boost::python::object dict) {
    return Detector(detector);
  }

  template <>
  boost::python::dict to_dict<Detector>(const Detector &obj) {
    boost::python::dict result;
    boost::python::list panels;
    for (std::size_t i = 0; i < obj.size(); ++i) {
      panels.append(to_dict(obj[i]));
    }
    result["panels"] = panels;
    result["hierarchy"] = to_dict(*obj.root());
    return result;
  }

  Detector* detector_from_dict(Detector *result, boost::python::dict obj) {
    boost::python::list panels =
      boost::python::extract<boost::python::list>(obj["panels"]);
    if (obj.contains("hierarchy")) {
      boost::python::dict hierarchy = boost::python::extract<boost::python::dict>(obj["hierarchy"]);
      scitbx::af::shared<bool> used(boost::python::len(panels), false);
      DXTBX_ASSERT(!hierarchy.contains("panel"));
      Panel *panel = from_dict<Panel>(hierarchy);
	  // TODO: needs fixing for Python3
      // std::swap<Panel>(*result->root(), *panel);
      for (std::size_t i = 0; i < boost::python::len(hierarchy["children"]); ++i) {
        boost::python::dict child = boost::python::extract<boost::python::dict>(
            hierarchy["children"][i]);
        node_from_dict(child, result->root(), panels, used.ref());
      }
      delete panel;
      for (std::size_t i = 0; i < used.size(); ++i) {
        DXTBX_ASSERT(used[i] == true);
      }
    } else {
      for (std::size_t i = 0; i < boost::python::len(panels); ++i) {
        boost::python::dict panel_dict = boost::python::extract<boost::python::dict>(panels[i]);
        Panel *panel = from_dict<Panel>(panel_dict);
        result->add_panel(*panel);
        delete panel;
      }
    }
    return result;
  }

  template <>
  Detector* from_dict<Detector>(boost::python::dict obj) {
    Detector *result = new Detector();
    return detector_from_dict(result, obj);
  }

  struct DetectorPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(boost::python::object obj) {
      unsigned int version = 3;
      const Detector &detector = boost::python::extract<const Detector&>(obj);
      boost::python::dict data;

      // Convert panel array into pythoh list of panels
      boost::python::list panels;
      for (std::size_t i = 0; i < detector.size(); ++i) {
        panels.append(detector[i]);
      }
      data["panels"] = panels;

      // Convert hierarchy into dict
      data["hierarchy"] = to_dict(*detector.root());

      // Return the tuple
      return boost::python::make_tuple(version, data);
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      Detector *detector = boost::python::extract<Detector*>(obj);
      DXTBX_ASSERT(len(state) == 2);
      unsigned int version = boost::python::extract<unsigned int>(state[0]);
      DXTBX_ASSERT(version == 3);

      // Extract data from state object
      boost::python::dict data =
        boost::python::extract<boost::python::dict>(state[1]);
      boost::python::list panels =
        boost::python::extract<boost::python::list>(data["panels"]);
      boost::python::dict hierarchy =
        boost::python::extract<boost::python::dict>(data["hierarchy"]);

      DXTBX_ASSERT(!hierarchy.contains("panel"));
      Panel *panel = from_dict<Panel>(hierarchy);
	  // TODO: needs fixing for Python3
      //std::swap<Panel>(*detector->root(), *panel);
      copy_node(detector->root(), hierarchy, panels);
      delete panel;

      for (std::size_t i = 0; i < detector->size(); ++i) {
        DXTBX_ASSERT(detector->at(i) != NULL);
      }
    }

    static
    void copy_node(
          Detector::Node::pointer self,
          boost::python::dict node,
          boost::python::list panels) {
      for (std::size_t i = 0; i < boost::python::len(node["children"]); ++i) {
        boost::python::dict child =
          boost::python::extract<boost::python::dict>(node["children"][i]);
        if (child.contains("panel")) {
          std::size_t index = boost::python::extract<std::size_t>(child["panel"]);
          Panel panel = boost::python::extract<Panel>(panels[index]);
          self->add_panel(panel, index);
        } else {
          Panel *panel = from_dict<Panel>(child);
          copy_node(self->add_group(*panel), child, panels);
          delete panel;
        }
      }
    }

    static bool getstate_manages_dict() { return true; }
  };

  /* struct DetectorPickleSuite : boost::python::pickle_suite { */

  /*   static */
  /*   boost::python::tuple getstate(boost::python::object obj) { */
  /*     unsigned int version = 2; */
  /*     const Detector &detector = boost::python::extract<const Detector&>(obj); */
  /*     boost::python::dict data = to_dict<Detector>(detector); */
  /*     return boost::python::make_tuple( */
  /*       version, */
  /*       data); */
  /*   } */

  /*   static */
  /*   void setstate(boost::python::object obj, boost::python::tuple state) { */
  /*     Detector *detector = boost::python::extract<Detector*>(obj); */
  /*     DXTBX_ASSERT(len(state) == 2); */
  /*     unsigned int version = boost::python::extract<unsigned int>(state[0]); */
  /*     DXTBX_ASSERT(version == 2); */
  /*     boost::python::dict data = boost::python::extract<boost::python::dict>(state[1]); */
  /*     detector_from_dict(detector, data); */
  /*   } */

  /*   static bool getstate_manages_dict() { return true; } */
  /* }; */

  /* struct DetectorPickleSuite : boost::python::pickle_suite { */

  /*   static */
  /*   boost::python::dict getnode_state(const Detector::Node* root) { */
  /*     boost::python::dict result; */
  /*     result["panel"] = PanelPickleSuite::getstate(boost::python::object(root)); */
  /*     boost::python::list children; */
  /*     for (std::size_t i = 0; i < root->size(); ++i) { */
  /*       Detector::const_node_pointer node = (*root)[i]; */
  /*       if (node->is_panel()) { */
  /*         std::size_t idx = node->index(); */
  /*         boost::python::dict panel; */
  /*         panel["panel"] = idx; */
  /*         children.append(panel); */
  /*       } else { */
  /*         children.append(getnode_state(node)); */
  /*       } */
  /*     } */
  /*     result["children"] = children; */
  /*     return result; */
  /*   } */

  /*   static */
  /*   boost::python::tuple getstate(boost::python::object obj) { */
  /*     using namespace boost::python; */

  /*     unsigned int version = 2; */
  /*     const Detector &detector = boost::python::extract<const Detector&>(obj); */
  /*     boost::python::dict data; */
  /*     boost::python::list panels; */
  /*     for (std::size_t i = 0; i < boost::python::len(obj); ++i) { */
  /*       panels.append(detector[i]); */
  /*     } */
  /*     data["panels"] = panels; */
  /*     data["hierarchy"] = DetectorPickleSuite::getnode_state(detector.root()); */
  /*     return boost::python::make_tuple( */
  /*       version, */
  /*       data); */
  /*   } */

  /*   void setnode_state( */
  /*       boost::python::dict obj, */
  /*       Detector::Node *node, */
  /*       boost::python::list panels, */
  /*       scitbx::af::ref<bool> used) { */
  /*     if (obj.contains("panel")) { */
  /*       std::size_t idx = boost::python::extract<std::size_t>(obj["panel"]); */
  /*       DXTBX_ASSERT(used.size() == boost::python::len(panels)); */
  /*       DXTBX_ASSERT(idx < boost::python::len(panels)); */
  /*       DXTBX_ASSERT(used[idx] == false); */
  /*       used[idx] = true; */
  /*       boost::python::dict panel_dict = boost::python::extract<boost::python::dict>(panels[idx]); */
  /*       Panel *panel = PanelPickleSuite::setstate(panel_dict); */
  /*       node->add_panel(*panel, idx); */
  /*       delete panel; */
  /*     } else { */
  /*       Panel *panel = from_dict<Panel>(obj); */
  /*       node = node->add_group(*panel); */
  /*       for (std::size_t i = 0; i < boost::python::len(obj["children"]); ++i) { */
  /*         boost::python::dict child = boost::python::extract<boost::python::dict>(obj["children"][i]); */
  /*         node_from_dict(child, node, panels, used); */
  /*       } */
  /*       delete panel; */
  /*     } */
  /*   } */

  /*   static */
  /*   void setstate(boost::python::object obj, boost::python::tuple state) { */
  /*     using namespace boost::python; */
  /*     Detector *detector = boost::python::extract<Detector*>(obj); */
  /*     DXTBX_ASSERT(len(state) == 2); */
  /*     unsigned int version = extract<unsigned int>(state[0]); */
  /*     DXTBX_ASSERT(version == 2); */
  /*     boost::python::dict data = boost::python::extract<boost::python::dict>(state[1]); */
  /*     boost::python::list panels = boost::python::extract<boost::python::list>(data["panels"]); */
  /*     boost::python::dict hierarchy = boost::python::extract<boost::python::dict>(obj["hierarchy"]); */
  /*     scitbx::af::shared<bool> used(boost::python::len(panels), false); */
  /*     setnode_state(hierarchy, result->root(), panels, used.ref()); */
  /*     for (std::size_t i = 0; i < used.size(); ++i) { */
  /*       DXTBX_ASSERT(used[i] == true); */
  /*     } */
  /*   } */
  /* }; */

  void export_detector()
  {
    using namespace boost::python;

    class_ <Detector::Node, bases<Panel> >("DetectorNode", no_init)
      .def("add_group",
          (Detector::Node::pointer(Detector::Node::*)())
            &Detector::Node::add_group,
          return_internal_reference<>())
      .def("add_group",
          (Detector::Node::pointer(Detector::Node::*)(const Panel &))
            &Detector::Node::add_group,
          return_internal_reference<>())
      .def("add_panel",
          (Detector::Node::pointer(Detector::Node::*)())
            &Detector::Node::add_panel,
          return_internal_reference<>())
      .def("add_panel",
          (Detector::Node::pointer(Detector::Node::*)(const Panel &))
            &Detector::Node::add_panel,
          return_internal_reference<>())
      .def("parent",
          (Detector::Node::pointer(Detector::Node::*)())
            &Detector::Node::parent,
          return_internal_reference<>())
      .def("root",
          (Detector::Node::pointer(Detector::Node::*)())
            &Detector::Node::root,
          return_internal_reference<>())
      .def("__getitem__",
          (Detector::Node::pointer(Detector::Node::*)(std::size_t))
            &Detector::Node::operator[],
          return_internal_reference<>())
      .def("empty",
          &Detector::Node::empty)
      .def("__len__",
          &Detector::Node::size)
      .def("index",
          &Detector::Node::index)
      .def("is_panel",
          &Detector::Node::is_panel)
      .def("is_group",
          &Detector::Node::is_group)
      .def("set_frame",
          &Detector::Node::set_frame, (
            arg("fast_axis"),
            arg("slow_axis"),
            arg("origin")))
      .def("set_local_frame",
          &Detector::Node::set_local_frame, (
            arg("fast_axis"),
            arg("slow_axis"),
            arg("origin")))
      .def("set_parent_frame",
          &Detector::Node::set_parent_frame, (
            arg("fast_axis"),
            arg("slow_axis"),
            arg("origin")))
      .def("__eq__",
          &Detector::Node::operator==)
      .def("__ne__",
          &Detector::Node::operator!=)
      .def("is_similar_to",
          &Detector::Node::is_similar_to,(
            arg("rhs"),
            arg("fast_axis_tolerance")=1e-6,
            arg("slow_axis_tolerance")=1e-6,
            arg("origin_tolerance")=1e-6,
            arg("static_only")=false,
            arg("ignore_trusted_range")=false))
      .def("__iter__",
        iterator<Detector::Node, return_internal_reference<> >())
      .def("children",
        iterator<Detector::Node, return_internal_reference<> >())
      ;

    // Export a Detector base class
    class_ <Detector, boost::shared_ptr<Detector> > ("Detector")
      .def(init<const Panel&>())
      .def("hierarchy",
          (Detector::node_pointer(Detector::*)())&Detector::root,
          return_internal_reference<>())
      .def("add_group",
          (Detector::node_pointer(Detector::*)())
          &Detector::add_group, return_internal_reference<>())
      .def("add_group",
          (Detector::node_pointer(Detector::*)(const Panel&))
          &Detector::add_group, return_internal_reference<>())
      .def("add_panel",
        (Detector::node_pointer(Detector::*)())&Detector::add_panel,
        return_internal_reference<>())
      .def("add_panel",
        (Detector::node_pointer(Detector::*)(const Panel&))&Detector::add_panel,
        return_internal_reference<>())
      .def("__len__",
        &Detector::size)
      .def("__setitem__",
        &detector_set_item)
      .def("__getitem__",
        &detector_get_item,
        return_internal_reference<>())
      .def("__iter__",
        iterator<Detector, return_internal_reference<> >())
      .def("__eq__", &Detector::operator==)
      .def("__ne__", &Detector::operator!=)
      .def("is_similar_to", &Detector::is_similar_to, (
            arg("other"),
            arg("fast_axis_tolerance")=1e-6,
            arg("slow_axis_tolerance")=1e-6,
            arg("origin_tolerance")=1e-6,
            arg("static_only")=false,
            arg("ignore_trusted_range")=false))
      .def("get_max_resolution",
        &Detector::get_max_resolution, (arg("s0")))
      .def("get_max_inscribed_resolution",
        &Detector::get_max_inscribed_resolution, (arg("s0")))
      .def("get_ray_intersection",
        &Detector::get_ray_intersection, (arg("s1")))
      .def("get_panel_intersection",
        &Detector::get_panel_intersection, (arg("s1")))
      //.def("do_panels_intersect",
      //  &Detector::do_panels_intersect)
      .def("get_names", &get_names)
      .def("rotate_around_origin",
          &rotate_around_origin, (
            arg("axis"),
            arg("angle"),
            arg("deg")=true))
      .def("__str__", &detector_to_string)
      .def("__deepcopy__", &detector_deepcopy)
      .def("__copy__", &detector_deepcopy)
      .def("to_dict", &to_dict<Detector>)
      .def("from_dict", &from_dict<Detector>,
        return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(DetectorPickleSuite())
      ;

    boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, vec2<double> >();
  }

}}} // namespace dxtbx::model::boost_python
