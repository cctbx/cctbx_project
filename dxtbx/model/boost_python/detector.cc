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
#include <boost/format.hpp>
#include <string>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <dxtbx/model/detector.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string flat_panel_detector_to_string(const FlatPanelDetector &detector)
  {
    boost::format fmt(
      "Detector:\n"
      "    type:          %1%\n"
      "    fast axis:     (%2%, %3%, %4%)\n"
      "    slow axis:     (%5%, %6%, %7%)\n"      
      "    normal:        (%8%, %9%, %10%)\n"
      "    origin:        (%11%, %12%)\n"
      "    pixel_size:    (%13%, %14%)\n"
      "    image_size:    (%15%, %16%)\n"
      "    trusted_range: (%17, %18%)");
              
    fmt % detector.get_type();
    fmt % detector.get_fast_axis()[0];
    fmt % detector.get_fast_axis()[1];
    fmt % detector.get_fast_axis()[2];
    fmt % detector.get_slow_axis()[0];
    fmt % detector.get_slow_axis()[1];
    fmt % detector.get_slow_axis()[2];
    fmt % detector.get_normal()[0];
    fmt % detector.get_normal()[1];
    fmt % detector.get_normal()[2];
    fmt % detector.get_origin()[0];
    fmt % detector.get_origin()[1];
    fmt % detector.get_pixel_size()[0];
    fmt % detector.get_pixel_size()[1];
    fmt % detector.get_image_size()[0];
    fmt % detector.get_image_size()[1];
    fmt % detector.get_trusted_range()[0];
    fmt % detector.get_trusted_range()[1];
    return fmt.str();
  }

  std::string multi_flat_panel_detector_to_string(
      const MultiFlatPanelDetector &detector)
  {
    std::string str;
    str += "MultiFlatPanelDetector:\n";
    for (std::size_t i = 0; i < detector.num_panels(); ++i) {
      str += "\n";
      str += flat_panel_detector_to_string(detector[i]);
    }
    return str;
  }

  void multi_flat_panel_detector_set_item(MultiFlatPanelDetector &d, 
      std::size_t i, const FlatPanelDetector &v) {
    d[i] = v;
  }

  void multi_flat_panel_detector_del_item(MultiFlatPanelDetector &d, 
      std::size_t i) {
    d.remove_panel(i);
  }

  FlatPanelDetector& multi_flat_panel_detector_get_item(
      MultiFlatPanelDetector &d, std::size_t i) {
    return d[i];
  }

  void export_detector() 
  {
    // Export a flex array - should probably move somewhere else
    scitbx::af::boost_python::flex_wrapper <int4>::plain("flex_int4");

    // Export the DetectorBase class
    class_ <DetectorBase> ("DetectorBase");

    // Export the FlatPanelDetector class
    class_ <FlatPanelDetector, bases <DetectorBase> > ("FlatPanelDetector")
      .def(init <std::string,
                 vec3 <double>,
                 vec3 <double>,
                 vec3 <double>,
                 vec2 <double>,
                 vec2 <std::size_t>,
                vec2 <int> > ((                 
          arg("type"),
          arg("fast_axis"),
          arg("slow_axis"),
          arg("origin"),
          arg("pixel_size"),
          arg("image_size"),
           arg("trusted_range"))))
      .add_property("type",
        &FlatPanelDetector::get_type,
        &FlatPanelDetector::set_type)    
      .add_property("fast_axis",
        &FlatPanelDetector::get_fast_axis,
        &FlatPanelDetector::set_fast_axis)
      .add_property("slow_axis",
        &FlatPanelDetector::get_slow_axis,
        &FlatPanelDetector::set_slow_axis)
      .add_property("normal",
        &FlatPanelDetector::get_normal,
        &FlatPanelDetector::set_normal)
      .add_property("origin",
        &FlatPanelDetector::get_origin,
        &FlatPanelDetector::set_origin)
      .add_property("pixel_size",
        &FlatPanelDetector::get_pixel_size,
        &FlatPanelDetector::set_pixel_size)
      .add_property("image_size",
        &FlatPanelDetector::get_image_size,
        &FlatPanelDetector::set_image_size)
      .add_property("trusted_range",
        &FlatPanelDetector::get_trusted_range,
        &FlatPanelDetector::set_trusted_range)
      .add_property("d_matrix",
        &FlatPanelDetector::get_d_matrix,
        &FlatPanelDetector::set_d_matrix)
      .add_property("inverse_d_matrix",
        &FlatPanelDetector::get_inverse_d_matrix,
        &FlatPanelDetector::set_inverse_d_matrix)
      .def("__eq__", &FlatPanelDetector::operator==)
      .def("__ne__", &FlatPanelDetector::operator!=)
      .def("__str__", &flat_panel_detector_to_string);

    // Export a MultiFlatPanelDetector class
    class_ <MultiFlatPanelDetector, 
            bases <DetectorBase> > ("MultiFlatPanelDetector")
      .def(init <std::string> ((
          arg("type"))))
      .def("add_panel",
        &MultiFlatPanelDetector::add_panel, (
          arg("panel")))
      .def("num_panels",
        &MultiFlatPanelDetector::num_panels)
      .def("__len__", 
        &MultiFlatPanelDetector::num_panels)
      .def("__setitem__", 
        &multi_flat_panel_detector_set_item)
      .def("__delitem__", 
        &multi_flat_panel_detector_del_item)
      .def("__getitem__", 
        &multi_flat_panel_detector_get_item, 
        return_internal_reference <> ())
      .def("__iter__", 
        iterator <
          MultiFlatPanelDetector, 
          return_internal_reference<> >())
      .def("__eq__", &MultiFlatPanelDetector::operator==)
      .def("__ne__", &MultiFlatPanelDetector::operator!=)
      .def("__str__", &multi_flat_panel_detector_to_string);
  }

}}} // namespace dials::model::boost_python
