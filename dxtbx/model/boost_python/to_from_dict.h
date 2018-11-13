/*
 * to_from_dict.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_BOOST_PYTHON_TO_FROM_DICT_H
#define DXTBX_MODEL_BOOST_PYTHON_TO_FROM_DICT_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/format/image.h>

namespace dxtbx { namespace model { namespace boost_python {

  using format::Image;

  template <typename T>
  boost::python::dict to_dict(const T &obj);

  template <typename T>
  T* from_dict(boost::python::dict obj);

  template <>
  boost::python::dict to_dict<BeamBase>(const BeamBase &obj);

  template <>
  boost::python::dict to_dict<Goniometer>(const Goniometer &obj);

  template <>
  boost::python::dict to_dict<VirtualPanel>(const VirtualPanel &obj);

  template <>
  boost::python::dict to_dict<Panel>(const Panel &obj);

  template <>
  boost::python::dict to_dict<Detector>(const Detector &obj);

  template <>
  boost::python::dict to_dict<Scan>(const Scan &obj);

  template <>
  boost::python::dict to_dict<CrystalBase>(const CrystalBase &obj);

  template <>
  BeamBase* from_dict<BeamBase>(boost::python::dict obj);

  template <>
  Goniometer* from_dict<Goniometer>(boost::python::dict obj);

  template <>
  VirtualPanel* from_dict<VirtualPanel>(boost::python::dict obj);

  template <>
  Panel* from_dict<Panel>(boost::python::dict obj);

  Panel* panel_from_dict_with_offset(boost::python::dict obj,
      scitbx::af::versa<double, scitbx::af::c_grid<2> > dx,
      scitbx::af::versa<double, scitbx::af::c_grid<2> > dy);

  template <>
  Detector* from_dict<Detector>(boost::python::dict obj);

  Detector* detector_from_dict_with_offset(
        boost::python::dict obj,
        const Image<double> &dx,
        const Image<double> &dy);

  template <>
  Scan* from_dict<Scan>(boost::python::dict obj);

  template <>
  CrystalBase* from_dict<CrystalBase>(boost::python::dict obj);

  template <typename K, typename V>
  boost::python::dict MaptoPythonDict(std::map<K, V> map) {
    typename std::map<K, V>::iterator iter;
      boost::python::dict dictionary;
      for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = boost::python::list(iter->second);
      }
      return dictionary;
  }

}}} // namespace dxtbx::model::boost_python

#endif /* DXTBX_MODEL_BOOST_PYTHON_TO_FROM_DICT_H */
