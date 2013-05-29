/*
 * parallax_correction.cc
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
#include <dxtbx/model/parallax_correction.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_parallax_correction()
  {
    def("parallax_correction", &parallax_correction, (
      arg("d"), arg("la"), arg("xy0"), arg("xy")));
    def("parallax_correction_inv", &parallax_correction_inv, (
      arg("d"), arg("la"), arg("xy0"), arg("xy")));
  }

}}} // namespace = dxtbx::model::boost_python
