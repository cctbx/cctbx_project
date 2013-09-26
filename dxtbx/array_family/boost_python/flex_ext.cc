/*
 * flex_ext.cc
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
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;

  void export_flex_panel();

  BOOST_PYTHON_MODULE(dxtbx_array_family_flex_ext)
  {
    export_flex_panel();
  }

}}} // namespace = dials::af::boost_python
