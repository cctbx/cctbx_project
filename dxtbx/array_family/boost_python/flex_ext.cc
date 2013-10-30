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
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/tiny_types.h>

namespace dxtbx { namespace af { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::int4;
  
  void export_flex_panel();

  BOOST_PYTHON_MODULE(dxtbx_array_family_flex_ext)
  {
    scitbx::af::boost_python::flex_wrapper <int4>::plain("flex_int4");  
  
    export_flex_panel();
  }

}}} // namespace = dials::af::boost_python
