/*
 * flex_panel.cc
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
#include <dxtbx/model/panel.h>

namespace dxtbx { namespace af { namespace boost_python {

  using namespace boost::python;
  using dxtbx::model::Panel;
  
  void export_flex_panel()
  {
    scitbx::af::boost_python::flex_wrapper <Panel, 
      return_internal_reference<> >::plain("panel");
  }

}}} // namespace dials::af::boost_python
