/*
 * command_line_ext.cc
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

namespace dxtbx { namespace command_line { namespace boost_python {

  using namespace boost::python;

  void export_to_ewald_sphere_helpers();

  BOOST_PYTHON_MODULE(dxtbx_command_line_ext)
  {
    export_to_ewald_sphere_helpers();
  }

}}} // namespace dxtbx::command_line::boost_python
