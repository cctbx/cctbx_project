/*
 * experiment_ext.cc
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

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_beam();
  void export_goniometer();
  void export_detector();
//  void export_scan();

  BOOST_PYTHON_MODULE(dxtbx_model_ext)
  {
    export_beam();
    export_goniometer();
    export_detector();
//    export_scan();
  }

}}} // namespace dxtbx::model::boost_python
