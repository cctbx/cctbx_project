// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created, based on fftbxmodule.cpp (rwgk)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/miller_bpl.h>
#include <cctbx/std_vector_bpl.h>

#include <cctbx/sgtbx/matrix.h>
#include <cctbx/xray_scatterer.h>

namespace {

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    using namespace cctbx;

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<sgtbx::RTMx>
    py_RTMx("cctbx.sgtbx", "RTMx");

    python::import_converters<
      sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >
    py_XrayScatterer("cctbx.sftbx", "XrayScatterer");

    python::wrap_std_vector<int>::run(this_module, "int");
    python::wrap_std_vector<double>::run(this_module, "double");
    python::wrap_std_vector<std::complex<double> >::run(this_module,
      "complex_double");

    python::wrap_std_vector<Miller::Index>::run(this_module,
      "Miller_Index");

    python::wrap_std_vector<sgtbx::RTMx>::run(this_module, "RTMx");

    python::wrap_std_vector<
      sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >::run(this_module,
      "XrayScatterer");
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE_INIT(std_vector)
{
  boost::python::module_builder this_module("std_vector");
  init_module(this_module);
}
