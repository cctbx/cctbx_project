// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/caasf.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(caasf_wk1995)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("caasf_wk1995");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<CAASF<5> > py_CAASF_5(this_module, "CAASF_5");
    class_builder<CAASF_WK1995> py_CAASF_WK1995(this_module, "CAASF_WK1995");
    python::export_converters(py_CAASF_WK1995);

    py_CAASF_WK1995.declare_base(py_CAASF_5, python::without_downcast);

    py_CAASF_WK1995.def(constructor<>());
    py_CAASF_WK1995.def(constructor<const std::string&>());
    py_CAASF_WK1995.def(constructor<const std::string&, bool>());
    py_CAASF_WK1995.def(&CAASF_WK1995::Table, "Table");
    py_CAASF_WK1995.def(&CAASF_WK1995::Label, "Label");
    py_CAASF_WK1995.def(&CAASF_WK1995::n_ab, "n_ab");
    py_CAASF_WK1995.def(&CAASF_WK1995::a, "a");
    py_CAASF_WK1995.def(&CAASF_WK1995::b, "b");
    py_CAASF_WK1995.def(&CAASF_WK1995::c, "c");
    py_CAASF_WK1995.def(&CAASF_WK1995::stol2, "stol2");
    py_CAASF_WK1995.def(&CAASF_WK1995::stol, "stol");
    py_CAASF_WK1995.def(&CAASF_WK1995::Q, "Q");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
