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
#include <cctbx/eltbx/efpfdp.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(fpfdp)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("fpfdp");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<fpfdp> py_fpfdp(this_module, "fpfdp");
    python::export_converters(py_fpfdp);

    py_fpfdp.def(constructor<float, float>());
    py_fpfdp.def(&fpfdp::isValid_fp, "isValid_fp");
    py_fpfdp.def(&fpfdp::isValid_fdp, "isValid_fdp");
    py_fpfdp.def(&fpfdp::isValid, "isValid");
    py_fpfdp.def(&fpfdp::fp, "fp");
    py_fpfdp.def(&fpfdp::fdp, "fdp");
    py_fpfdp.def(&fpfdp::operator(), "__call__");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
