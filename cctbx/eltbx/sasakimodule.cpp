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
#include <cctbx/eltbx/sasaki.h>

using namespace cctbx;
using namespace cctbx::eltbx;

BOOST_PYTHON_MODULE_INIT(sasaki)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("sasaki");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  python::import_converters<fpfdp>
  py_fpfdp("cctbx_boost.eltbx.fpfdp", "fpfdp");
  class_builder<Sasaki> py_Sasaki(this_module, "Sasaki");

  py_Sasaki.def(constructor<>());
  py_Sasaki.def(constructor<const std::string&>());
  py_Sasaki.def(constructor<const std::string&, bool>());
  py_Sasaki.def(&Sasaki::Label, "Label");
  py_Sasaki.def(&Sasaki::Z, "Z");
  py_Sasaki.def(&Sasaki::operator(), "__call__");
}
