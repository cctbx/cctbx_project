// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/python/class_builder.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/tinypse.h>

using namespace cctbx;
using namespace cctbx::eltbx;

BOOST_PYTHON_MODULE_INIT(tinypse)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("tinypse");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  class_builder<TinyPSE> py_TinyPSE(this_module, "TinyPSE");

  py_TinyPSE.def(constructor<>());
  py_TinyPSE.def(constructor<const std::string&>());
  py_TinyPSE.def(constructor<const std::string&, bool>());
  py_TinyPSE.def(constructor<int>());
  py_TinyPSE.def(&TinyPSE::Z, "Z");
  py_TinyPSE.def(&TinyPSE::Symbol, "Symbol");
  py_TinyPSE.def(&TinyPSE::Name, "Name");
  py_TinyPSE.def(&TinyPSE::Weight, "Weight");
}
