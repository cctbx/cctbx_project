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
#include <cctbx/eltbx/icsd_radii.h>

using namespace cctbx;
using namespace cctbx::eltbx;

BOOST_PYTHON_MODULE_INIT(icsd_radii)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("icsd_radii");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  class_builder<ICSD_Radius> py_ICSD_Radius(this_module, "ICSD_Radius");

  py_ICSD_Radius.def(constructor<>());
  py_ICSD_Radius.def(constructor<const std::string&>());
  py_ICSD_Radius.def(constructor<const std::string&, bool>());
  py_ICSD_Radius.def(&ICSD_Radius::Label, "Label");
  py_ICSD_Radius.def(&ICSD_Radius::Radius, "Radius");
}
