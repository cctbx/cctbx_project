// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 22: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/bpl_utils.h>
#include <cctbx/basic/boost_array_bpl.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/coordinates_bpl.h>
#include <cctbx/xray_scatterer.h>

namespace {


}

BOOST_PYTHON_MODULE_INIT(sftbx)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("sftbx");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    UnitCell_converters("uctbx", "UnitCell");
    python::import_converters<sgtbx::SpaceGroup>
    SpaceGroup_converters("sgtbx", "SpaceGroup");
    python::import_converters<eltbx::CAASF_WK1995>
    CAASF_WK1995_converters("eltbx.caasf_wk1995", "CAASF_WK1995");

    class_builder<cctbx::XrayScatterer<double, 5> >
    py_XrayScatterer(this_module, "XrayScatterer");

    py_XrayScatterer.def(constructor<>());
    py_XrayScatterer.def(constructor<
      const std::string&,
      const eltbx::CAASF_WK1995&,
      const cctbx::fractional<double>&,
      const double&,
      const double&,
      const std::complex<double>&>());
    py_XrayScatterer.def(
      &cctbx::XrayScatterer<double, 5>::DetermineMultiplicity,
                                       "DetermineMultiplicity");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
