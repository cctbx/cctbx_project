// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/bpl_utils.h>
#include <cctbx/basic/boost_array_bpl.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/uctbx.h>
#include <cctbx/basic/define_range.h>

using namespace uctbx;

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  // Convert a Python tuple to a uc_params object.
  //
  uc_params from_python(PyObject* p, boost::python::type<const uc_params&>)
  {
#   include <cctbx/basic/from_bpl_import.h>
    tuple tup = bpl_utils::tuple_from_python_list_or_tuple(p);
    if (tup.size() == 9) {
      Mx33 G = from_python(p, python::type<const Mx33&>());
      return UnitCell(G).getParameters();
    }
    if (tup.size() > 6) {
      PyErr_SetString(PyExc_ValueError,
        "expecting up to six unit cell parameters"
        " or nine coefficients of the metrical matrix.");
      throw python::error_already_set();
    }
    uc_params ucp = uc_params();
    rangei(tup.size()) {
      ucp[i] = from_python(tup[i].get(), python::type<double>());
    }
    return ucp;
  }

  PyObject* to_python(const uc_params& ucp) {
    return to_python(static_cast<boost::array<double, 6> >(ucp));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

namespace {

# include <cctbx/basic/from_bpl_import.h>

  // Python representation, e.g. output of "print UnitCell()"
  python::string UnitCell_repr(const UnitCell& uc)
  {
    uc_params ucp = uc.getParameters();
#ifdef USE_OSTRINGSTREAM
    std::ostringstream ost;
    ost << ucp.Len(0) << " "
        << ucp.Len(1) << " "
        << ucp.Len(2) << " "
        << ucp.Ang(0) << " "
        << ucp.Ang(1) << " "
        << ucp.Ang(2);
    return python::string(ost.str());
#else
    char buf[128];
    sprintf(buf, "%.6g %.6g %.6g %.6g %.6g %.6g",
      ucp.Len(0), ucp.Len(1), ucp.Len(2),
      ucp.Ang(0), ucp.Ang(1), ucp.Ang(2));
    return python::string(buf);
#endif
  }

  // Support for pickle.
  tuple UnitCell_getinitargs(const UnitCell& uc) {
    tuple initargs(1);
    initargs.set_item(0, ref(to_python(uc.getParameters(false))));
    return initargs;
  }

  // Substitute default parameters.
  uc_params UnitCell_getParameters_0(const UnitCell& uc) {
    return uc.getParameters();
  }
  Mx33 UnitCell_getMetricalMatrix_0(const UnitCell& uc) {
    return uc.getMetricalMatrix();
  }
  UnitCell UnitCell_ChangeBasis_1(const UnitCell& uc, const Mx33& InvCBMxR) {
    return uc.ChangeBasis(InvCBMxR);
  }
  UnitCell UnitCell_ChangeBasis_2(const UnitCell& uc,
                                  const Mx33& InvCBMxR, double RBF) {
    return uc.ChangeBasis(InvCBMxR, RBF);
  }
}

BOOST_PYTHON_MODULE_INIT(uctbx)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("uctbx");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<UnitCell> UnitCell_class(this_module, "UnitCell");
    python::export_converters(UnitCell_class);

    UnitCell_class.def(constructor<const uc_params&>());

    UnitCell_class.def(UnitCell_repr, "__repr__");
    UnitCell_class.def(UnitCell_getinitargs, "__getinitargs__");
    UnitCell_class.def(&UnitCell::getParameters, "getParameters");
    UnitCell_class.def(UnitCell_getParameters_0, "getParameters");
    UnitCell_class.def(&UnitCell::getVolume, "getVolume");
    UnitCell_class.def(&UnitCell::getMetricalMatrix, "getMetricalMatrix");
    UnitCell_class.def(UnitCell_getMetricalMatrix_0, "getMetricalMatrix");
    UnitCell_class.def(&UnitCell::getFractionalizationMatrix,
                                 "getFractionalizationMatrix");
    UnitCell_class.def(&UnitCell::getOrthogonalizationMatrix,
                                 "getOrthogonalizationMatrix");
    UnitCell_class.def(&UnitCell::MaxMillerIndices, "MaxMillerIndices");
    UnitCell_class.def(UnitCell_ChangeBasis_1, "ChangeBasis");
    UnitCell_class.def(UnitCell_ChangeBasis_2, "ChangeBasis");
    UnitCell_class.def(&UnitCell::Q, "Q");
    UnitCell_class.def(&UnitCell::s, "s");
    UnitCell_class.def(&UnitCell::d, "d");
    UnitCell_class.def(&UnitCell::fractionalize, "fractionalize");
    UnitCell_class.def(&UnitCell::orthogonalize, "orthogonalize");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
