// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/bpl_utils.h>
#include <cctbx/coordinates_bpl.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/uctbx.h>
#include <cctbx/basic/define_range.h>

using namespace cctbx;
using namespace cctbx::uctbx;

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  // Convert a Python tuple to a uc_params object.
  //
  uc_params from_python(PyObject* p, boost::python::type<const uc_params&>)
  {
#   include <cctbx/basic/from_bpl_import.h>
    tuple tup = bpl_utils::tuple_from_python_list_or_tuple(p);
    if (tup.size() == 9) {
      af::double9 G = from_python(p, python::type<const af::double9&>());
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
    return to_python(static_cast<cctbx::af::double6 >(ucp));
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
  af::double9 UnitCell_getMetricalMatrix_0(const UnitCell& uc) {
    return uc.getMetricalMatrix();
  }
  double
  UnitCell_Length2(const UnitCell& uc,
                   const fractional<double>& Xf) {
    return uc.Length2(Xf);
  }
  double
  UnitCell_Length(const UnitCell& uc,
                  const fractional<double>& Xf) {
    return uc.Length(Xf);
  }
  double
  UnitCell_Distance2(const UnitCell& uc,
                     const fractional<double>& Xf,
                     const fractional<double>& Yf) {
    return uc.Distance2(Xf, Yf);
  }
  double
  UnitCell_Distance(const UnitCell& uc,
                    const fractional<double>& Xf,
                    const fractional<double>& Yf) {
    return uc.Distance(Xf, Yf);
  }
  double
  UnitCell_modShortLength2(const UnitCell& uc,
                           const fractional<double>& Xf) {
    return uc.modShortLength2(Xf);
  }
  double
  UnitCell_modShortLength(const UnitCell& uc,
                          const fractional<double>& Xf) {
    return uc.modShortLength(Xf);
  }
  double
  UnitCell_modShortDistance2(const UnitCell& uc,
                             const fractional<double>& Xf,
                             const fractional<double>& Yf) {
    return uc.modShortDistance2(Xf, Yf);
  }
  double
  UnitCell_modShortDistance(const UnitCell& uc,
                            const fractional<double>& Xf,
                            const fractional<double>& Yf) {
    return uc.modShortDistance(Xf, Yf);
  }
  UnitCell
  UnitCell_ChangeBasis_1(const UnitCell& uc, const af::double9& InvCBMxR) {
    return uc.ChangeBasis(InvCBMxR);
  }
  fractional<double>
  UnitCell_fractionalize(const UnitCell& uc,
                         const cartesian<double>& Xc) {
    return uc.fractionalize(Xc);
  }
  cartesian<double>
  UnitCell_orthogonalize(const UnitCell& uc,
                         const fractional<double>& Xf) {
    return uc.orthogonalize(Xf);
  }

  bool UnitCell_isEqual_1(const UnitCell& uc, const UnitCell& other) {
    return uc.isEqual(other);
  }

  double
  UnitCell_two_theta_s_2(
    const UnitCell& uc,
    const Miller::Index& MIx,
    double wavelength)
  {
    return uc.two_theta(MIx, wavelength);
  }

  af::shared<double>
  UnitCell_two_theta_a_2(
    const UnitCell& uc,
    af::shared<Miller::Index> MIx,
    double wavelength)
  {
    return uc.two_theta(MIx, wavelength);
  }

  af::double3
  modShortDifference(const af::double3& Xf, const af::double3& Yf) {
    return fractional<double>(Xf - Yf).modShort();
  }

}

BOOST_PYTHON_MODULE_INIT(uctbx)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("uctbx");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  python::import_converters<af::shared<double> >
  py_shared_double("cctbx_boost.arraytbx.shared", "double");

  python::import_converters<af::shared<Miller::Index> >
  py_shared_Miller_Index("cctbx_boost.arraytbx.shared", "Miller_Index");

  class_builder<UnitCell> UnitCell_class(this_module, "UnitCell");
  python::export_converters(UnitCell_class);

  UnitCell_class.def(constructor<>());
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
  UnitCell_class.def(UnitCell_Length2, "Length2");
  UnitCell_class.def(UnitCell_Length, "Length");
  UnitCell_class.def(UnitCell_Distance2, "Distance2");
  UnitCell_class.def(UnitCell_Distance, "Distance");
  UnitCell_class.def(UnitCell_modShortLength2, "modShortLength2");
  UnitCell_class.def(UnitCell_modShortLength, "modShortLength");
  UnitCell_class.def(UnitCell_modShortDistance2, "modShortDistance2");
  UnitCell_class.def(UnitCell_modShortDistance, "modShortDistance");
  UnitCell_class.def(&UnitCell::MaxMillerIndices, "MaxMillerIndices");
  UnitCell_class.def(UnitCell_ChangeBasis_1, "ChangeBasis");
  UnitCell_class.def(
    (UnitCell (UnitCell::*)(const af::double9&, double) const)
    &UnitCell::ChangeBasis, "ChangeBasis");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::Q, "Q");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&) const)
    &UnitCell::Q, "Q");
  UnitCell_class.def(
    (double (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::max_Q, "max_Q");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::stol2, "stol2");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&) const)
    &UnitCell::stol2, "stol2");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::two_stol, "two_stol");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&) const)
    &UnitCell::two_stol, "two_stol");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::stol, "stol");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&) const)
    &UnitCell::stol, "stol");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(const af::shared<Miller::Index>&) const)
    &UnitCell::d, "d");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&) const)
    &UnitCell::d, "d");
  UnitCell_class.def(
    (af::shared<double> (UnitCell::*)(
      af::shared<Miller::Index>, double, bool) const)
    &UnitCell::two_theta, "two_theta");
  UnitCell_class.def(UnitCell_two_theta_s_2, "two_theta");
  UnitCell_class.def(
    (double (UnitCell::*)(const Miller::Index&, double, bool) const)
    &UnitCell::two_theta, "two_theta");
  UnitCell_class.def(UnitCell_two_theta_a_2, "two_theta");
  UnitCell_class.def(UnitCell_fractionalize, "fractionalize");
  UnitCell_class.def(UnitCell_orthogonalize, "orthogonalize");
  UnitCell_class.def(&UnitCell::getLongestVector2, "getLongestVector2");
  UnitCell_class.def(UnitCell_isEqual_1, "isEqual");
  UnitCell_class.def(&UnitCell::isEqual, "isEqual");

  this_module.def(modShortDifference, "modShortDifference");
}
