// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 11: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/bpl_utils.h>
#include <cctbx/basic/boost_array_bpl.h>
#include <cctbx/miller_bpl.h>
#include <cctbx/adptbx.h>

using namespace adptbx;

namespace {

# include <cctbx/basic/from_bpl_import.h>

  double
  py_U_as_B_iso(double Uiso) {
    return U_as_B(Uiso);
  }
  double
  py_B_as_U_iso(double Biso) {
    return B_as_U(Biso);
  }
  boost::array<double, 6>
  py_U_as_B_ansio(const boost::array<double, 6>& Uaniso) {
    return U_as_B(Uaniso);
  }
  boost::array<double, 6>
  py_B_as_U_ansio(const boost::array<double, 6>& Baniso) {
    return B_as_U(Baniso);
  }

  boost::array<double, 6>
  py_Uuvrs_as_Ustar(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Uuvrs) {
    return Uuvrs_as_Ustar(uc, Uuvrs);
  }
  boost::array<double, 6>
  py_Ustar_as_Uuvrs(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Ustar) {
    return Ustar_as_Uuvrs(uc, Ustar);
  }

  boost::array<double, 6>
  py_Ucart_as_Ustar(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Ucart) {
    return Ucart_as_Ustar(uc, Ucart);
  }
  boost::array<double, 6>
  py_Ustar_as_Ucart(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Ustar) {
    return Ustar_as_Ucart(uc, Ustar);
  }

  double
  py_Uuvrs_as_Uiso(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Uuvrs)
  {
    return Uuvrs_as_Uiso(uc, Uuvrs);
  }
  boost::array<double, 6>
  py_Uiso_as_Uuvrs(const uctbx::UnitCell& uc,
                   const double& Uiso)
  {
    return Uiso_as_Uuvrs(uc, Uiso);
  }

  double
  py_DebyeWallerFactorBiso_2(double stol2, double Biso) {
    return DebyeWallerFactorBiso(stol2, Biso);
  }
  double
  py_DebyeWallerFactorUiso_2(double stol2, double Uiso) {
    return DebyeWallerFactorUiso(stol2, Uiso);
  }
  double
  py_DebyeWallerFactorBiso_3(const uctbx::UnitCell& uc,
                             const Miller::Index& MIx,
                             double Biso) {
    return DebyeWallerFactorBiso(uc, MIx, Biso);
  }
  double
  py_DebyeWallerFactorUiso_3(const uctbx::UnitCell& uc,
                             const Miller::Index& MIx,
                             double Uiso) {
    return DebyeWallerFactorUiso(uc, MIx, Uiso);
  }

  double
  py_DebyeWallerFactorUstar(const Miller::Index& MIx,
                            const boost::array<double, 6>& Ustar) {
    return DebyeWallerFactorUstar(MIx, Ustar);
  }
  double
  py_DebyeWallerFactorUuvrs(const uctbx::UnitCell& uc,
                            const Miller::Index& MIx,
                            const boost::array<double, 6>& Uuvrs) {
    return DebyeWallerFactorUuvrs(uc, MIx, Uuvrs);
  }
  double
  py_DebyeWallerFactorUcart(const uctbx::UnitCell& uc,
                            const Miller::Index& MIx,
                            const boost::array<double, 6>& Ucart) {
    return DebyeWallerFactorUcart(uc, MIx, Ucart);
  }

  boost::array<double, 3>
  py_EigenValues(const boost::array<double, 6>& adp) {
    return EigenValues(adp);
  }

  tuple
  py_EigenVectors(const boost::array<double, 6>& adp) {
    boost::array<boost::array<double, 3>, 3>
    EVs = EigenVectors(adp);
    tuple result(3);
    for(std::size_t i=0;i<result.size();i++) {
      result.set_item(i, boost::python::ref(to_python(EVs[i])));
    }
    return result;
  }

}

BOOST_PYTHON_MODULE_INIT(adptbx)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("adptbx");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<uctbx::UnitCell>
    UnitCell_converters("uctbx", "UnitCell");

    this_module.def(py_U_as_B_iso, "U_as_B");
    this_module.def(py_B_as_U_iso, "B_as_U");
    this_module.def(py_U_as_B_ansio, "U_as_B");
    this_module.def(py_B_as_U_ansio, "B_as_U");

    this_module.def(py_Uuvrs_as_Ustar, "Uuvrs_as_Ustar");
    this_module.def(py_Ustar_as_Uuvrs, "Ustar_as_Uuvrs");
    this_module.def(py_Ucart_as_Ustar, "Ucart_as_Ustar");
    this_module.def(py_Ustar_as_Ucart, "Ustar_as_Ucart");

    this_module.def(py_Uuvrs_as_Uiso, "Uuvrs_as_Uiso");
    this_module.def(py_Uiso_as_Uuvrs, "Uiso_as_Uuvrs");

    this_module.def(py_DebyeWallerFactorBiso_2, "DebyeWallerFactorBiso");
    this_module.def(py_DebyeWallerFactorUiso_2, "DebyeWallerFactorUiso");
    this_module.def(py_DebyeWallerFactorBiso_3, "DebyeWallerFactorBiso");
    this_module.def(py_DebyeWallerFactorUiso_3, "DebyeWallerFactorUiso");
    this_module.def(py_DebyeWallerFactorUstar, "DebyeWallerFactorUstar");
    this_module.def(py_DebyeWallerFactorUuvrs, "DebyeWallerFactorUuvrs");
    this_module.def(py_DebyeWallerFactorUcart, "DebyeWallerFactorUcart");

    this_module.def(py_EigenValues, "EigenValues");
    this_module.def(py_EigenVectors, "EigenVectors");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
