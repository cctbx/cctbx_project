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

using namespace cctbx;
using namespace cctbx::adptbx;

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

  boost::array<double, 6>
  py_Ucart_as_Uuvrs(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Ucart) {
    return Ucart_as_Uuvrs(uc, Ucart);
  }
  boost::array<double, 6>
  py_Uuvrs_as_Ucart(const uctbx::UnitCell& uc,
                    const boost::array<double, 6>& Uuvrs) {
    return Uuvrs_as_Ucart(uc, Uuvrs);
  }

  boost::array<double, 6>
  py_Ustar_as_beta(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Ustar) {
    return Ustar_as_beta(Ustar);
  }
  boost::array<double, 6>
  py_beta_as_Ustar(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& beta) {
    return beta_as_Ustar(beta);
  }

  boost::array<double, 6>
  py_Ucart_as_beta(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Ucart) {
    return Ucart_as_beta(uc, Ucart);
  }
  boost::array<double, 6>
  py_beta_as_Ucart(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& beta) {
    return beta_as_Ucart(uc, beta);
  }

  boost::array<double, 6>
  py_Uuvrs_as_beta(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Uuvrs) {
    return Uuvrs_as_beta(uc, Uuvrs);
  }
  boost::array<double, 6>
  py_beta_as_Uuvrs(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& beta) {
    return beta_as_Uuvrs(uc, beta);
  }

  double
  py_Ucart_as_Uiso(const boost::array<double, 6>& Ucart) {
    return Ucart_as_Uiso(Ucart);
  }
  boost::array<double, 6>
  py_Uiso_as_Ucart(const double& Uiso) {
    return Uiso_as_Ucart(Uiso);
  }

  double
  py_Ustar_as_Uiso(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Ustar) {
    return Ustar_as_Uiso(uc, Ustar);
  }
  boost::array<double, 6>
  py_Uiso_as_Ustar(const uctbx::UnitCell& uc,
                   const double& Uiso) {
    return Uiso_as_Ustar(uc, Uiso);
  }

  double
  py_Uuvrs_as_Uiso(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& Uuvrs) {
    return Uuvrs_as_Uiso(uc, Uuvrs);
  }
  boost::array<double, 6>
  py_Uiso_as_Uuvrs(const uctbx::UnitCell& uc,
                   const double& Uiso) {
    return Uiso_as_Uuvrs(uc, Uiso);
  }

  double
  py_beta_as_Uiso(const uctbx::UnitCell& uc,
                   const boost::array<double, 6>& beta) {
    return beta_as_Uiso(uc, beta);
  }
  boost::array<double, 6>
  py_Uiso_as_beta(const uctbx::UnitCell& uc,
                   const double& Uiso) {
    return Uiso_as_beta(uc, Uiso);
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
  py_DebyeWallerFactor_beta(const Miller::Index& MIx,
                            const boost::array<double, 6>& beta) {
    return DebyeWallerFactor_beta(MIx, beta);
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
  py_Eigenvalues(const boost::array<double, 6>& adp) {
    return Eigenvalues(adp);
  }

  bool
  py_isPositiveDefinite_adp_eigenvalues(
    const boost::array<double, 3>& adp_eigenvalues) {
    return isPositiveDefinite(adp_eigenvalues);
  }

  bool
  py_isPositiveDefinite_adp(
    const boost::array<double, 6>& adp) {
    return isPositiveDefinite(adp);
  }

  void
  py_CheckPositiveDefinite_adp_eigenvalues(
    const boost::array<double, 3>& adp_eigenvalues) {
    CheckPositiveDefinite(adp_eigenvalues);
  }

  void
  py_CheckPositiveDefinite_adp(
    const boost::array<double, 6>& adp) {
    CheckPositiveDefinite(adp);
  }

  // We need this wrapper only to make Visual C++ 6 happy.
  boost::array<double, 3>
  py_Eigensystem_vectors(const Eigensystem<double>& ES, std::size_t i) {
    return ES.vectors(i);
  }

}

BOOST_PYTHON_MODULE_INIT(adptbx)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("adptbx");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  python::import_converters<uctbx::UnitCell>
  UnitCell_converters("cctbx.uctbx", "UnitCell");

  class_builder<Eigensystem<double> >
  py_Eigensystem(this_module, "Eigensystem");

  this_module.def(py_U_as_B_iso, "U_as_B");
  this_module.def(py_B_as_U_iso, "B_as_U");
  this_module.def(py_U_as_B_ansio, "U_as_B");
  this_module.def(py_B_as_U_ansio, "B_as_U");

  this_module.def(py_Uuvrs_as_Ustar, "Uuvrs_as_Ustar");
  this_module.def(py_Ustar_as_Uuvrs, "Ustar_as_Uuvrs");
  this_module.def(py_Ucart_as_Ustar, "Ucart_as_Ustar");
  this_module.def(py_Ustar_as_Ucart, "Ustar_as_Ucart");
  this_module.def(py_Ucart_as_Uuvrs, "Ucart_as_Uuvrs");
  this_module.def(py_Uuvrs_as_Ucart, "Uuvrs_as_Ucart");
  this_module.def(py_Ustar_as_beta, "Ustar_as_beta");
  this_module.def(py_beta_as_Ustar, "beta_as_Ustar");
  this_module.def(py_Ucart_as_beta, "Ucart_as_beta");
  this_module.def(py_beta_as_Ucart, "beta_as_Ucart");
  this_module.def(py_Uuvrs_as_beta, "Uuvrs_as_beta");
  this_module.def(py_beta_as_Uuvrs, "beta_as_Uuvrs");

  this_module.def(py_Ucart_as_Uiso, "Ucart_as_Uiso");
  this_module.def(py_Uiso_as_Ucart, "Uiso_as_Ucart");
  this_module.def(py_Ustar_as_Uiso, "Ustar_as_Uiso");
  this_module.def(py_Uiso_as_Ustar, "Uiso_as_Ustar");
  this_module.def(py_Uuvrs_as_Uiso, "Uuvrs_as_Uiso");
  this_module.def(py_Uiso_as_Uuvrs, "Uiso_as_Uuvrs");
  this_module.def(py_beta_as_Uiso, "beta_as_Uiso");
  this_module.def(py_Uiso_as_beta, "Uiso_as_beta");

  this_module.def(py_DebyeWallerFactorBiso_2, "DebyeWallerFactorBiso");
  this_module.def(py_DebyeWallerFactorUiso_2, "DebyeWallerFactorUiso");
  this_module.def(py_DebyeWallerFactorBiso_3, "DebyeWallerFactorBiso");
  this_module.def(py_DebyeWallerFactorUiso_3, "DebyeWallerFactorUiso");
  this_module.def(py_DebyeWallerFactorUstar, "DebyeWallerFactorUstar");
  this_module.def(py_DebyeWallerFactor_beta, "DebyeWallerFactor_beta");
  this_module.def(py_DebyeWallerFactorUuvrs, "DebyeWallerFactorUuvrs");
  this_module.def(py_DebyeWallerFactorUcart, "DebyeWallerFactorUcart");

  this_module.def(py_Eigenvalues, "Eigenvalues");
  this_module.def(py_isPositiveDefinite_adp_eigenvalues,
                    "isPositiveDefinite");
  this_module.def(py_isPositiveDefinite_adp,
                    "isPositiveDefinite");
  this_module.def(py_CheckPositiveDefinite_adp_eigenvalues,
                    "CheckPositiveDefinite");
  this_module.def(py_CheckPositiveDefinite_adp,
                    "CheckPositiveDefinite");

  py_Eigensystem.def(constructor<>());
  py_Eigensystem.def(constructor<const boost::array<double, 6>&>());
  py_Eigensystem.def(constructor<const boost::array<double, 6>&, double>());
  py_Eigensystem.def(py_Eigensystem_vectors, "vectors");
  py_Eigensystem.def(&Eigensystem<double>::values, "values");
}
