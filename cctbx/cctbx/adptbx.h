// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 11: Created (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for the handling of anisotropic displacement parameters.
 */

#ifndef CCTBX_ADPTBX_H
#define CCTBX_ADPTBX_H

#include <cctbx/uctbx.h>

//! ADP Toolbox namespace.
namespace adptbx {

  using namespace cctbx;

  const double   TwoPiSquared = 2. * constants::pi * constants::pi;
  const double EightPiSquared = 8. * constants::pi * constants::pi;

  template <class FloatType>
  inline uctbx::Mx33
  Xaniso_as_SymMx33(const boost::array<FloatType, 6>& Xaniso)
  {
    uctbx::Mx33 M;
    M[0] = Xaniso[0];
    M[1] = Xaniso[3];
    M[2] = Xaniso[4];
    M[3] = Xaniso[3];
    M[4] = Xaniso[1];
    M[5] = Xaniso[5];
    M[6] = Xaniso[4];
    M[7] = Xaniso[5];
    M[8] = Xaniso[2];
    return M;
  }

  template <class FloatType>
  struct return_type {};

  template <class FloatType>
  inline boost::array<FloatType, 6>
  SymMx33_as_Xaniso(const uctbx::Mx33& M, return_type<FloatType>)
  {
    boost::array<FloatType, 6> Xaniso;
    Xaniso[0] = M[0];
    Xaniso[1] = M[4];
    Xaniso[2] = M[8];
    Xaniso[3] = M[1];
    Xaniso[4] = M[2];
    Xaniso[5] = M[5];
    return Xaniso;
  }

  template <class FloatType>
  inline boost::array<FloatType, 6>
  A_Xaniso_At(const uctbx::Mx33& A, const boost::array<FloatType, 6>& Xaniso)
  {
    uctbx::Mx33 X = Xaniso_as_SymMx33(Xaniso);
    uctbx::Mx33 AX;
    MatrixLite::multiply<double>(A.elems, X.elems, 3, 3, 3, AX.elems);
    uctbx::Mx33 At;
    MatrixLite::transpose<double>(A.elems, 3, 3, At.elems);
    uctbx::Mx33 AXAt;
    MatrixLite::multiply<double>(AX.elems, At.elems, 3, 3, 3, AXAt.elems);
    return SymMx33_as_Xaniso(AXAt, return_type<FloatType>());
  }

  //! Convert isotropic adp U -> B.
  inline double
  U_as_B(double Uiso) {
    return Uiso * EightPiSquared;
  }
  //! Convert isotropic adp B -> U.
  inline double
  B_as_U(double Biso) {
    return Biso / EightPiSquared;
  }
  //! Convert anisotropic adp U -> B.
  template <class FloatType>
  boost::array<FloatType, 6>
  U_as_B(const boost::array<FloatType, 6>& Uaniso) {
    return EightPiSquared * Uaniso;
  }
  //! Convert anisotropic adp B -> U.
  template <class FloatType>
  boost::array<FloatType, 6>
  B_as_U(const boost::array<FloatType, 6>& Baniso) {
    return (1. / EightPiSquared) * Baniso;
  }

  //! Convert anisotropic adp Uuvrs -> Ustar.
  template <class FloatType>
  boost::array<FloatType, 6>
  Uuvrs_as_Ustar(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Uuvrs) {
    const uctbx::Vec3& R_Len = uc.getLen(true);
    boost::array<FloatType, 6> Ustar;
    Ustar[0] = Uuvrs[0] * (R_Len[0] * R_Len[0]);
    Ustar[1] = Uuvrs[1] * (R_Len[1] * R_Len[1]);
    Ustar[2] = Uuvrs[2] * (R_Len[2] * R_Len[2]);
    Ustar[3] = Uuvrs[3] * (R_Len[0] * R_Len[1]);
    Ustar[4] = Uuvrs[4] * (R_Len[0] * R_Len[2]);
    Ustar[5] = Uuvrs[5] * (R_Len[1] * R_Len[2]);
    return Ustar;
  }
  //! Convert anisotropic adp Ustar -> Uuvrs.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ustar_as_Uuvrs(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ustar) {
    const uctbx::Vec3& R_Len = uc.getLen(true);
    boost::array<FloatType, 6> Uuvrs;
    Uuvrs[0] = Ustar[0] / (R_Len[0] * R_Len[0]);
    Uuvrs[1] = Ustar[1] / (R_Len[1] * R_Len[1]);
    Uuvrs[2] = Ustar[2] / (R_Len[2] * R_Len[2]);
    Uuvrs[3] = Ustar[3] / (R_Len[0] * R_Len[1]);
    Uuvrs[4] = Ustar[4] / (R_Len[0] * R_Len[2]);
    Uuvrs[5] = Ustar[5] / (R_Len[1] * R_Len[2]);
    return Uuvrs;
  }

  //! Convert anisotropic adp Ucart -> Ustar.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ucart_as_Ustar(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ucart) {
    return A_Xaniso_At(uc.getFractionalizationMatrix(), Ucart);
  }
  //! Convert anisotropic adp Ustar -> Ucart.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ustar_as_Ucart(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ustar) {
    return A_Xaniso_At(uc.getOrthogonalizationMatrix(), Ustar);
  }

  //! Convert Uuvrs -> Uiso.
  // From Xtal 3.7.1 source code.
  template <class FloatType>
  inline FloatType
  Uuvrs_as_Uiso(const uctbx::UnitCell& uc,
                const boost::array<FloatType, 6>& Uuvrs)
  {
    const uctbx::Vec3&   Len = uc.getLen(false);
    const uctbx::Vec3& R_Len = uc.getLen(true);
    const uctbx::Vec3& cosAng = uc.get_cosAng(false);
    FloatType Uiso = 0.;
    FloatType LRL[3];
    for(std::size_t i=0;i<3;i++) {
      LRL[i] = Len[i] * R_Len[i];
      Uiso += Uuvrs[i] * LRL[i] * LRL[i];
    }
    Uiso += Uuvrs[3] * 2. * LRL[0] * LRL[1] * cosAng[2];
    Uiso += Uuvrs[4] * 2. * LRL[0] * LRL[2] * cosAng[1];
    Uiso += Uuvrs[5] * 2. * LRL[1] * LRL[2] * cosAng[0];
    return Uiso / 3.;
  }
  //! Convert Uiso -> Uuvrs.
  // From Xtal 3.7.1 source code.
  template <class FloatType>
  inline boost::array<FloatType, 6>
  Uiso_as_Uuvrs(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    const uctbx::Vec3& R_cosAng = uc.get_cosAng(true);
    boost::array<FloatType, 6> Uuvrs;
    Uuvrs.assign(Uiso);
    Uuvrs[3] *= R_cosAng[2];
    Uuvrs[4] *= R_cosAng[1];
    Uuvrs[5] *= R_cosAng[0];
    return Uuvrs;
  }

  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and Biso.
  inline double
  DebyeWallerFactorBiso(double stol2,
                        double Biso)
  {
    return std::exp(-Biso * stol2);
  }
  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and Uiso.
  inline double
  DebyeWallerFactorUiso(double stol2,
                        double Uiso)
  {
    return DebyeWallerFactorBiso(stol2, U_as_B(Uiso));
  }
  //! Isotropic Debye-Waller factor given Miller index and Biso.
  inline double
  DebyeWallerFactorBiso(const uctbx::UnitCell& uc,
                        const Miller::Index& MIx,
                        double Biso)
  {
    return DebyeWallerFactorBiso(uc.Q(MIx) / 4., Biso);
  }
  //! Isotropic Debye-Waller factor given Miller index and Uiso.
  inline double
  DebyeWallerFactorUiso(const uctbx::UnitCell& uc,
                        const Miller::Index& MIx,
                        double Uiso)
  {
    return DebyeWallerFactorBiso(uc, MIx, U_as_B(Uiso));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Ustar.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUstar(const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Ustar)
  {
    return std::exp(-TwoPiSquared * (
        (MIx[0] * MIx[0]) * Ustar[0]
      + (MIx[1] * MIx[1]) * Ustar[1]
      + (MIx[2] * MIx[2]) * Ustar[2]
      + (2 * MIx[0] * MIx[1]) * Ustar[3]
      + (2 * MIx[0] * MIx[2]) * Ustar[4]
      + (2 * MIx[1] * MIx[2]) * Ustar[5]));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Uuvrs.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUuvrs(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Uuvrs)
  {
    return DebyeWallerFactorUstar(MIx, Uuvrs_as_Ustar(uc, Uuvrs));
  }
  //! Anisotropic Debye-Waller factor given Miller index and Ucart.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUcart(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Ucart)
  {
    return DebyeWallerFactorUstar(MIx, Ucart_as_Ustar(uc, Ucart));
  }

} // namespace adptbx

#endif // CCTBX_ADPTBX_H
