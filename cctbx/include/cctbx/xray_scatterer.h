// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 22: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_XRAY_SCATTERER_H
#define CCTBX_XRAY_SCATTERER_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/coordinates.h>
#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/eltbx/caasf.h>
#include <cctbx/adptbx.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/tiny_algebra.h>
#include <cctbx/array_family/reductions.h>

namespace cctbx {
  //! Structure Factor Toolbox namespace.
  namespace sftbx {

  //! Structure factor without Debye-Waller factor.
  /*! Sum of exp(2 pi j H S X) over all symmetry operations S.
      j is the imaginary number.
   */
  template <class FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  const Miller::Index& H,
                  const fractional<FloatType>& X)
  {
    using constants::pi;
    std::complex<FloatType> F(FloatType(0));
    for (int s = 0; s < SgOps.nSMx(); s++) {
      Miller::Index HR = H * SgOps[s].Rpart();
      FloatType HRX = HR * X;
      sgtbx::TrVec T = SgOps[s].Tpart();
      for (int i = 0; i < SgOps.fInv(); i++) {
        if (i) {
          HRX = -HRX;
          T = SgOps.InvT() - T;
        }
        for (int l = 0; l < SgOps.nLTr(); l++) {
          FloatType HT = FloatType(H * (T + SgOps.LTr(l))) / SgOps.TBF();
          FloatType phase = 2. * pi * (HRX + HT);
          F += std::complex<FloatType>(std::cos(phase), std::sin(phase));
        }
      }
    }
    return F;
  }

  //! Structure factor with isotropic Debye-Waller factor given Uiso.
  /*! Sum of exp(2 pi j H S X) over all symmetry operations S,
      multiplied by cctbx::adptbx::DebyeWallerFactorUiso().
      j is the imaginary number.
      <p>
      This function is provided for symmetry with the
      structure factor calculation given anisotropic
      displacement parameters.
   */
  template <class FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  const uctbx::UnitCell& uc,
                  const Miller::Index& H,
                  const fractional<FloatType>& X,
                  double Uiso)
  {
    return
     StructureFactor(SgOps, H, X) * adptbx::DebyeWallerFactorUiso(uc, H, Uiso);
  }

  //! Structure factor with isotropic Debye-Waller factor given Uiso.
  /*! Sum of exp(2 pi j H S X) over all symmetry operations S,
      multiplied by cctbx::adptbx::DebyeWallerFactorUiso().
      j is the imaginary number.
      <p>
      This function is provided for symmetry with the
      structure factor calculation given anisotropic
      displacement parameters.
   */
  template <class FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  double stol2,
                  const Miller::Index& H,
                  const fractional<FloatType>& X,
                  double Uiso)
  {
    return
     StructureFactor(SgOps, H, X) * adptbx::DebyeWallerFactorUiso(stol2, Uiso);
  }

  //! Structure factor with anisotropic Debye-Waller factor given Ustar.
  /*! Sum of cctbx::adptbx::DebyeWallerFactorUstar() * exp(2 pi j H S X)
      over all symmetry operations S.
      j is the imaginary number.
   */
  template <class FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  const Miller::Index& H,
                  const fractional<FloatType>& X,
                  const af::tiny<FloatType, 6>& Ustar)
  {
    using constants::pi;
    std::complex<FloatType> F(FloatType(0));
    for (int s = 0; s < SgOps.nSMx(); s++) {
      Miller::Index HR = H * SgOps[s].Rpart();
      FloatType HRX = HR * X;
      sgtbx::TrVec T = SgOps[s].Tpart();
      std::complex<FloatType> Fs(0., 0.);
      for (int i = 0; i < SgOps.fInv(); i++) {
        if (i) {
          HRX = -HRX;
          T = SgOps.InvT() - T;
        }
        for (int l = 0; l < SgOps.nLTr(); l++) {
          FloatType HT = FloatType(H * (T + SgOps.LTr(l))) / SgOps.TBF();
          FloatType phase = 2. * pi * (HRX + HT);
          Fs += std::complex<FloatType>(std::cos(phase), std::sin(phase));
        }
      }
      F += Fs * adptbx::DebyeWallerFactorUstar(HR, Ustar);
    }
    return F;
  }

  //! XXX
  template <typename FloatType>
  af::tiny<FloatType, 3>
  StructureFactor_dX(
    const sgtbx::SpaceGroup& SgOps,
    const Miller::Index& H,
    const fractional<FloatType>& X,
    const std::complex<FloatType>& phase_indep_coeff,
    const std::complex<FloatType>& dTarget_dFcalc)
  {
    using constants::pi;
    af::tiny<FloatType, 3> result;
    result.fill(FloatType(0));
    for (int s = 0; s < SgOps.nSMx(); s++) {
      Miller::Index HR = H * SgOps[s].Rpart();
      FloatType HRX = HR * X;
      sgtbx::TrVec T = SgOps[s].Tpart();
      for (int i = 0; i < SgOps.fInv(); i++) {
        if (i) {
          HR = -HR;
          HRX = -HRX;
          T = SgOps.InvT() - T;
        }
        for (int l = 0; l < SgOps.nLTr(); l++) {
          FloatType HT = FloatType(H * (T + SgOps.LTr(l))) / SgOps.TBF();
          FloatType phase = 2. * pi * (HRX + HT);
          std::complex<FloatType> fX = std::complex<FloatType>(
            std::cos(phase), std::sin(phase));
          fX *= phase_indep_coeff;
          FloatType c =   fX.real() * dTarget_dFcalc.imag()
                        - fX.imag() * dTarget_dFcalc.real();
          for(int j=0;j<3;j++) {
            result[j] += HR[j] * c;
          }
        }
      }
    }
    result *= 2. * pi;
    return result;
  }

  /*! \brief This class groups the information about an atom that
      is needed for a structure factor calculation.
   */
  /*! Constructors are provided for scatterers with isotropic
      and anisotropic displacement parameters (temperature factors).
      <p>
      The anisotropic displacement parameters have to be
      "Ustar." Converters between different conventions for
      the representation of anisotropic displacement
      parameters are provided by the cctbx::adptbx.
      <p>
      The ApplySymmetry() function has to be called before
      the scatterer can be used in a structure factor
      calculation.
   */
  template <typename FloatType, typename CAASF_Type>
  class XrayScatterer
  {
    public:
      typedef FloatType float_type;
      //! Default constructor. Data members are not initialized!
      XrayScatterer() {}
      //! Constructor with isotropic Debye-Waller factor.
      XrayScatterer(const std::string& Label,
                    const CAASF_Type& CAASF,
                    const std::complex<FloatType>& fpfdp,
                    const fractional<FloatType>& Coordinates,
                    const FloatType& Occ,
                    const FloatType& Uiso)
        : m_Label(Label),
          m_CAASF(CAASF),
          m_fpfdp(fpfdp),
          m_Coordinates(Coordinates),
          m_Occ(Occ),
          m_Anisotropic(false),
          m_M(0),
          m_w(0)
      {
        m_U.fill(0.);
        m_U[0] = Uiso;
      }
      //! Constructor with anisotropic Debye-Waller factor.
      XrayScatterer(const std::string& Label,
                    const CAASF_Type& CAASF,
                    const std::complex<FloatType>& fpfdp,
                    const fractional<FloatType>& Coordinates,
                    const FloatType& Occ,
                    const af::tiny<FloatType, 6>& Uaniso)
        : m_Label(Label),
          m_CAASF(CAASF),
          m_fpfdp(fpfdp),
          m_Coordinates(Coordinates),
          m_Occ(Occ),
          m_Anisotropic(true),
          m_U(Uaniso),
          m_M(0),
          m_w(0)
      {}
      //! Access the Label.
      const std::string& Label() const { return m_Label; }
      //! Access the analytical approximation to the scattering factor.
      /*! See eltbx::CAASF_IT1992 and eltbx::CAASF_WK1995.
       */
      const CAASF_Type& CAASF() const { return m_CAASF; }
      //! Access f-prime and f-double-prime.
      /*! f-prime is the dispersive contribution to the scattering
          factor. f-double-prime is the anomalous contribution.
       */
      const std::complex<FloatType>& fpfdp() const { return m_fpfdp; }
      //! Access the fractional coordinates.
      const fractional<FloatType>& Coordinates() const {
        return m_Coordinates;
      }
      //! Access the occupancy factor.
      const FloatType& Occ() const { return m_Occ; }
      //! Test if the scatterer is anisotropic.
      bool isAnisotropic() const { return m_Anisotropic; }
      //! Access to the isotropic displacement parameter.
      /*! It is the responsibility of the caller to ensure
          that the scatterer is isotropic.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by the cctbx::adptbx.
       */
      const FloatType& Uiso() const { return m_U[0]; }
      //! Access to the anisotropic displacement parameters.
      /*! It is the responsibility of the caller to ensure
          that the scatterer is anisotropic.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by the cctbx::adptbx.
       */
      const af::tiny<FloatType, 6>& Uaniso() { return m_U; }
      //! Redefine the coordinates.
      void set_Coordinates(const fractional<FloatType>& Coordinates) {
        m_Coordinates = Coordinates;
      }
      //! Compute multiplicity and average anisotropic displacement parameters.
      /*! The given unit cell and space group are used to determine
          the site symmetry of the scatterer see
          cctbx::sgtbx::SpecialPositionSnapParameters and
          cctbx::sgtbx::SiteSymmetry for information on the
          MinMateDistance parameter.
          <p>
          For scatterers with anisotropic displacement
          parameters, the average Ustar is determined using
          cctbx::sgtbx::SiteSymmetry::AverageUstar. If
          Ustar_tolerance is greater than zero, an
          exception is thrown if the discrepancy between
          components of Ustar before and after the
          application of the site symmetry is greater than
          the given tolerance.
          <p>
          If TestPositiveDefiniteness = true,
          for scatterers with anisotropic displacement
          parameters it is tested if the symmetry-averaged
          Ustar tensor is positive definite. An exception
          is thrown if this is not the case.
          <p>
          See also:
            cctbx::sgtbx::SpecialPositionSnapParameters,
            cctbx::sgtbx::SiteSymmetry,
            cctbx::adptbx::CheckPositiveDefinite
       */
      void ApplySymmetry(const uctbx::UnitCell& UC,
                         const sgtbx::SpaceGroup& SgOps,
                         double MinMateDistance = 0.5,
                         double Ustar_tolerance = 0.1,
                         bool TestPositiveDefiniteness = true)
      {
        sgtbx::SpecialPositionSnapParameters
        SnapParameters(UC, SgOps, true, MinMateDistance);
        sgtbx::SiteSymmetry
        SS(SnapParameters, m_Coordinates);
        m_Coordinates = SS.SnapPosition();
        m_M = SS.M();
        m_w = m_Occ * m_M / SgOps.OrderZ();
        if (m_Anisotropic) {
          if (Ustar_tolerance > 0.) {
            SS.CheckUstar(m_U, Ustar_tolerance);
          }
          m_U = SS.AverageUstar(m_U);
          if (TestPositiveDefiniteness) {
            adptbx::CheckPositiveDefinite(m_U);
          }
        }
      }
      //! Access the multiplicity computed by ApplySymmetry().
      int M() const { return m_M; }
      //! Access the "weight" computed by ApplySymmetry().
      /*! The weight is defined as Occ() * (M() / SgOps.OrderZ()),
          with the SgOps passed to ApplySymmetry(). This weight
          is used in the structure factor calculation.
       */
      const FloatType& w() const { return m_w; }
      //! XXX
      void CheckMultiplicity() const {
        if (m_M == 0) {
          throw error(
            "ApplySymmetry() has not been called for this scatterer.");
        }
      }
      //! XXX
      fractional<FloatType>
      difference(const fractional<FloatType>& other_coordinates) const {
        return m_Coordinates - other_coordinates;
      }
      //! XXX
      FloatType
      distance2(const uctbx::UnitCell& ucell,
                const fractional<FloatType>& other_coordinates) const {
        return ucell.orthogonalize(difference(other_coordinates)).Length2();
      }
      //! XXX
      FloatType
      distance(const uctbx::UnitCell& ucell,
               const fractional<FloatType>& other_coordinates) const {
        return std::sqrt(distance2(ucell, other_coordinates));
      }
      /*! \brief Contribution of the (one) scatterer to the
          structure factor with the Miller index H.
       */
      /*! Q is a d-spacing measure (cctbx::uctbx::UnitCell::Q()).
       */
      std::complex<FloatType>
      StructureFactor(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q) const
      {
        if (!m_Anisotropic) {
          return
              m_w * (m_CAASF.Q(Q) + m_fpfdp)
            * sftbx::StructureFactor(SgOps, Q / 4., H, m_Coordinates, m_U[0]);
        }
        return
            m_w * (m_CAASF.Q(Q) + m_fpfdp)
          * sftbx::StructureFactor(SgOps, H, m_Coordinates, m_U);
      }
      //! XXX
      af::tiny<FloatType, 3>
      StructureFactor_dX(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q,
        const std::complex<FloatType>& dTarget_dFcalc) const
      {
        if (!m_Anisotropic) {
          return sftbx::StructureFactor_dX(
            SgOps, H, m_Coordinates,
              m_w
            * adptbx::DebyeWallerFactorUiso(Q / 4., m_U[0])
            * (m_CAASF.Q(Q) + m_fpfdp),
            dTarget_dFcalc);
        }
        throw cctbx_not_implemented();
      }
      /*! \brief Contribution of the (one) scatterer to an array of
          structure factors.
       */
      /*! The Miller indices H and d-spacing measures Q
          (cctbx::uctbx::UnitCell::Q())
          are passed as arrays. The contribution of the scatterer
          to each of the structure factors is added to the complex
          valued array Fcalc.
          <p>
          Requirements for all array types:
          <ul>
          <li>Must support the size() method.
          <li>Elements must be accessible with operator[].
          </ul>
          <p>
          See also: cctbx::StructureFactorArray()
       */
      template <typename MillerIndexArrayType,
                typename QArrayType,
                typename FcalcArrayType>
      void
      StructureFactorArray(
        const sgtbx::SpaceGroup& SgOps,
        const MillerIndexArrayType& H,
        const QArrayType& Q,
        FcalcArrayType Fcalc) const
      {
        CheckMultiplicity();
        cctbx_assert(Q.size() == H.size());
        cctbx_assert(Fcalc.size() == H.size());
        for (std::size_t i = 0; i < H.size(); i++) {
          Fcalc[i] += StructureFactor(SgOps, H[i], Q[i]);
        }
      }
      //! XXX
      template <typename MillerIndexArrayType,
                typename QArrayType,
                typename DerivativesArrayType>
      void
      StructureFactor_dX_Array(
        const sgtbx::SpaceGroup& SgOps,
        const MillerIndexArrayType& H,
        const QArrayType& Q,
        const DerivativesArrayType& dTarget_dFcalc,
        af::tiny<FloatType, 3>& dF_dX) const
      {
        CheckMultiplicity();
        cctbx_assert(Q.size() == H.size());
        cctbx_assert(dTarget_dFcalc.size() == H.size());
        for (std::size_t i = 0; i < H.size(); i++) {
          dF_dX += StructureFactor_dX(
            SgOps, H[i], Q[i], dTarget_dFcalc[i]);
        }
      }
    private:
      std::string m_Label;
      CAASF_Type m_CAASF;
      std::complex<FloatType> m_fpfdp;
      fractional<FloatType> m_Coordinates;
      FloatType m_Occ;
      bool m_Anisotropic;
      af::tiny<FloatType, 6> m_U;
      int m_M;
      FloatType m_w;
  };

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
# define CCTBX_tn_S_A_T_v_t_f_t typename SiteArrayType::value_type::float_type
#else
# define CCTBX_tn_S_A_T_v_t_f_t double
#endif

  //! XXX
  template <typename SiteArrayType>
  fractional<CCTBX_tn_S_A_T_v_t_f_t>
  least_squares_shift(const uctbx::UnitCell& ucell,
                      const SiteArrayType& sites1,
                      const SiteArrayType& sites2)
  {
    typedef CCTBX_tn_S_A_T_v_t_f_t float_type;
    cctbx_assert(sites1.size() == sites2.size());
    cartesian<float_type> sum_delta_c(0.,0.,0.);
    for(std::size_t i=0;i<sites1.size();i++) {
      fractional<float_type>
      delta_f = sites2[i].Coordinates() - sites1[i].Coordinates();
      cartesian<float_type> delta_c = ucell.orthogonalize(delta_f);
      sum_delta_c += delta_c;
    }
    cartesian<float_type> shift_c = sum_delta_c / float_type(sites1.size());
    return ucell.fractionalize(shift_c);
  }

  //! XXX
  template <typename SiteArrayType>
  CCTBX_tn_S_A_T_v_t_f_t
  rms_coordinates(
    const uctbx::UnitCell& ucell,
    const SiteArrayType& sites1,
    const SiteArrayType& sites2,
    const fractional<CCTBX_tn_S_A_T_v_t_f_t>& shift
        = fractional<CCTBX_tn_S_A_T_v_t_f_t>(0.,0.,0.))
  {
    typedef CCTBX_tn_S_A_T_v_t_f_t float_type;
    cctbx_assert(sites1.size() == sites2.size());
    std::vector<float_type> len2;
    len2.reserve(sites1.size());
    for(std::size_t i=0;i<sites1.size();i++) {
      fractional<float_type> x = sites2[i].Coordinates() - shift;
      len2.push_back(sites1[i].distance2(ucell, x));
    }
    return std::sqrt(af::mean(af::make_ref(len2)));
  }

  /*! \brief Compute structure factors for an array of
       Miller indices and an array of sites.
   */
  /*! The contribution for each site is added to each Fcalc.
      <p>
      The array of pre-computed d-spacing measures Q avoids the
      overhead of recomputing these values for each pass
      through the structure factor array (for each site).
      <p>
      Requirements for all array types:
      <ul>
      <li>Must support the size() method.
      <li>Elements must be accessible with operator[].
      </ul>
      <p>
      See also: XrayScatterer::StructureFactorArray()
   */
  template <typename MillerIndexArrayType,
            typename QArrayType,
            typename SiteArrayType,
            typename FcalcArrayType>
  void
  StructureFactorArray(const sgtbx::SpaceGroup& SgOps,
                       const MillerIndexArrayType& H,
                       const QArrayType& Q,
                       const SiteArrayType& Sites,
                       FcalcArrayType Fcalc)
  {
    for (std::size_t i = 0; i < Sites.size(); i++) {
      Sites[i].StructureFactorArray(SgOps, H, Q, Fcalc);
    }
  }

  /*! \brief Compute structure factors for an array of
       Miller indices and an array of sites.
   */
  /*! The contribution for each site is added to each Fcalc.
      <p>
      This functions pre-computed the d-spacing measures Q
      and calls the alternative overload.
      If the d-spacing measures Q are already available to
      the caller, use the alternative overload directly.
      <p>
      Requirements for all array types:
      <ul>
      <li>Must support the size() method.
      <li>Elements must be accessible with operator[].
      </ul>
      <p>
      See also: XrayScatterer::StructureFactorArray()
   */
  template <typename MillerIndexArrayType,
            typename SiteArrayType,
            typename FcalcArrayType>
  void
  StructureFactorArray(const uctbx::UnitCell& UC,
                       const sgtbx::SpaceGroup& SgOps,
                       const MillerIndexArrayType& H,
                       const SiteArrayType& Sites,
                       FcalcArrayType Fcalc)
  {
    af::shared<double> Q(H.size()); // FUTURE: avoid default initialization
    for (std::size_t i = 0; i < H.size(); i++) {
      Q[i] = UC.Q(H[i]);
    }
    StructureFactorArray(SgOps, H, Q, Sites, Fcalc);
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename QArrayType,
            typename DerivativesArrayType,
            typename SiteArrayType,
            typename DerivativesXArrayType>
  void
  StructureFactor_dX_Array(
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const QArrayType& Q,
    const DerivativesArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesXArrayType dF_dX)
  {
    cctbx_assert(Sites.size() == dF_dX.size());
    for (std::size_t i = 0; i < Sites.size(); i++) {
      Sites[i].StructureFactor_dX_Array(
        SgOps, H, Q, dTarget_dFcalc, dF_dX[i]);
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename DerivativesArrayType,
            typename SiteArrayType,
            typename DerivativesXArrayType>
  void
  StructureFactor_dX_Array(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const DerivativesArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesXArrayType dF_dX)
  {
    af::shared<double> Q(H.size()); // FUTURE: avoid default initialization
    for (std::size_t i = 0; i < H.size(); i++) {
      Q[i] = UC.Q(H[i]);
    }
    StructureFactor_dX_Array(
      SgOps, H, Q.const_ref(), dTarget_dFcalc, Sites, dF_dX);
  }

}} // namespace cctbx::sftbx

#endif // CCTBX_XRAY_SCATTERER_H
