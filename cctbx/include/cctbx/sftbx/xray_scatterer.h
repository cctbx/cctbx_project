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

  namespace detail {

    template <typename FloatType>
    struct null_accumulator
    {
      void inner_sum_reset(const Miller::Index&) {}
      void inner_sum_add(const std::complex<FloatType>&,
                         const Miller::Index&) {}
      void outer_sum_add() {}
      void finalize() {}
    };

    template <typename FloatType,
              typename AccumulatorType1,
              typename AccumulatorType2 = null_accumulator<FloatType>,
              typename AccumulatorType3 = null_accumulator<FloatType>,
              typename AccumulatorType4 = null_accumulator<FloatType>,
              typename AccumulatorType5 = null_accumulator<FloatType> >
    struct generic_sum_over_symmetry
    {
      static
      void
      run(const sgtbx::SpaceGroup& SgOps,
          const Miller::Index& H,
          const fractional<FloatType>& X,
          AccumulatorType1* accu1,
          AccumulatorType2* accu2 = 0,
          AccumulatorType3* accu3 = 0,
          AccumulatorType4* accu4 = 0,
          AccumulatorType5* accu5 = 0)
      {
        using constants::pi;
        for (int s = 0; s < SgOps.nSMx(); s++) {
          Miller::Index HR = H * SgOps[s].Rpart();
          FloatType HRX = HR * X;
          sgtbx::TrVec T = SgOps[s].Tpart();
                     accu1->inner_sum_reset(HR);
          if (accu2) accu2->inner_sum_reset(HR);
          if (accu3) accu3->inner_sum_reset(HR);
          if (accu4) accu4->inner_sum_reset(HR);
          if (accu5) accu5->inner_sum_reset(HR);
          for (int i = 0; i < SgOps.fInv(); i++) {
            if (i) {
              HR = -HR;
              HRX = -HRX;
              T = SgOps.InvT() - T;
            }
            for (int l = 0; l < SgOps.nLTr(); l++) {
              FloatType HT = FloatType(H * (T + SgOps.LTr(l))) / SgOps.TBF();
              FloatType phase = FloatType(2 * pi) * (HRX + HT);
              std::complex<FloatType> e2piiHXs(
                std::cos(phase), std::sin(phase));
                         accu1->inner_sum_add(e2piiHXs, HR);
              if (accu2) accu2->inner_sum_add(e2piiHXs, HR);
              if (accu3) accu3->inner_sum_add(e2piiHXs, HR);
              if (accu4) accu4->inner_sum_add(e2piiHXs, HR);
              if (accu5) accu5->inner_sum_add(e2piiHXs, HR);
            }
          }
                     accu1->outer_sum_add();
          if (accu2) accu2->outer_sum_add();
          if (accu3) accu3->outer_sum_add();
          if (accu4) accu4->outer_sum_add();
          if (accu5) accu5->outer_sum_add();
        }
                   accu1->finalize();
        if (accu2) accu2->finalize();
        if (accu3) accu3->finalize();
        if (accu4) accu4->finalize();
        if (accu5) accu5->finalize();
      }
    };

    template <typename FloatType>
    struct structure_factor_plain_accumulator : null_accumulator<FloatType>
    {
      structure_factor_plain_accumulator()
      : outer_sum(0)
      {}
      void inner_sum_add(const std::complex<FloatType>& e2piiHXs,
                         const Miller::Index&) {
        outer_sum += e2piiHXs;
      }
      std::complex<FloatType> outer_sum;
    };

    template <typename FloatType>
    struct structure_factor_aniso_accumulator : null_accumulator<FloatType>
    {
      structure_factor_aniso_accumulator(const af::tiny<FloatType, 6>& Ustar)
      : Ustar_(Ustar),
        outer_sum(0)
      {}
      void inner_sum_reset(const Miller::Index& HR) {
        HR_ = HR;
        inner_sum = std::complex<FloatType>(0);
      }
      void inner_sum_add(const std::complex<FloatType>& e2piiHXs,
                         const Miller::Index&) {
        inner_sum += e2piiHXs;
      }
      void outer_sum_add() {
        outer_sum += inner_sum * adptbx::DebyeWallerFactorUstar(HR_, Ustar_);
      }
      const af::tiny<FloatType, 6>& Ustar_;
      Miller::Index HR_;
      std::complex<FloatType> inner_sum;
      std::complex<FloatType> outer_sum;
    };

    template <typename FloatType>
    struct dTarget_dX_accumulator : null_accumulator<FloatType>
    {
      dTarget_dX_accumulator(
        const std::complex<FloatType>& phase_indep_coeff,
        const std::complex<FloatType>& dTarget_dFcalc)
      : pic_(phase_indep_coeff),
        dT_dF_(dTarget_dFcalc)
      {
        outer_sum.fill(FloatType(0));
      }
      void inner_sum_add(const std::complex<FloatType>& e2piiHXs,
                         const Miller::Index& HR) {
        std::complex<FloatType> f = pic_ * e2piiHXs;
        FloatType c =   f.real() * dT_dF_.imag()
                      - f.imag() * dT_dF_.real();
        for(int i=0;i<3;i++) {
          outer_sum[i] += FloatType(HR[i]) * c;
        }
      }
      void finalize() {
        outer_sum *= FloatType(2 * constants::pi);
      }
      const std::complex<FloatType>& pic_;
      const std::complex<FloatType>& dT_dF_;
      af::tiny<FloatType, 3> outer_sum;
    };

    template <typename FloatType>
    struct dTarget_dS_accumulator : null_accumulator<FloatType>
    {
      dTarget_dS_accumulator(
        const std::complex<FloatType>& phase_indep_coeff,
        const std::complex<FloatType>& dTarget_dFcalc)
      : pic_(phase_indep_coeff),
        dT_dF_(dTarget_dFcalc),
        outer_sum(0)
      {}
      void inner_sum_add(const std::complex<FloatType>& e2piiHXs,
                         const Miller::Index& HR) {
        std::complex<FloatType> f = pic_ * e2piiHXs;
        outer_sum +=   f.real() * dT_dF_.real()
                     + f.imag() * dT_dF_.imag();
      }
      const std::complex<FloatType>& pic_;
      const std::complex<FloatType>& dT_dF_;
      FloatType outer_sum;
    };

  } // namespace detail

  //! Structure factor without Debye-Waller factor.
  /*! Sum of exp(2 pi j H S X) over all symmetry operations S.
      j is the imaginary number.
   */
  template <typename FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  const Miller::Index& H,
                  const fractional<FloatType>& X)
  {
    typedef detail::structure_factor_plain_accumulator<FloatType> accu_type;
    accu_type accu;
    detail::generic_sum_over_symmetry<FloatType, accu_type>::run(
      SgOps, H, X, &accu);
    return accu.outer_sum;
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
  template <typename FloatType>
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
  template <typename FloatType>
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
  template <typename FloatType>
  std::complex<FloatType>
  StructureFactor(const sgtbx::SpaceGroup& SgOps,
                  const Miller::Index& H,
                  const fractional<FloatType>& X,
                  const af::tiny<FloatType, 6>& Ustar)
  {
    typedef detail::structure_factor_aniso_accumulator<FloatType> accu_type;
    accu_type accu(Ustar);
    detail::generic_sum_over_symmetry<FloatType, accu_type>::run(
      SgOps, H, X, &accu);
    return accu.outer_sum;
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
      //! Constructor with isotropic displacement parameter.
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
      //! Constructor with anisotropic displacement parameters.
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
      //! Redefinition of the coordinates.
      /*! It is the user's responsibility to ensure that the
          multiplicity is still correct or updated.
       */
      void set_Coordinates(const fractional<FloatType>& Coordinates) {
        m_Coordinates = Coordinates;
      }
      //! Redefinition of the occupancy.
      /*! This overload may only be used before ApplySymmetry() is called.
       */
      void set_Occ(const FloatType& Occ) {
        cctbx_assert(m_M == 0);
        m_Occ = Occ;
      }
      //! Redefinition of the occupancy.
      /*! This overload must be used after ApplySymmetry() was called.
          w() will be updated along with Occ().
       */
      void set_Occ(const FloatType& Occ, const sgtbx::SpaceGroup& SgOps) {
        m_Occ = Occ;
        m_w = m_Occ * m_M / SgOps.OrderZ();
      }
      //! Redefinition of the isotropic displacement parameter.
      /*! Requires isAnisotropic() == <code>false</code>.
       */
      void set_Uiso(const FloatType& Uiso) {
        cctbx_assert(!m_Anisotropic);
        m_U[0] = Uiso;
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
      sgtbx::SiteSymmetry
      ApplySymmetry(const uctbx::UnitCell& UC,
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
        return SS;
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
      StructureFactor_dT_dX(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q,
        const std::complex<FloatType>& dTarget_dFcalc) const
      {
        if (!m_Anisotropic) {
          typedef detail::dTarget_dX_accumulator<FloatType> accu_type;
          accu_type accu(
              m_w
            * adptbx::DebyeWallerFactorUiso(Q / 4., m_U[0])
            * (m_CAASF.Q(Q) + m_fpfdp),
            dTarget_dFcalc);
          detail::generic_sum_over_symmetry<FloatType, accu_type>::run(
            SgOps, H, m_Coordinates, &accu);
          return accu.outer_sum;
        }
        throw cctbx_not_implemented();
      }

      //! XXX
      FloatType
      StructureFactor_dT_dOcc(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q,
        const std::complex<FloatType>& dTarget_dFcalc) const
      {
        if (!m_Anisotropic) {
          typedef detail::dTarget_dS_accumulator<FloatType> accu_type;
          accu_type accu(
              FloatType(m_M) / SgOps.OrderZ()
            * adptbx::DebyeWallerFactorUiso(Q / 4., m_U[0])
            * (m_CAASF.Q(Q) + m_fpfdp),
            dTarget_dFcalc);
          detail::generic_sum_over_symmetry<FloatType, accu_type>::run(
            SgOps, H, m_Coordinates, &accu);
          return accu.outer_sum;
        }
        throw cctbx_not_implemented();
      }

      //! XXX
      FloatType
      StructureFactor_dT_dUiso(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q,
        const std::complex<FloatType>& dTarget_dFcalc) const
      {
        if (!m_Anisotropic) {
          typedef detail::dTarget_dS_accumulator<FloatType> accu_type;
          accu_type accu(
              m_w
            * adptbx::DebyeWallerFactorUiso(Q / 4., m_U[0])
            * (m_CAASF.Q(Q) + m_fpfdp),
            dTarget_dFcalc);
          detail::generic_sum_over_symmetry<FloatType, accu_type>::run(
            SgOps, H, m_Coordinates, &accu);
          return FloatType(-adptbx::U_as_B(1.) * Q / 4.) * accu.outer_sum;
        }
        throw cctbx_not_implemented();
      }

      //! XXX
      void
      StructureFactor_dT_dX_dUiso(
        const sgtbx::SpaceGroup& SgOps,
        const Miller::Index& H,
        double Q,
        const std::complex<FloatType>& dTarget_dFcalc,
        af::tiny<FloatType, 3>& dT_dX,
        FloatType& dT_dUiso) const
      {
        if (!m_Anisotropic) {
          std::complex<FloatType> phase_indep_coeff(
              m_w
            * adptbx::DebyeWallerFactorUiso(Q / 4., m_U[0])
            * (m_CAASF.Q(Q) + m_fpfdp));
          typedef detail::dTarget_dX_accumulator<FloatType> accu_dX_type;
          typedef detail::dTarget_dS_accumulator<FloatType> accu_dU_type;
          accu_dX_type accu_dX(phase_indep_coeff, dTarget_dFcalc);
          accu_dU_type accu_dU(phase_indep_coeff, dTarget_dFcalc);
          detail::generic_sum_over_symmetry<
            FloatType, accu_dX_type, accu_dU_type>::run(
              SgOps, H, m_Coordinates, &accu_dX, &accu_dU);
          dT_dX += accu_dX.outer_sum;
          dT_dUiso +=   FloatType(-adptbx::U_as_B(1.) * Q / 4.)
                      * accu_dU.outer_sum;
          return;
        }
        throw cctbx_not_implemented();
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

  /*! \brief Computation of structure factors for an array of
       Miller indices and an array of sites.
   */
  /*! The contribution of each site is added to each Fcalc.
      <p>
      Requirements for all array types:
      <ul>
      <li>Must support the size() method.
      <li>Elements must be accessible with operator[].
      </ul>
      <p>
      See also: cctbx::sftbx::XrayScatterer::StructureFactor()
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
    cctbx_assert(H.size() == Fcalc.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      double Q = UC.Q(H[i]);
      for (std::size_t j = 0; j < Sites.size(); j++) {
        Fcalc[i] += Sites[j].StructureFactor(SgOps, H[i], Q);
      }
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename DerivativesFcalcArrayType,
            typename SiteArrayType,
            typename DerivativesXArrayType>
  void
  StructureFactor_dT_dX_Array(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const DerivativesFcalcArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesXArrayType dT_dX)
  {
    cctbx_assert(Sites.size() == dT_dX.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      double Q = UC.Q(H[i]);
      for (std::size_t j = 0; j < Sites.size(); j++) {
        dT_dX[j] += Sites[j].StructureFactor_dT_dX(
          SgOps, H[i], Q, dTarget_dFcalc[i]);
      }
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename DerivativesFcalcArrayType,
            typename SiteArrayType,
            typename DerivativesOccArrayType>
  void
  StructureFactor_dT_dOcc_Array(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const DerivativesFcalcArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesOccArrayType dT_dOcc)
  {
    cctbx_assert(Sites.size() == dT_dOcc.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      double Q = UC.Q(H[i]);
      for (std::size_t j = 0; j < Sites.size(); j++) {
        dT_dOcc[j] += Sites[j].StructureFactor_dT_dOcc(
          SgOps, H[i], Q, dTarget_dFcalc[i]);
      }
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename DerivativesFcalcArrayType,
            typename SiteArrayType,
            typename DerivativesUisoArrayType>
  void
  StructureFactor_dT_dUiso_Array(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const DerivativesFcalcArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesUisoArrayType dT_dUiso)
  {
    cctbx_assert(Sites.size() == dT_dUiso.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      double Q = UC.Q(H[i]);
      for (std::size_t j = 0; j < Sites.size(); j++) {
        dT_dUiso[j] += Sites[j].StructureFactor_dT_dUiso(
          SgOps, H[i], Q, dTarget_dFcalc[i]);
      }
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename DerivativesFcalcArrayType,
            typename SiteArrayType,
            typename DerivativesXArrayType,
            typename DerivativesUisoArrayType>
  void
  StructureFactor_dT_dX_dUiso_Array(
    const uctbx::UnitCell& UC,
    const sgtbx::SpaceGroup& SgOps,
    const MillerIndexArrayType& H,
    const DerivativesFcalcArrayType& dTarget_dFcalc,
    const SiteArrayType& Sites,
    DerivativesXArrayType dT_dX,
    DerivativesUisoArrayType dT_dUiso)
  {
    cctbx_assert(Sites.size() == dT_dX.size());
    cctbx_assert(Sites.size() == dT_dUiso.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      double Q = UC.Q(H[i]);
      for (std::size_t j = 0; j < Sites.size(); j++) {
        Sites[j].StructureFactor_dT_dX_dUiso(
          SgOps, H[i], Q, dTarget_dFcalc[i],
          dT_dX[j], dT_dUiso[j]);
      }
    }
  }

}} // namespace cctbx::sftbx

#endif // CCTBX_XRAY_SCATTERER_H
