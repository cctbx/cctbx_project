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

namespace cctbx {
  //! Structure Factor Toolbox namespace.
  namespace sftbx {

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
      //! Default constructor. Data members are not initialized!
      XrayScatterer() {}
      //! Constructor with isotropic Debye-Waller factor.
      XrayScatterer(const std::string& Label,
                    const CAASF_Type& CAASF,
                    const std::complex<FloatType>& fpfdp,
                    const cctbx::fractional<FloatType>& Coordinates,
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
                    const cctbx::fractional<FloatType>& Coordinates,
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
      {
      }
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
      const cctbx::fractional<FloatType>& Coordinates() const {
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
      void ApplySymmetry(const cctbx::uctbx::UnitCell& UC,
                         const cctbx::sgtbx::SpaceGroup& SgOps,
                         double MinMateDistance = 0.5,
                         double Ustar_tolerance = 0.1,
                         bool TestPositiveDefiniteness = true)
      {
        cctbx::sgtbx::SpecialPositionSnapParameters
        SnapParameters(UC, SgOps, true, MinMateDistance);
        cctbx::sgtbx::SiteSymmetry
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
            cctbx::adptbx::CheckPositiveDefinite(m_U);
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
      /*! \brief Contribution of the (one) scatterer to the
          structure factor with the Miller index H.
       */
      /*! Q is a d-spacing measure (cctbx::uctbx::UnitCell::Q()).
       */
      std::complex<FloatType>
      StructureFactor(const cctbx::sgtbx::SpaceGroup& SgOps,
                      const cctbx::Miller::Index& H,
                      double Q) const
      {
        if (!m_Anisotropic) {
          return
              m_w * (m_CAASF.Q(Q) + m_fpfdp)
            * SgOps.StructureFactor(Q / 4., H, m_Coordinates, m_U[0]);
        }
        return
            m_w * (m_CAASF.Q(Q) + m_fpfdp)
          * SgOps.StructureFactor(H, m_Coordinates, m_U);
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
                typename doubleArrayType,
                typename StdComplexArrayType>
      void
      StructureFactorArray(const cctbx::sgtbx::SpaceGroup& SgOps,
                            const MillerIndexArrayType& H,
                            const doubleArrayType& Q,
                            StdComplexArrayType Fcalc) const
      {
        if (m_M == 0) {
          throw cctbx::error(
            "ApplySymmetry() has not been called for this scatterer.");
        }
        cctbx_assert(Q.size() == H.size());
        cctbx_assert(Fcalc.size() == H.size());
        for (std::size_t i = 0; i < H.size(); i++) {
          Fcalc[i] += StructureFactor(SgOps, H[i], Q[i]);
        }
      }
    private:
      std::string m_Label;
      CAASF_Type m_CAASF;
      std::complex<FloatType> m_fpfdp;
      cctbx::fractional<FloatType> m_Coordinates;
      FloatType m_Occ;
      bool m_Anisotropic;
      af::tiny<FloatType, 6> m_U;
      int m_M;
      FloatType m_w;
  };

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
            typename doubleArrayType,
            typename XrayScattererArrayType,
            typename StdComplexArrayType>
  void
  StructureFactorArray(const cctbx::sgtbx::SpaceGroup& SgOps,
                       const MillerIndexArrayType& H,
                       const doubleArrayType& Q,
                       const XrayScattererArrayType& Sites,
                       StdComplexArrayType Fcalc)
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
            typename XrayScattererArrayType,
            typename StdComplexArrayType>
  void
  StructureFactorArray(const cctbx::uctbx::UnitCell& UC,
                       const cctbx::sgtbx::SpaceGroup& SgOps,
                       const MillerIndexArrayType& H,
                       const XrayScattererArrayType& Sites,
                       StdComplexArrayType Fcalc)
  {
    af::shared<double> Q(H.size()); // FUTURE: avoid default initialization
    for (std::size_t i = 0; i < H.size(); i++) {
      Q[i] = UC.Q(H[i]);
    }
    StructureFactorArray(SgOps, H, Q, Sites, Fcalc);
  }

}} // namespace cctbx::sftbx

#endif // CCTBX_XRAY_SCATTERER_H
