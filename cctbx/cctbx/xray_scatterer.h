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

namespace cctbx {

  template <class FloatType, class CAASF_Type>
  class XrayScatterer : public fractional<FloatType>
  {
    public:
      XrayScatterer() {}
      XrayScatterer(const std::string& Label,
                    const CAASF_Type& CAASF,
                    const std::complex<FloatType>& fpfdp,
                    const fractional<FloatType>& Coordinates,
                    const FloatType& Occ,
                    const FloatType& Uiso)
        : m_Label(Label),
          m_CAASF(CAASF),
          m_fpfdp(fpfdp),
          fractional<FloatType>(Coordinates),
          m_Occ(Occ),
          m_Anisotropic(false),
          m_M(0),
          m_w(0)
      {
        m_U.assign(0.);
        m_U[0] = Uiso;
      }
      XrayScatterer(const std::string& Label,
                    const CAASF_Type& CAASF,
                    const std::complex<FloatType>& fpfdp,
                    const fractional<FloatType>& Coordinates,
                    const FloatType& Occ,
                    const boost::array<FloatType, 6>& Uaniso)
        : m_Label(Label),
          m_CAASF(CAASF),
          m_fpfdp(fpfdp),
          fractional<FloatType>(Coordinates),
          m_Occ(Occ),
          m_Anisotropic(true),
          m_U(Uaniso),
          m_M(0),
          m_w(0)
      {
      }
      inline const std::string& Label() const { return m_Label; }
      inline const CAASF_Type& CAASF() const { return m_CAASF; }
      inline const std::complex<FloatType>& fpfdp() const { return m_fpfdp; }
      inline const fractional<FloatType>& Coordinates() const {
        return static_cast<const fractional<FloatType>&>(*this);
      }
      inline bool isAnisotropic() const { return m_Anisotropic; }
      inline const FloatType& Uiso() const { return m_U[0]; }
      inline const boost::array<FloatType, 6>& Uaniso() { return m_U; }
      inline int M() const { return m_M; }
      inline const FloatType& w() const { return m_w; }
      void ApplySymmetry(const uctbx::UnitCell& UC,
                         const sgtbx::SpaceGroup& SgOps,
                         double MinMateDistance = 0.5,
                         double Ustar_tolerance = 0.1)
      {
        sgtbx::SpecialPositionSnapParameters
        SnapParameters(UC, SgOps, true, MinMateDistance);
        sgtbx::SiteSymmetry
        SS(SnapParameters, *this);
        for(std::size_t i=0;i<3;i++) elems[i] = SS.SnapPosition()[i];
        m_M = SS.M();
        m_w = m_Occ * m_M / SgOps.OrderZ();
        if (m_Anisotropic) {
          if (Ustar_tolerance > 0.) {
            SS.CheckUstar(m_U, Ustar_tolerance);
          }
          m_U = SS.AverageUstar(m_U);
          adptbx::CheckPositiveDefinite(m_U);
        }
      }
      inline std::complex<FloatType>
      StructureFactor(const sgtbx::SpaceGroup& SgOps,
                      const Miller::Index& H,
                      double Q) const
      {
        if (!m_Anisotropic) {
          return
              m_w * (m_CAASF.Q(Q) + m_fpfdp)
            * SgOps.StructureFactor(Q / 4., H, *this, m_U[0]);
        }
        return
            m_w * (m_CAASF.Q(Q) + m_fpfdp)
          * SgOps.StructureFactor(H, *this, m_U);
      }
      template <class MillerVectorType,
                class doubleVectorType,
                class StdComplexVectorType>
      inline void
      StructureFactorVector(const sgtbx::SpaceGroup& SgOps,
                            const MillerVectorType& H,
                            const doubleVectorType& Q,
                            StdComplexVectorType& Fcalc) const
      {
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
      FloatType m_Occ;
      bool m_Anisotropic;
      boost::array<FloatType, 6> m_U;
      int m_M;
      FloatType m_w;
  };

  template <class MillerVectorType,
            class doubleVectorType,
            class XrayScattererVectorType,
            class StdComplexVectorType>
  inline void
  StructureFactorVector(const sgtbx::SpaceGroup& SgOps,
                        const MillerVectorType& H,
                        const doubleVectorType& Q,
                        const XrayScattererVectorType& Sites,
                        StdComplexVectorType& Fcalc)
  {
    for (std::size_t i = 0; i < Sites.size(); i++) {
      Sites[i].StructureFactorVector(SgOps, H, Q, Fcalc);
    }
  }

  template <class MillerVectorType,
            class XrayScattererVectorType,
            class StdComplexVectorType>
  inline void
  StructureFactorVector(const uctbx::UnitCell& UC,
                        const sgtbx::SpaceGroup& SgOps,
                        const MillerVectorType& H,
                        const XrayScattererVectorType& Sites,
                        StdComplexVectorType& Fcalc)
  {
    std::vector<double> Q(H.size());
    for (std::size_t i = 0; i < H.size(); i++) {
      Q[i] = UC.Q(H[i]);
    }
    StructureFactorVector(SgOps, H, Q, Sites, Fcalc);
  }

} // namespace cctbx

#endif // CCTBX_XRAY_SCATTERER_H
