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
#include <cctbx/eltbx/caasf.h>
#include <cctbx/adptbx.h>

namespace cctbx {

  template <class FloatType, std::size_t CAASF_N>
  class XrayScatterer : public fractional<FloatType>
  {
    public:
      XrayScatterer() {}
      XrayScatterer(const std::string& Label,
                    const eltbx::CAASF<CAASF_N>& CAASF,
                    const fractional<FloatType>& Coordinates,
                    const FloatType& Occ,
                    const FloatType& Uiso,
                    const std::complex<FloatType>& fpfdp)
        : m_Label(Label),
          m_CAASF(CAASF),
          fractional<FloatType>(Coordinates),
          m_Occ(Occ),
          m_Anisotropic(false),
          m_fpfdp(fpfdp),
          m_M(0),
          m_w(0)
      {
        m_U.assign(0.);
        m_U[0] = Uiso;
      }
      void DetermineMultiplicity(const uctbx::UnitCell& UC,
                                 const sgtbx::SpaceGroup& SgOps,
                                 double MinMateDistance = 0.5)
      {
        sgtbx::SpecialPositionSnapParameters
        SnapParameters(UC, SgOps, true, MinMateDistance);
        sgtbx::SiteSymmetry
        SS(SnapParameters, *this);
        for(std::size_t i=0;i<3;i++) elems[i] = SS.SnapPosition()[i];
        m_M = SS.M();
        m_w = m_Occ * m_M / SgOps.OrderZ();
      }
    private:
      std::string m_Label;
      eltbx::CAASF<CAASF_N> m_CAASF;
      FloatType m_Occ;
      bool m_Anisotropic;
      boost::array<FloatType, 6> m_U;
      std::complex<FloatType> m_fpfdp;
      FloatType m_M;
      FloatType m_w;
  };

  template <class FloatType, std::size_t CAASF_N>
  std::complex<FloatType>
  inline StructureFactor(const uctbx::UnitCell& UC,
                         const sgtbx::SpaceGroup& SgOps,
                         const XrayScatterer<FloatType, CAASF_N>& Site)
  {
    return std::complex<FloatType>(0);
  }

} // namespace cctbx

#endif // CCTBX_XRAY_SCATTERER_H
