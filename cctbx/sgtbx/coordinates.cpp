// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 12: SpecialPosition -> SiteSymmetry (R.W. Grosse-Kunstleve)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cmath>
#include <algorithm>
#include <cctbx/sgtbx/coordinates.h>

namespace cctbx { namespace sgtbx {

  namespace detail {

    void SetUniqueOps(const SpaceGroup& SgOps,
                      const RTMx& SpecialOp,
                      std::vector<RTMx>& UniqueOps)
    {
      for(int i=0;i<SgOps.OrderZ();i++) {
        RTMx SS = SgOps(i).multiply(SpecialOp).modPositive();
        if (std::find(UniqueOps.begin(),
                      UniqueOps.end(), SS) == UniqueOps.end()) {
          UniqueOps.push_back(SS);
        }
      }
    }

    class CloseMate {
      public:
        CloseMate(const RTMx& M, double CartDelta2) {
          this->M = M;
          this->CartDelta2 = CartDelta2;
        }
        RTMx M;
        double CartDelta2;
    };

    class CmpCloseMates {
      public:
        bool operator()(const CloseMate& a, const CloseMate& b) {
          if (a.CartDelta2 < b.CartDelta2) return true;
          return false;
        }
    };

    void RT_PointGroup::reset(const RTMx& M)
    {
      Matrices.clear();
      Matrices.push_back(M);
    }

    void RT_PointGroup::add(const RTMx& M)
    {
      for (vec_type::const_iterator
           Mi = Matrices.begin(); Mi != Matrices.end(); Mi++) {
        if (Mi->Rpart() == M.Rpart()) {
          if (Mi->Tpart() != M.Tpart()) invalid = true;
          return;
        }
      }
      Matrices.push_back(M);
    }

    void RT_PointGroup::expand(const RTMx& M)
    {
      RTMx TrialM = M;
      int i = Matrices.size();
      int j = 1;
      for (;;) {
        add(TrialM);
        if (invalid) return;
        if (j > i) {
          i++;
          j = 1;
        }
        if (i == Matrices.size()) break;
        TrialM = Matrices[j] * Matrices[i];
        j++;
      }
    }

    bool RT_PointGroup::try_expand(const RTMx& M)
    {
      vec_type::size_type size_before_expand = Matrices.size();
      expand(M);
      if (invalid) {
        Matrices.resize(size_before_expand);
        invalid = false;
        return false;
      }
      return true;
    }

    RTMx RT_PointGroup::accumulate() const
    {
      RTMx result = Matrices[0];
      for(int i=1;i<Matrices.size();i++) result += Matrices[i];
      result.pseudo_divide(Matrices.size());
      return result.cancel();
    }

  } // namespace detail

  void SiteSymmetry::BuildSpecialOp()
  {
    const uctbx::UnitCell& uc = m_Parameters->m_UnitCell;
    const SpaceGroup& SgOps = m_Parameters->m_SgOps;
    int TBF = SgOps.TBF();
    m_ShortestDistance2 = uc.getLongestVector2();
    std::vector<detail::CloseMate> CloseMates;
    CloseMates.push_back(detail::CloseMate(SgOps[0], 0.));
    for (int iOp = 1; iOp < SgOps.OrderP(); iOp++) {
      RTMx M = SgOps(iOp);
      RotMx CumR = M.Rpart().accumulate();
      fractional<double> Mate0 = M * m_SnapPosition;
      for (int iLTr = 0; iLTr < SgOps.nLTr(); iLTr++) {
        fractional<double> Mate = Mate0 + SgOps.LTr(iLTr);
        fractional<double> Delta0 = Mate - m_SnapPosition;
        Delta0 = Delta0.modShort();
        TrVec
        UShifts = detail::getUnitShifts((m_SnapPosition + Delta0) - Mate);
        UShifts = UShifts.scale(TBF);
        TrVec MT = M.Tpart() + SgOps.LTr(iLTr) + UShifts;
        bool special = false;
        double SD2 = m_ShortestDistance2;
        for (UShifts[0] = -TBF; UShifts[0] <= TBF; UShifts[0] += TBF)
        for (UShifts[1] = -TBF; UShifts[1] <= TBF; UShifts[1] += TBF)
        for (UShifts[2] = -TBF; UShifts[2] <= TBF; UShifts[2] += TBF) {
          fractional<double> Delta = Delta0 + UShifts;
          double CartDelta2 = uc.Length2(Delta);
          if (SD2 > CartDelta2) SD2 = CartDelta2;
          if (CartDelta2 <= m_Parameters->m_MinMateDistance2) {
            TrVec MTU = MT + UShifts;
            TrVec IntrinsicPart = CumR * MTU;
            if (static_cast<Vec3>(IntrinsicPart) == 0) {
              CloseMates.push_back(
                detail::CloseMate(RTMx(M.Rpart(), MTU), CartDelta2));
              special = true;
            }
          }
        }
        if (!special) m_ShortestDistance2 = SD2;
      }
    }
    cctbx_assert(CloseMates.size() > 0);
    if (CloseMates.size() > 1) {
      std::sort(CloseMates.begin() + 1, CloseMates.end(),
                detail::CmpCloseMates());
    }
    m_PointGroup.reset(CloseMates[0].M);
    for(std::size_t i=1;i<CloseMates.size();i++) {
      if (!m_PointGroup.try_expand(CloseMates[i].M)) {
        if (m_ShortestDistance2 > CloseMates[i].CartDelta2)
            m_ShortestDistance2 = CloseMates[i].CartDelta2;
      }
    }
    cctbx_assert(SgOps.OrderZ() % m_PointGroup.Matrices.size() == 0);
    m_M = SgOps.OrderZ() / m_PointGroup.Matrices.size();
    m_SpecialOp = m_PointGroup.accumulate();
    m_SnapPosition = m_SpecialOp * m_OriginalPosition;
  }

  SiteSymmetry::SiteSymmetry(const SpecialPositionSnapParameters& params,
                             const fractional<double>& X,
                             bool auto_expand)

    : m_Parameters(&params),
      m_OriginalPosition(X),
      m_SnapPosition(X),
      m_ShortestDistance2(-1.),
      m_M(0),
      m_SpecialOp(0, 0)
  {
    RTMx LastSpecialOp(1, 1);
    for (;;) {
      BuildSpecialOp();
      if (m_SpecialOp == LastSpecialOp) break;
      LastSpecialOp = m_SpecialOp;
    }
    if (m_Parameters->m_MustBeWellBehaved && !isWellBehaved()) {
      throw error("SiteSymmetry: MinMateDistance too large.");
    }
    if (auto_expand) expand();
  }

  tables::MatrixGroup::Code SiteSymmetry::PointGroupType() const
  {
    SpaceGroup SiteSgOps;
    for(std::size_t i=0;i<m_PointGroup.Matrices.size();i++) {
      SiteSgOps.expandSMx(m_PointGroup.Matrices[i]);
    }
    cctbx_assert(SiteSgOps.nLTr() == 1);
    return SiteSgOps.getPointGroupType();
  }

}} // namespace cctbx::sgtbx
