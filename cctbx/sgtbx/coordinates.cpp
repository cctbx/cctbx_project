// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cmath>
#include <algorithm>
#include <cctbx/sgtbx/coordinates.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  namespace detail {

    uctbx::Vec3 ModPositive(const uctbx::Vec3& X)
    {
      uctbx::Vec3 result;
      rangei(3) {
        result[i] = std::fmod(X[i], 1.);
        while (result[i] <  0.) result[i] += 1.;
        while (result[i] >= 1.) result[i] -= 1.;
      }
      return result;
    }

    uctbx::Vec3 ModShort(const uctbx::Vec3& X)
    {
      uctbx::Vec3 result;
      rangei(3) {
        result[i] = std::fmod(X[i], 1.);
        if      (result[i] <= -.5) result[i] += 1.;
        else if (result[i] >   .5) result[i] -= 1.;
      }
      return result;
    }

    TrVec getUnitShifts(const uctbx::Vec3& Delta, int TBF)
    {
      TrVec result(TBF);
      rangei(3) {
        if (Delta[i] >= 0.) result[i] = TBF * static_cast<int>(Delta[i] + 0.5);
        else                result[i] = TBF * static_cast<int>(Delta[i] - 0.5);
      }
      return result;
    }

    void SetUniqueOps(const SgOps& sgo,
                      const RTMx& SpecialOp,
                      std::vector<RTMx>& UniqueOps)
    {
      rangei(sgo.OrderZ()) {
        RTMx SS = sgo(i).multiply(SpecialOp);
        SS.ModPositive();
        if (std::find(UniqueOps.begin(),
                      UniqueOps.end(), SS) == UniqueOps.end()) {
          UniqueOps.push_back(SS);
        }
      }
    }

    double ModLength2(const uctbx::UnitCell& uc, const uctbx::Vec3& Diff) {
      return uc.Length2(ModShort(Diff));
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
        inline bool operator()(const CloseMate& a, const CloseMate& b) {
          if (a.CartDelta2 < b.CartDelta2) return true;
          return false;
        }
    };

    class RT_PointGroup {
      public:
        bool invalid;
        typedef std::vector<RTMx> vec_type;
        vec_type Matrices;
        RT_PointGroup(const RTMx& M) : invalid(false), Matrices(1, M) {}
        void add(const RTMx& M);
        void expand(const RTMx& M);
        bool try_expand(const RTMx& M);
        RTMx accumulate() const;
    };

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

  void SpecialPosition::BuildSpecialOp()
  {
    const uctbx::UnitCell& uc = m_Parameters.m_UnitCell;
    const SgOps& sgo = m_Parameters.m_SgOps;
    int TBF = sgo.TBF();
    m_ShortestDistance2 = uc.getLongestVector2();
    std::vector<detail::CloseMate> CloseMates;
    CloseMates.push_back(detail::CloseMate(sgo[0], 0.));
    for (int iOp = 1; iOp < sgo.OrderP(); iOp++) {
      RTMx M = sgo(iOp);
      RotMx CumR = M.Rpart().accumulate();
      uctbx::Vec3 Mate0 = M * m_SnapPosition;
      for (int iLTr = 0; iLTr < sgo.nLTr(); iLTr++) {
        uctbx::Vec3 Mate = Mate0 + sgo.LTr(iLTr);
        uctbx::Vec3 Delta0 = detail::ModShort(Mate - m_SnapPosition);
        TrVec
        UShifts = detail::getUnitShifts((m_SnapPosition + Delta0) - Mate,TBF);
        TrVec MT = M.Tpart() + sgo.LTr(iLTr) + UShifts;
        bool special = false;
        double SD2 = m_ShortestDistance2;
        for (UShifts[0] = -TBF; UShifts[0] <= TBF; UShifts[0] += TBF)
        for (UShifts[1] = -TBF; UShifts[1] <= TBF; UShifts[1] += TBF)
        for (UShifts[2] = -TBF; UShifts[2] <= TBF; UShifts[2] += TBF) {
          uctbx::Vec3 Delta = Delta0 + UShifts;
          double CartDelta2 = uc.Length2(Delta);
          if (SD2 > CartDelta2) SD2 = CartDelta2;
          if (CartDelta2 <= m_Parameters.m_MinMateDistance2) {
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
    detail::RT_PointGroup PG(CloseMates[0].M);
    for(int i=1;i<CloseMates.size();i++) {
      if (!PG.try_expand(CloseMates[i].M)) {
        if (m_ShortestDistance2 > CloseMates[i].CartDelta2)
            m_ShortestDistance2 = CloseMates[i].CartDelta2;
      }
    }
    cctbx_assert(sgo.OrderZ() % PG.Matrices.size() == 0);
    m_M = sgo.OrderZ() / PG.Matrices.size();
    m_SpecialOp = PG.accumulate();
    m_SnapPosition = m_SpecialOp * m_OriginalPosition;
  }

  SpecialPosition::SpecialPosition(const SpecialPositionSnapParameters& params,
                                   const uctbx::Vec3& X,
                                   bool auto_expand)

    : m_Parameters(params),
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
    if (m_Parameters.m_MustBeWellBehaved && !isWellBehaved()) {
      throw error("SpecialPosition: MinMateDistance too large.");
    }
    if (auto_expand) expand();
  }

  SymEquivCoordinates::SymEquivCoordinates(
    const SpecialPositionSnapParameters& params,
    const uctbx::Vec3& X)
  {
    SpecialPosition SP(params, X, true);
    for(std::size_t i=0;i<SP.M();i++) {
      m_Coordinates.push_back(SP[i] * SP.SnapPosition());
    }
  }

  SymEquivCoordinates::SymEquivCoordinates(const SpecialPosition& SP)
  {
    SP.CheckExpanded();
    for(std::size_t i=0;i<SP.M();i++) {
      m_Coordinates.push_back(SP[i] * SP.SnapPosition());
    }
  }

  SymEquivCoordinates::SymEquivCoordinates(const WyckoffMapping& WM,
                                           const uctbx::Vec3& X)
  {
    const WyckoffPosition& WP = WM.WP();
    WP.CheckExpanded();
    uctbx::Vec3 Xr = WM.snap_to_representative(X);
    for(std::size_t i=0;i<WP.M();i++) {
      m_Coordinates.push_back(WP[i] * Xr);
    }
  }

  SymEquivCoordinates::SymEquivCoordinates(const WyckoffPosition& WP,
                                           const uctbx::Vec3& X)
  {
    WP.CheckExpanded();
    for(std::size_t i=0;i<WP.M();i++) {
      m_Coordinates.push_back(WP[i] * X);
    }
  }

  SymEquivCoordinates::SymEquivCoordinates(
    const SpecialPositionTolerances& params,
    const uctbx::Vec3& X)
  {
    m_Coordinates.push_back(X);
    for(int i=1;i<params.m_SgOps.OrderZ();i++) {
      uctbx::Vec3 SX = params.m_SgOps(i) * X;
      double Delta2 = getShortestDistance2(params.m_UnitCell, SX);
      if (Delta2 >= params.m_Tolerance2) {
        if (Delta2 < params.m_MinimumDistance2) {
          throw error(
          "Special position not well defined."
          " Use SpecialPositionSnapParameters.");
        }
        else {
          m_Coordinates.push_back(SX);
        }
      }
    }
    if (params.m_SgOps.OrderZ() % m_Coordinates.size() != 0) {
      throw error(
      "Numerical instability. Use SpecialPositionSnapParameters.");
    }
  }

  SymEquivCoordinates::SymEquivCoordinates(const SgOps& sgo,
                                           const uctbx::Vec3& X)
  {
    m_Coordinates.push_back(X);
    for(int i=1;i<sgo.OrderZ();i++) {
      uctbx::Vec3 SX = sgo(i) * X;
      m_Coordinates.push_back(SX);
    }
  }

  double
  SymEquivCoordinates::getShortestDistance2(const uctbx::UnitCell& uc,
                                            const uctbx::Vec3& Y) const
  {
    double result = detail::ModLength2(uc, Y - m_Coordinates[0]);
    for(std::size_t i=1;i<m_Coordinates.size();i++) {
      double Delta2 = detail::ModLength2(uc, Y - m_Coordinates[i]);
      if (result > Delta2)
          result = Delta2;
    }
    return result;
  }

  std::complex<double>
  SymEquivCoordinates::StructureFactor(const Miller::Index& H) const
  {
    using cctbx::constants::pi;
    std::complex<double> F(0., 0.);
    rangei(M()) {
      double phase = 2. * pi * (H * m_Coordinates[i]);
      F += std::complex<double>(std::cos(phase), std::sin(phase));
    }
    return F;
  }

} // namespace sgtbx
