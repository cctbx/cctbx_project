// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <algorithm>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace cctbx { namespace sgtbx {

  TrVec TrOps::TidyT(const TrVec& T) const
  {
    int lcmBF = lcm(m_Vects[0].BF(), T.BF());
    int fBF = lcmBF / m_Vects[0].BF();
    TrVec Tlcm = T.scale(lcmBF / T.BF());
    TrVec Tbest = Tlcm.modShort();
    for (int i = 1; i < nVects(); i++) {
      TrVec Ttrial = (Tlcm + m_Vects[i].scale(fBF)).modShort();
      if (CmpiVect(3)(Ttrial.vec().begin(), Tbest.vec().begin())) {
        Tbest = Ttrial;
      }
    }
    return Tbest.newBaseFactor(T.BF()).modPositive();
  }

  class CmpLTr {
    public:
      bool operator()(const TrVec& a, const TrVec& b) {
        rangei(3) {
          if (a.vec()[i] < b.vec()[i]) return true;
          if (a.vec()[i] > b.vec()[i]) return false;
        }
        return false;
      }
  };

  class CmpSMx {
    public:
      bool operator()(const RTMx& a, const RTMx& b) {
        RotMxInfo RI_a = a.Rpart().getInfo();
        RotMxInfo RI_b = b.Rpart().getInfo();
        if (std::abs(RI_a.Rtype()) > std::abs(RI_b.Rtype())) return true;
        if (std::abs(RI_a.Rtype()) < std::abs(RI_b.Rtype())) return false;
        if (RI_a.Rtype() > RI_b.Rtype()) return true;
        if (RI_a.Rtype() < RI_b.Rtype()) return false;
        if (CmpiVect(3)(RI_a.EV().elems, RI_b.EV().elems)) return true;
        if (CmpiVect(3)(RI_b.EV().elems, RI_a.EV().elems)) return false;
        if (RI_a.SenseOfRotation() > RI_b.SenseOfRotation()) return true;
        if (RI_a.SenseOfRotation() < RI_b.SenseOfRotation()) return false;
        if (CmpiVect(3)(
          a.Tpart().vec().begin(), b.Tpart().vec().begin())) return true;
        if (CmpiVect(3)(
          b.Tpart().vec().begin(), a.Tpart().vec().begin())) return false;
        int i;
        for(i=0;i<9;i++) {
          if (a.Rpart()[i] < b.Rpart()[i]) return true;
          if (a.Rpart()[i] > b.Rpart()[i]) return false;
        }
        for(i=0;i<3;i++) {
          if (a.Tpart().vec()[i] < b.Tpart().vec()[i]) return true;
          if (a.Tpart().vec()[i] > b.Tpart().vec()[i]) return false;
        }
        return false;
      }
  };

  void SpaceGroup::makeTidy()
  {
    if (m_isTidy) return;
    if (isCentric()) {
      m_InvT = InvT(true);
      for (int i = 1; i < m_nSMx; i++) {
        if (m_SMx[i].Rpart().det() < 0) {
          m_SMx[i] = m_SMx[i].pre_multiply_InvT(m_InvT);
        }
      }
    }
    for (int i = 1; i < m_nSMx; i++) {
      m_SMx[i] = RTMx(m_SMx[i].Rpart(), m_LTr.TidyT(m_SMx[i].Tpart()));
    }
    if (nLTr() > 2)
      std::sort(m_LTr.m_Vects.begin() + 1, m_LTr.m_Vects.end(), CmpLTr());
    if (m_nSMx > 2)
      std::sort(m_SMx.begin() + 1, m_SMx.begin() + m_nSMx, CmpSMx());
    m_isTidy = true;
  }

}} // namespace cctbx::sgtbx
