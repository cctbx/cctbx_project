// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <algorithm>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  class CmpiVect {
    private:
      int m_n;
    public:
      inline CmpiVect(int n) : m_n(n) { }
      bool operator()(const int *a, const int *b) const {
        int i;
        int n0a = 0; for(i=0;i<m_n;i++) if (a[i] == 0) n0a++;
        int n0b = 0; for(i=0;i<m_n;i++) if (b[i] == 0) n0b++;
        if (n0a > n0b) return true;
        if (n0a < n0b) return false;
        for(i=0;i<m_n;i++) {
          if (a[i] != 0 && b[i] == 0) return true;
          if (a[i] == 0 && b[i] != 0) return false;
        }
        for(i=0;i<m_n;i++) {
          if (abs(a[i]) < abs(b[i])) return true;
          if (abs(a[i]) > abs(b[i])) return false;
        }
        for(i=0;i<m_n;i++) {
          if (a[i] > b[i]) return true;
          if (a[i] < b[i]) return false;
        }
        return false;
      }
  };

  TrVec TrOps::TidyT(const TrVec& T) const
  {
    int lcmBF = lcm(m_Vects[0].BF(), T.BF());
    int fBF = lcmBF / m_Vects[0].BF();
    TrVec Tlcm = T.scale(lcmBF / T.BF());
    TrVec Tbest = Tlcm.ModShort();
    for (int i = 1; i < nVects(); i++) {
      TrVec Ttrial = (Tlcm + m_Vects[i].scale(fBF)).ModShort();
      if (CmpiVect(3)(Ttrial.elems, Tbest.elems)) Tbest = Ttrial;
    }
    return Tbest.newBaseFactor(T.BF()).ModPositive();
  }

  class CmpLTr {
    public:
      bool operator()(const TrVec& a, const TrVec& b) {
        rangei(3) {
          if (a[i] < b[i]) return true;
          if (a[i] > b[i]) return false;
        }
        return false;
      }
  };

  class CmpSMx {
    public:
      bool operator()(const RTMx& a, const RTMx& b) {
        RotMxInfo RI_a = a.Rpart().getInfo();
        RotMxInfo RI_b = b.Rpart().getInfo();
        if (abs(RI_a.Rtype()) > abs(RI_b.Rtype())) return true;
        if (abs(RI_a.Rtype()) < abs(RI_b.Rtype())) return false;
        if (RI_a.Rtype() > RI_b.Rtype()) return true;
        if (RI_a.Rtype() < RI_b.Rtype()) return false;
        if (CmpiVect(3)(RI_a.EV().elems, RI_b.EV().elems)) return true;
        if (CmpiVect(3)(RI_b.EV().elems, RI_a.EV().elems)) return false;
        if (RI_a.SenseOfRotation() > RI_b.SenseOfRotation()) return true;
        if (RI_a.SenseOfRotation() < RI_b.SenseOfRotation()) return false;
        if (CmpiVect(3)(a.Tpart().elems, b.Tpart().elems)) return true;
        if (CmpiVect(3)(b.Tpart().elems, a.Tpart().elems)) return false;
        int i;
        for(i=0;i<9;i++) {
          if (a.Rpart()[i] < b.Rpart()[i]) return true;
          if (a.Rpart()[i] > b.Rpart()[i]) return false;
        }
        for(i=0;i<3;i++) {
          if (a.Tpart()[i] < b.Tpart()[i]) return true;
          if (a.Tpart()[i] > b.Tpart()[i]) return false;
        }
        return false;
      }
  };

  void SgOps::makeTidy()
  {
    if (m_isTidy) return;
    if (m_fInv == 2) {
      m_InvT = m_LTr.TidyT(m_InvT);
      for (int i = 1; i < m_nSMx; i++) {
        if (m_SMx[i].Rpart().getRtype() < 0) {
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

} // namespace sgtbx
