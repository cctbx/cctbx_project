// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 31: Redesign: AsymIndex (rwgk)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace cctbx { namespace sgtbx {

  bool SpaceGroup::isSysAbsent(const Miller::Index& H) const
  {
    // systematically absent reflection: if HR == H and HT != 0 mod 1
    // restricted phase: if HR == -H: phi(H) = pi*HT + n*pi

    int HT_Restriction = -1; // no restriction

    int iSMx;
    for(iSMx=0;iSMx<m_nSMx;iSMx++)
    {
      const RotMx& R = m_SMx[iSMx].Rpart();
      const TrVec& T = m_SMx[iSMx].Tpart();
      TrVec TS(0);
      TrVec TR(0);
      Miller::Index HR = H * R;
      if      (H == HR) {
        TS = T;
        if (isCentric()) TR = m_InvT - T;
      }
      else if (H == HR.FriedelMate()) {
        TR = T;
        if (isCentric()) TS = m_InvT - T;
      }
      if (TS.isValid()) {
        int iLTr;
        for(iLTr=0;iLTr<nLTr();iLTr++) {
          if ((H * (TS + m_LTr[iLTr])) % TS.BF() != 0) {
            return true;
          }
        }
      }
      if (TR.isValid()) {
        int iLTr;
        for(iLTr=0;iLTr<nLTr();iLTr++) {
          int HT = HT_mod_1(H, TR + m_LTr[iLTr]);
          if      (HT_Restriction < 0) HT_Restriction = HT;
          else if (HT_Restriction != HT) {
            return true;
          }
        }
      }
    }
    return false;
  }

  af::shared<bool>
  SpaceGroup::isSysAbsent(const af::shared<Miller::Index>& H) const
  {
    af::shared<bool> result(H.size()); // FUTURE: avoid initialization
    for(std::size_t i=0;i<H.size();i++) {
      result[i] = isSysAbsent(H[i]);
    }
    return result;
  }

  bool SpaceGroup::isCentric(const Miller::Index& H) const
  {
    if (isCentric()) return true;
    rangei(m_nSMx) {
      if (H * m_SMx[i].Rpart() == H.FriedelMate()) return true;
    }
    return false;
  }

  af::shared<bool>
  SpaceGroup::isCentric(const af::shared<Miller::Index>& H) const
  {
    af::shared<bool> result(H.size()); // FUTURE: avoid initialization
    for(std::size_t i=0;i<H.size();i++) {
      result[i] = isCentric(H[i]);
    }
    return result;
  }

  PhaseRestriction
  SpaceGroup::getPhaseRestriction(const Miller::Index& H) const
  {
    // restricted phase: if HR == -H: phi(H) = pi*HT + n*pi
    if (isCentric()) {
      return PhaseRestriction(HT_mod_1(H, m_InvT), TBF());
    }
    rangei(m_nSMx) {
      if (H * m_SMx[i].Rpart() == H.FriedelMate()) {
        return PhaseRestriction(HT_mod_1(H, m_SMx[i].Tpart()), TBF());
      }
    }
    return PhaseRestriction(-1, TBF()); // no restriction
  }

  int SpaceGroup::epsilon(const Miller::Index& H) const
  {
    int result = 0;
    rangei(m_nSMx) {
      Miller::Index HR = H * m_SMx[i].Rpart();
      if (HR == H || (isCentric() && HR == H.FriedelMate()))
        result++;
    }
    cctbx_assert(result != 0 && m_nSMx % result == 0);
    return result;
  }

  af::shared<int>
  SpaceGroup::epsilon(const af::shared<Miller::Index>& H) const
  {
    af::shared<int> result(H.size()); // FUTURE: avoid initialization
    for(std::size_t i=0;i<H.size();i++) {
      result[i] = epsilon(H[i]);
    }
    return result;
  }

  int SpaceGroup::multiplicity(const Miller::Index& H, bool FriedelFlag) const
  {
    if (H.is000()) return 1;
    int Centro = (isCentric() || FriedelFlag);
    int M = 0;
    int R = 0;
    rangei(m_nSMx) {
      Miller::Index HR = H * m_SMx[i].Rpart();
      if      (HR == H) M++;
      else if (HR == H.FriedelMate()) R++;
    }
    cctbx_assert(M != 0 && m_nSMx % M == 0 && (R == 0 || R == M));
    M = m_nSMx / M;
    if (Centro && R == 0) M *= 2;
    return M;
  }

  af::shared<int>
  SpaceGroup::multiplicity(const af::shared<Miller::Index>& H,
                           bool FriedelFlag) const
  {
    af::shared<int> result(H.size()); // FUTURE: avoid initialization
    for(std::size_t i=0;i<H.size();i++) {
      result[i] = multiplicity(H[i], FriedelFlag);
    }
    return result;
  }

  bool PhaseRestriction::isValidPhase(double Period,
                                      double phi, double tolerance) const
  {
    if (m_HT < 0) return true;
    double delta = std::fmod(phi - HT(Period), Period);
    if      (delta >  tolerance) delta -= Period;
    else if (delta < -tolerance) delta += Period;
    if (delta <= tolerance) return true;
    return false;
  }

  SymEquivMillerIndices
  SpaceGroup::getEquivMillerIndices(const Miller::Index& H) const
  {
    SymEquivMillerIndices SEMI(TBF(), OrderP());
    int iInv;
    for(iInv=0;iInv<m_fInv;iInv++) {
      int iSMx;
      for(iSMx=0;iSMx<m_nSMx;iSMx++) {
        RTMx M = operator()(0, iInv, iSMx);
        Miller::Index HR = H * M.Rpart();
        bool found = false;
        for (int i = 0; i < SEMI.N(); i++) {
          if (SEMI[i].HR() == HR) {
            found = true;
            break;
          }
        }
        if (!found) {
          SEMI.add(Miller::SymEquivIndex(HR, HT_mod_1(H, M.Tpart()), TBF()));
        }
      }
    }
    cctbx_assert((m_nSMx * m_fInv) % SEMI.N() == 0);
    return SEMI;
  }

  void SymEquivMillerIndices::add(const Miller::SymEquivIndex& SEI)
  {
    m_List.push_back(SEI);
    if (m_List.size() > 1) {
      if (SEI.HR() == m_List[0].HR().FriedelMate()) {
        cctbx_assert(m_HT_Restriction < 0 || m_HT_Restriction == SEI.HT());
        m_HT_Restriction = SEI.HT();
      }
    }
  }

  Miller::Index SymEquivMillerIndices::operator()(int iMate, int iList) const
  {
    if (   iMate < 0 || iMate >= fMates(true)
        || iList < 0 || iList >= N()) {
      throw error_index();
    }
    if (iMate == 0) return m_List[iList].HR();
    return m_List[iList].HR().FriedelMate();
  }

  Miller::Index SymEquivMillerIndices::operator()(int iIL) const
  {
    // iIL = iMate * N + iList
    if (iIL < 0 || iIL >= M(true)) {
      throw error_index();
    }
    int iList = iIL % N();
    int iMate = iIL / N();
    return operator()(iMate, iList);
  }

}} // namespace cctbx::sgtbx
