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
#include <cctbx/sgtbx/miller_ref_asu.h>
#include <cctbx/basic/define_range.h>

namespace cctbx { namespace sgtbx {

  PhaseInfo::PhaseInfo(
    SpaceGroup const& sgops,
    Miller::Index const& h,
    bool no_test_sys_absent)
    : m_HT(-1), m_TBF(sgops.TBF()), m_SysAbsWasTested(!no_test_sys_absent)
  {
    // systematically absent reflection: if HR == H and HT != 0 mod 1
    // restricted phase: if HR == -H: phi(H) = pi*HT + n*pi
    if (no_test_sys_absent) {
      // Fast determination of phase restriction without considering
      // conditions for systematically absent reflections.
      if (sgops.isCentric()) {
        m_HT = HT_mod_1(h, sgops.InvT());
      }
      else {
        for(int i=0;i<sgops.nSMx();i++) {
          if (h * sgops[i].Rpart() == -h) {
            m_HT = HT_mod_1(h, sgops[i].Tpart());
            break;
          }
        }
      }
      return;
    }
    // Simulatenous determination of phase restriction and evaluation
    // conditions for systematically absent reflections.
    for(int iSMx=0;iSMx<sgops.nSMx();iSMx++) {
      const RotMx& R = sgops[iSMx].Rpart();
      const TrVec& T = sgops[iSMx].Tpart();
      TrVec TS(0);
      TrVec TR(0);
      Miller::Index HR = h * R;
      if      (h == HR) {
        TS = T;
        if (sgops.isCentric()) TR = sgops.InvT() - T;
      }
      else if (h == -HR) {
        TR = T;
        if (sgops.isCentric()) TS = sgops.InvT() - T;
      }
      if (TS.isValid()) {
        for(int iLTr=0;iLTr<sgops.nLTr();iLTr++) {
          if ((h * (TS + sgops.LTr(iLTr))) % TS.BF() != 0) {
            m_HT = -2;
            return;
          }
        }
      }
      if (TR.isValid()) {
        for(int iLTr=0;iLTr<sgops.nLTr();iLTr++) {
          int HT = HT_mod_1(h, TR + sgops.LTr(iLTr));
          if      (m_HT < 0) m_HT = HT;
          else if (m_HT != HT) {
            m_HT = -2;
            return;
          }
        }
      }
    }
  }

  af::shared<bool>
  SpaceGroup::isSysAbsent(af::shared<Miller::Index> H) const
  {
    af::shared<bool> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(isSysAbsent(H[i]));
    }
    return result;
  }

  bool SpaceGroup::isCentric(const Miller::Index& H) const
  {
    if (isCentric()) return true;
    rangei(m_nSMx) {
      if (H * m_SMx[i].Rpart() == -H) return true;
    }
    return false;
  }

  af::shared<bool>
  SpaceGroup::isCentric(const af::shared<Miller::Index>& H) const
  {
    af::shared<bool> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(isCentric(H[i]));
    }
    return result;
  }

  int SpaceGroup::epsilon(const Miller::Index& H) const
  {
    int result = 0;
    rangei(m_nSMx) {
      Miller::Index HR = H * m_SMx[i].Rpart();
      if (HR == H || (isCentric() && HR == -H))
        result++;
    }
    cctbx_assert(result != 0 && m_nSMx % result == 0);
    return result;
  }

  af::shared<int>
  SpaceGroup::epsilon(const af::shared<Miller::Index>& H) const
  {
    af::shared<int> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(epsilon(H[i]));
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
      else if (HR == -H) R++;
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
    af::shared<int> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(multiplicity(H[i], FriedelFlag));
    }
    return result;
  }

  bool
  PhaseInfo::isValidPhase_(double Period, double phi, double tolerance) const
  {
    if (m_HT < 0) return true;
    double delta = std::fmod(phi - HT_angle(Period), Period);
    if (delta >  tolerance) delta -= Period;
    if (delta < -tolerance) delta += Period;
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
          SEMI.add(
            Miller::SymEquivIndex(HR, HT_mod_1(H, M.Tpart()), TBF(), false));
        }
      }
    }
    cctbx_assert((m_nSMx * m_fInv) % SEMI.N() == 0);
    cctbx_assert(!SEMI.isCentric() || SEMI.N() % 2 == 0);
    SEMI.sort_in_hemispheres();
    return SEMI;
  }

  void SymEquivMillerIndices::add(const Miller::SymEquivIndex& SEI)
  {
    m_List.push_back(SEI);
    if (m_List.size() > 1) {
      if (SEI.HR() == -m_List[0].HR()) {
        cctbx_assert(m_HT_Restriction < 0 || m_HT_Restriction == SEI.HT());
        m_HT_Restriction = SEI.HT();
      }
    }
  }

  void SymEquivMillerIndices::sort_in_hemispheres()
  {
    if (!isCentric()) return;
    std::vector<Miller::SymEquivIndex> plus;
    std::vector<Miller::SymEquivIndex> minus;
    for(std::size_t i=0;i<m_List.size();i++) {
      if (isInReferenceReciprocalSpaceASU_1b(m_List[i].HR())) {
        plus.push_back(m_List[i]);
      }
      else {
        minus.push_back(m_List[i]);
      }
    }
    m_List.clear();
    m_List.insert(m_List.end(), plus.begin(), plus.end());
    m_List.insert(m_List.end(), minus.begin(), minus.end());
  }

  Miller::SymEquivIndex
  SymEquivMillerIndices::operator()(int iMate, int iList) const
  {
    if (   iMate < 0 || iMate >= fMates(true)
        || iList < 0 || iList >= N()) {
      throw error_index();
    }
    return m_List[iList].Mate(iMate);
  }

  SymEquivMillerIndices::iIL_decomposition
  SymEquivMillerIndices::decompose_iIL(int iIL) const
  {
    // iIL = iMate * N + iList
    if (iIL < 0 || iIL >= M(true)) {
      throw error_index();
    }
    return iIL_decomposition(iIL / N(), iIL % N());
  }

  Miller::SymEquivIndex
  SymEquivMillerIndices::operator()(int iIL) const
  {
    iIL_decomposition d = decompose_iIL(iIL);
    return operator()(d.iMate, d.iList);
  }

}} // namespace cctbx::sgtbx
