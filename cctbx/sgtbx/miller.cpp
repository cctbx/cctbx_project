// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  bool SgOps::isSysAbsent(const Miller::Index& H) const
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

  bool SgOps::isCentric(const Miller::Index& H) const
  {
    if (isCentric()) return true;
    rangei(m_nSMx) {
      if (H * m_SMx[i].Rpart() == H.FriedelMate()) return true;
    }
    return false;
  }

  PhaseRestriction SgOps::getPhaseRestriction(const Miller::Index& H) const
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

  int SgOps::epsilon(const Miller::Index& H) const
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

  int SgOps::multiplicity(const Miller::Index& H, bool FriedelSym) const
  {
    if (H.is000()) return 1;
    int Centro = (isCentric() || FriedelSym);
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
  SgOps::getEquivMillerIndices(const Miller::Index& H) const
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
      throw error("Index out of range.");
    }
    if (iMate == 0) return m_List[iList].HR();
    return m_List[iList].HR().FriedelMate();
  }

  Miller::Index SymEquivMillerIndices::operator()(int iIL) const
  {
    // iIL = iMate * N + iList
    if (iIL < 0 || iIL >= M(true)) {
      throw error("Index out of range.");
    }
    int iList = iIL % N();
    int iMate = iIL / N();
    return operator()(iMate, iList);
  }

  namespace {

    int OneMxCutPRange(const RotMx& R)
    {
      int m = 0;
      int ic;
      for(ic=0;ic<3;ic++) {
        int s = 0;
        int ir;
        for(ir=0;ir<3;ir++) s += std::abs(R(ir, ic));
        if (m < s)
            m = s;
      }
      return m + 1;
    }

    int SetCheckCutPRange(const SgOps& sgo)
    {
      int Range = 0;
      int iSMx;
      for(iSMx=0;iSMx<sgo.nSMx();iSMx++) {
                int m = OneMxCutPRange(sgo[iSMx].Rpart());
        if (Range < m)
            Range = m;
      }
      return Range;
    }

    int CheckCutParam(const SgOps& sgo, bool FriedelSym, Miller::Vec3 CutP,
                      int Range, int FullBlock)
    {
      Miller::Vec3 AdjRange;
      rangei(3) AdjRange[i] = Range;
      int iBV;
      for(iBV=0;iBV<3;iBV++) {
        Miller::Vec3 Step;
        rangei(3) Step[i] = 1;
        if (FullBlock == 0) Step[iBV] = 2 * Range;
        Miller::Index H;
        for (H[0] = -AdjRange[0]; H[0] <= AdjRange[0]; H[0] += Step[0])
        for (H[1] = -AdjRange[1]; H[1] <= AdjRange[1]; H[1] += Step[1])
        for (H[2] = -AdjRange[2]; H[2] <= AdjRange[2]; H[2] += Step[2])
        {
          // search for equivalent hkl in an active octant
          SymEquivMillerIndices
          SEMI = sgo.getEquivMillerIndices(H);
          int iEq;
          for (iEq = 0; iEq < SEMI.N(); iEq++) {
            int i;
            for (i = 0; i < 3; i++)
              if (CutP[i] == 0 && SEMI[iEq].HR()[i] < 0) break;
            if (i == 3) break;
            if (SEMI.fMates(FriedelSym) != 2) continue;
            for (i = 0; i < 3; i++)
              if (CutP[i] == 0 && SEMI[iEq].HR()[i] > 0) break;
            if (i == 3) break;
          }
          if (iEq == SEMI.N()) return 0; // CutParam does not work
        }
        if (FullBlock != 0) break;
        AdjRange[iBV]--;
      }
      return 1; /* CutParam works fine */
    }

  } // namespace <anonymous>

  Miller::Vec3 SgOps::getCutParameters(bool FriedelSym) const
  {
    const int nTrials = 8;
    static const Miller::Vec3 ListTrialCutP[nTrials] = {
      { 0,  0,  0},
      { 0, -1,  0},
      {-1,  0,  0},
      { 0,  0, -1},
      {-1, -1,  0},
      { 0, -1, -1},
      {-1,  0, -1},
      {-1, -1, -1},
    };
    int Range = SetCheckCutPRange(*this);
    int iTrial;
    for (iTrial = 0; iTrial < nTrials - 1; iTrial++) {
      if (CheckCutParam(*this, FriedelSym, ListTrialCutP[iTrial],
                        Range, 0) != 0) {
        return ListTrialCutP[iTrial];
      }
    }
    return ListTrialCutP[nTrials - 1];
  }

  void SymEquivMillerIndices::setMasterIndex(Miller::MasterIndex& Master) const
  {
    int iList;
    for(iList=0;iList<N();iList++) {
      const Miller::SymEquivIndex& HS = m_List[iList];
      Miller::Index TrialH = HS.HR();
      int iMate;
      for(iMate = 0; iMate < fMates(true); iMate++) {
        if (iMate) TrialH = TrialH.FriedelMate();
        Master.Update(TrialH, iMate, HS);
      }
    }
  }

  Miller::MasterIndex
  SymEquivMillerIndices::getMasterIndex(bool Pretty) const
  {
    Miller::MasterIndex Master(Pretty);
    setMasterIndex(Master);
    cctbx_assert(Master.TBF() != 0);
    return Master;
  }

  Miller::MasterIndex
  SymEquivMillerIndices::getMasterIndex(const Miller::Vec3& CutP,
                                        bool Pretty) const
  {
    Miller::MasterIndex Master(CutP, Pretty);
    setMasterIndex(Master);
    if (Master.TBF() == 0) {
      throw error("No MasterIndex found due to improper CutParameters.");
    }
    return Master;
  }

  Miller::MasterIndex
  SgOps::getMasterIndex(const Miller::Index& H,
                        bool Pretty) const
  {
    SymEquivMillerIndices SEMI = getEquivMillerIndices(H);
    return SEMI.getMasterIndex(Pretty);
  }

  Miller::MasterIndex
  SgOps::getMasterIndex(const Miller::Index& H, const Miller::Vec3& CutP,
                        bool Pretty) const
  {
    SymEquivMillerIndices SEMI = getEquivMillerIndices(H);
    return SEMI.getMasterIndex(CutP, Pretty);
  }

} // namespace sgtbx
