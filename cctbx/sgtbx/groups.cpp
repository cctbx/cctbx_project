// $Id$

#include <algorithm>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  bool TrOps::add(const TrVec& NewTr)
  {
    if (! NewTr.isValid()) return false;
    TrVec NTr(NewTr.ModPositive());
    if (std::find(m_Vects.begin(), m_Vects.end(), NTr) != m_Vects.end()) {
      return false;
    }
    m_Vects.push_back(NTr);
    return true;
  }

  bool TrOps::expand(const TrVec& NewTr)
  {
    const TrVec *pNewTr = &NewTr;
    int old_size = m_Vects.size();
    int i = old_size;
    int j = 1;
    for (;;) {
      add(*pNewTr);
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == m_Vects.size()) break;
      TrVec TrialTr = m_Vects[j] + m_Vects[j];
      pNewTr = &TrialTr;
      j++;
    }
    return m_Vects.size() != old_size;
  }

  void SgOps::reset()
  {
    nLSL = 1;
    nSSL = 1;
    m_fInv = 1;
    m_nSMx = 1;
    m_LTr.reset();
    m_InvT = TrVec();
    rangei(m_SMx.max_size()) m_SMx[i] = RTMx(RotMx(0, 0), TrVec(0));
    m_SMx[0] = RTMx();
    m_isTidy = true;
  }

  void SgOps::addInv(const TrVec& NewInvT)
  {
    if (m_fInv == 2) { // there is a centre of inversion already
      if (m_LTr.add(m_InvT - NewInvT)) m_isTidy = false;
      return;
    }
    m_InvT = NewInvT.ModPositive();
    m_fInv = 2;
    if (! m_NoExpand) {
      for (int i = 1; i < m_nSMx; i++) {
        if (m_LTr.add(      m_SMx[i].Rpart() * m_InvT
                      + 2 * m_SMx[i].Tpart() - m_InvT)) {
          m_isTidy = false;
        }
      }
    }
  }

  void SgOps::addSMx(const RTMx& NewSMx)
  {
    RotMx mR = -NewSMx.Rpart();
    for (int iSMx = 0; iSMx < m_nSMx; iSMx++) {
      if (m_SMx[iSMx].Rpart() == NewSMx.Rpart()) {
        if (m_LTr.add(m_SMx[iSMx].Tpart() - NewSMx.Tpart())) m_isTidy = false;
        return;
      }
      if (m_SMx[iSMx].Rpart() == mR) {
        addInv(m_SMx[iSMx].Tpart() + NewSMx.Tpart());
        return;
      }
    }
    int d = NewSMx.Rpart().det();
    if (m_nSMx >= m_SMx.max_size() || (d != -1 && d != 1))
      throw error("Non-crystallographic rotation matrix encountered.");
    m_SMx[m_nSMx] = NewSMx;
    m_SMx[m_nSMx].normalize();
    m_nSMx++;
    if (! m_NoExpand && m_fInv == 2) {
      m_LTr.add(      m_SMx[m_nSMx - 1].Rpart() * m_InvT
                + 2 * m_SMx[m_nSMx - 1].Tpart() - m_InvT);
    }
    m_isTidy = false;
  }

  void SgOps::expandLTr(const TrVec& NewLTr)
  {
    if (m_NoExpand) {
      if (m_LTr.add(NewLTr)) m_isTidy = false;
      return;
    }
    for (int iSMx = nSSL; iSMx < m_nSMx; iSMx++) {
      for (int iLTr = 1; iLTr < nLSL; iLTr++) {
        if (m_LTr.add(m_SMx[iSMx].Rpart() * m_LTr[iLTr])) m_isTidy = false;
      }
    }
    nSSL = m_nSMx;
    TrVec TrialLTr = NewLTr;
    int i = nLSL;
    int j = 1;
    for (;;)
    {
      if (m_LTr.add(TrialLTr)) m_isTidy = false;
      for (int iSMx = 1; iSMx < m_nSMx; iSMx++) {
        for (int iLTr = nLSL; iLTr < m_LTr.nVects(); iLTr++) {
          if (m_LTr.add(m_SMx[iSMx].Rpart() * m_LTr[iLTr])) m_isTidy = false;
        }
      }
      nLSL = m_LTr.nVects();
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == m_LTr.nVects()) break;
      TrialLTr = m_LTr[j] + m_LTr[i];
      j++;
    }
  }

  void SgOps::expandInv(const TrVec& NewInvT)
  {
    addInv(NewInvT);
    expandLTr(TrVec(0));
  }

  void SgOps::expandSMx(const RTMx& NewSMx)
  {
    if (m_NoExpand) {
      addSMx(NewSMx);
      return;
    }
    const RTMx *pNewSMx = &NewSMx;
    int i = m_nSMx;
    int j = 1;
    for (;;) {
      addSMx(*pNewSMx);
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == m_nSMx) break;
      RTMx TrialSMx = m_SMx[j] * m_SMx[i];
      pNewSMx = &TrialSMx;
      j++;
    }
    expandLTr(TrVec(0));
  }

  int SgOps::expandConventionalCentringType(char Symbol)
  {
    const lattice::CentringTypeMap*
    mapping = lattice::getConventionalCentringType(Symbol);
    if (mapping == 0) {
      throw error("Illegal symbol for centring type of cell.");
    }
    rangei(mapping->nTrs) {
      expandLTr(mapping->Trs[i]);
    }
    return mapping->nTrs;
  }

  RTMx SgOps::operator()(int iLTr, int iInv, int iSMx) const
  {
    if (   iLTr < 0 || iLTr >= m_LTr.nVects()
        || iInv < 0 || iInv >= m_fInv
        || iSMx < 0 || iSMx >= m_nSMx) {
      throw error("Index out of range.");
    }
    if (iInv == 0) return m_SMx[iSMx] + m_LTr[iLTr];
    return -m_SMx[iSMx] + m_InvT + m_LTr[iLTr];
  }

  RTMx SgOps::operator()(int iLIS) const
  {
    // iLIS = ((iLTr * fInv) + iInv) * nSMx + iSMx
    if (iLIS < 0 || iLIS >= OrderZ()) {
      throw error("Index out of range.");
    }
    int iSMx = iLIS % m_nSMx;
    int iInv = (iLIS / m_nSMx) % m_fInv;
    int iLTr = iLIS / (m_fInv * m_nSMx);
    return operator()(iLTr, iInv, iSMx);
  }

  bool operator==(const SgOps& lhs, const SgOps& rhs)
  {
    if (lhs.nLTr() != rhs.nLTr()) return false;
    if (lhs.fInv() != rhs.fInv()) return false;
    if (lhs.nSMx() != rhs.nSMx()) return false;
    SgOps tidy_lhs = lhs; tidy_lhs.makeTidy();
    SgOps tidy_rhs = rhs; tidy_rhs.makeTidy();
    if (tidy_lhs.m_InvT != tidy_rhs.m_InvT) return false;
    if (tidy_lhs.m_LTr.Vects() != tidy_rhs.m_LTr.Vects()) return false;
    rangei(tidy_lhs.nSMx()) {
      if (tidy_lhs[i] != tidy_rhs[i]) return false;
    }
    return true;
  }

} // namespace sgtbx
