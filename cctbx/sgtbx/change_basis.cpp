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

#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace cctbx { namespace sgtbx {

  TrVec ChOfBasisOp::operator()(const TrVec& T, int SignI) const
  {
    // (C|V)( I|T)(C^-1|W)=( C|CT+V)(C^-1|W)=( I|CT+V+CW)=( I|C(T+W)+V)
    // (C|V)(-I|T)(C^-1|W)=(-C|CT+V)(C^-1|W)=(-I|CT+V-CW)=(-I|C(T-W)+V)
    TrVec Tf = T.newBaseFactor(m_InvM.Tpart().BF());
    TrVec TW;
    if (SignI >= 0) TW = Tf + m_InvM.Tpart();
    else            TW = Tf - m_InvM.Tpart();
    return (  m_M.Rpart() * TW
            + m_M.Tpart().scale(m_M.Rpart().BF())).newBaseFactor(T.BF());
  }

  RotMx ChOfBasisOp::operator()(const RotMx& R) const
  {
    cctbx_assert(R.BF() == 1);
    return (m_M.Rpart() * R * m_InvM.Rpart()).newBaseFactor(1);
  }

  RTMx ChOfBasisOp::operator()(const RTMx& RT) const
  {
    // This function only works for Seitz matrices:
    cctbx_assert(RT.RBF() == 1);
    cctbx_assert(m_M.TBF() % RT.TBF() == 0);
    return
      (m_M * (RT.scale(1, m_M.TBF()/RT.TBF()) * m_InvM)).newBaseFactors(RT);
  }

  RTMx ChOfBasisOp::apply(const RTMx& RT) const {
    return m_M.multiply(RT.multiply(m_InvM));
  }

  TrOps TrOps::ChangeBasis(const ChOfBasisOp& CBOp) const
  {
    TrOps result;
    int i;
    for(i=0;i<3;i++) {
      TrVec BV;
      BV[i] = BV.BF();
      result.expand(CBOp(BV, 1));
    }
    for (i = 1; i < nVects(); i++) {
      result.expand(CBOp(m_Vects[i], 1));
    }
    return result;
  }

  SpaceGroup SpaceGroup::ChangeBasis(const ChOfBasisOp& CBOp) const
  {
    SpaceGroup result(m_NoExpand);
    result.m_LTr = m_LTr.ChangeBasis(CBOp);
    if (isCentric()) {
      result.expandInv(CBOp(m_InvT, -1));
    }
    for (int i = 1; i < m_nSMx; i++) {
      result.expandSMx(CBOp(m_SMx[i]));
    }
    return result;
  }

}} // namespace cctbx::sgtbx
