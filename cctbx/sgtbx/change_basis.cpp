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

  TrVec ChOfBasisOp::operator()(const TrVec& T, int SignI) const
  {
    // (C|V)( I|T)(C^-1|W)=( C|CT+V)(C^-1|W)=( I|CT+V+CW)=( I|C(T+W)+V)
    // (C|V)(-I|T)(C^-1|W)=(-C|CT+V)(C^-1|W)=(-I|CT+V-CW)=(-I|C(T-W)+V)
    TrVec Tf = T.newBaseFactor(InvMx.Tpart().BF());
    TrVec TW;
    if (SignI >= 0) TW = Tf + InvMx.Tpart();
    else            TW = Tf - InvMx.Tpart();
    return (  Mx.Rpart() * TW
            + Mx.Tpart().scale(Mx.Rpart().BF())).newBaseFactor(T.BF());
  }

  RTMx ChOfBasisOp::operator()(const RTMx& M) const
  {
    // This function only works for Seitz matrices:
    cctbx_assert(M.Rpart().BF() == 1);
    cctbx_assert(Mx.Tpart().BF() % M.Tpart().BF() == 0);
    RTMx
    mx_rt_invmx = Mx * (M.scale(1, Mx.Tpart().BF() / M.Tpart().BF()) * InvMx);
    return mx_rt_invmx.newBaseFactors(M);
  }

  MatrixLite::dtype::Vec3
  ChOfBasisOp::operator()(const MatrixLite::dtype::Vec3& X) const {
    MatrixLite::dtype::Vec3 result;
    const double RBF = Mx.RBF();
    const double TBF = Mx.TBF();
    for(int i=0;i<3;i++) {
      result[i] = (  Mx.Rpart()[3 * i + 0] * X[0]
                   + Mx.Rpart()[3 * i + 1] * X[1]
                   + Mx.Rpart()[3 * i + 2] * X[2]) / RBF
                  + Mx.Tpart()[i] / TBF;
    }
    return result;
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

  SgOps SgOps::ChangeBasis(const ChOfBasisOp& CBOp) const
  {
    SgOps result(m_NoExpand);
    result.m_LTr = m_LTr.ChangeBasis(CBOp);
    if (isCentric()) {
      result.expandInv(CBOp(m_InvT, -1));
    }
    for (int i = 1; i < m_nSMx; i++) {
      result.expandSMx(CBOp(m_SMx[i]));
    }
    return result;
  }

} // namespace sgtbx
