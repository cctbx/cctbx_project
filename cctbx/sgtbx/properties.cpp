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
#include <cctbx/sgtbx/reference.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  int RotMx::getRtype() const
  {
    int deter = det();
    if (deter == -1 || deter == 1) {
      switch (trace()) {
        case -3:                  return -1;
        case -2:                  return -6;
        case -1: if (deter == -1) return -4;
                 else             return  2;
        case  0: if (deter == -1) return -3;
                 else             return  3;
        case  1: if (deter == -1) return -2;
                 else             return  4;
        case  2:                  return  6;
        case  3:                  return  1;
      }
    }
    return 0;
  }

  int SenseOfRotation(const RotMx& R, int Rtype, const Vec3& EV)
  {
    // M.B. Boisen, Jr. & G.V. Gibbs
    // Mathematical Crystallography, Revised Edition 1990
    // pp. 348-349, 354-356

    int f = 1; if (Rtype < 0) f = -1;
    int trace = f * R.trace();
    if (trace == 3 || trace == -1) return 0; /* 1-fold or 2-fold */
    if (EV[1] == 0 && EV[2] == 0) {
      if (EV[0] * f * R[7] > 0)
        return 1;
    }
    else {
      if (f * (R[3] * EV[2] - R[6] * EV[1]) > 0)
        return 1;
    }
    return -1;
  }

  Vec3 SolveHomRE2(const Mx33& REMx)
  {
    // REMx must be in row echelon form with Rank 2.

    int IxIndep[1];
    cctbx_assert(iRESetIxIndep(REMx.elems, 2, 3, IxIndep, 1) == 1);
    Vec3 EV;
    rangei(3) EV[i] = 0;
    EV[IxIndep[0]] = 1;
    cctbx_assert(iREBacksubst(REMx.elems, 0, 2, 3, EV.elems, 0) >= 1);
    if (SignHemisphere(EV) < 0) {
      rangei(3) EV[i] *= -1;
    }
    return EV;
  }

  RotMxInfo RotMx::getInfo() const
  {
    RotMxInfo result;
    result.m_Rtype = getRtype();
    if (result.m_Rtype == 0) {
      throw error("Cannot determine type of rotation matrix.");
    }
    RotMx ProperR = *this;
    int ProperOrder = result.m_Rtype;
    if (ProperOrder < 0) {
      ProperOrder *= -1;
      ProperR = -ProperR;
    }
    if (ProperOrder > 1) {
      RotMx RmI = ProperR - RotMx(ProperR.BF());
      if (iRowEchelonFormT(RmI.elems, 3, 3, 0, 0) != 2) {
        throw error("Cannot determine Eigenvector of rotation matrix.");
      }
      result.m_EV = SolveHomRE2(RmI);
      result.m_SenseOfRotation
             = SenseOfRotation(*this, result.m_Rtype, result.m_EV);
    }
    return result;
  }

  int OrderOfRtype(int Rtype) {
    if (Rtype > 0) return  Rtype;
    if (Rtype % 2) return -Rtype * 2;
                   return -Rtype;
  }

  RotMx RotMx::accumulate(int Rtype) const
  {
    if (Rtype == 0) Rtype = getRtype();
    cctbx_assert(Rtype != 0);
    if (Rtype == 1) return *this;
    int Order = OrderOfRtype(Rtype);
    RotMx result(m_BF);
    RotMx MxA(*this), MxB;
    RotMx* RR = &MxA;
    RotMx* RRR = &MxB;
    for (int iO = 1;;) {
      result += (*RR);
      if (++iO == Order)
        break;
      *RRR = (*this) * (*RR);
      RotMx* Swp = RR; RR = RRR; RRR = Swp;
    }
    return result;
  }

  bool SpaceGroup::isChiral() const
  {
    if (isCentric()) return false;
    for (int i = 1; i < m_nSMx; i++) {
      if (m_SMx[i].Rpart().getRtype() < 0) return false;
    }
    return true;
  }

  bool SpaceGroup::isEnantiomorphic() const
  {
    if (isCentric()) return false;
    SpaceGroupType SgType = getSpaceGroupType();
    std::vector<RTMx>
    AddlG = ReferenceSettings::GetNormAddlG(SgType.SgNumber(),
                                            false, true, false);
    if (AddlG.size() == 1) return false;
    cctbx_assert(AddlG.size() == 0);
    return true;
  }

  ChOfBasisOp SpaceGroup::getChangeOfHandOp() const
  {
    if (isCentric()) return ChOfBasisOp(1, STBF);
    SpaceGroupType SgType = getSpaceGroupType();
    std::vector<RTMx>
    AddlG = SgType.getAddlGeneratorsOfEuclideanNormalizer(true, false);
    if (AddlG.size() == 1) return ChOfBasisOp(AddlG[0]);
    cctbx_assert(AddlG.size() == 0);
    return ChOfBasisOp(SgType.CBOp().swap()(RTMx(RotMx(1, -1), STBF)));
  }

} // namespace sgtbx
