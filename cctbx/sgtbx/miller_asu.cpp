// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 31: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/sgtbx/reference.h>

namespace cctbx {
  namespace sgtbx {

    ReciprocalSpaceASU::ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo)
      : m_CBOp(SgInfo.CBOp()),
        m_isReferenceASU(SgInfo.CBOp().M().Rpart().isUnit()),
        m_ReferenceASU(LookupReferenceReciprocalSpaceASU(
          tables::ReferenceSettings::MatrixGroupCodes[SgInfo.SgNumber()]))
    {}

    void
    MillerIndexGenerator::InitializeLoop(const Miller::Index& ReferenceHmax)
    {
      af::int3 CutP = m_ASU.ReferenceASU()->getCutParameters();
      Miller::Index ReferenceHbegin;
      Miller::Index ReferenceHend;
      for(std::size_t i=0;i<3;i++) {
        ReferenceHbegin[i] = ReferenceHmax[i] * CutP[i];
        ReferenceHend[i] = ReferenceHmax[i] + 1;
      }
      m_loop = af::nested_loop<Miller::Index>(ReferenceHbegin, ReferenceHend);
      m_next_is_minus_previous = false;
    }

    MillerIndexGenerator::MillerIndexGenerator(const uctbx::UnitCell& uc,
                                               const SpaceGroupInfo& SgInfo,
                                               bool FriedelFlag,
                                               double Resolution_d_min)
      : m_UnitCell(uc),
        m_SgOps(SgInfo.SgOps()),
        m_FriedelFlag(FriedelFlag),
        m_ASU(SgInfo)
    {
      if (Resolution_d_min <= 0.) {
        throw error("Resolution limit must be greater than zero.");
      }
      m_Qhigh = 1. / (Resolution_d_min * Resolution_d_min);
      uctbx::UnitCell
      ReferenceUnitCell = m_UnitCell.ChangeBasis(SgInfo.CBOp().InvM().Rpart());
      InitializeLoop(ReferenceUnitCell.MaxMillerIndices(Resolution_d_min));
    }

    MillerIndexGenerator::MillerIndexGenerator(const SpaceGroupInfo& SgInfo,
                                               bool FriedelFlag,
                                               const Miller::Index& MaxIndex)
      : m_UnitCell(),
        m_SgOps(SgInfo.SgOps()),
        m_FriedelFlag(FriedelFlag),
        m_ASU(SgInfo),
        m_Qhigh(-1.)
    {
      InitializeLoop(Miller::Index(af::abs(MaxIndex)));
    }

    bool MillerIndexGenerator::set_sys_abs_test(const Miller::Index& h)
    {
      m_phase_info = PhaseInfo(m_SgOps, h, false);
      return m_phase_info.isSysAbsent();
    }

    Miller::Index MillerIndexGenerator::next_under_friedel_symmetry()
    {
      const int RBF = m_ASU.CBOp().M().RBF();
      for (; m_loop.over() == 0;) {
        Miller::Index ReferenceH = m_loop();
        m_loop.incr();
        if (m_ASU.ReferenceASU()->isInASU(ReferenceH)) {
          if (m_ASU.isReferenceASU()) {
            if (m_Qhigh < 0.) {
              if (!ReferenceH.is000() && !set_sys_abs_test(ReferenceH)) {
                return ReferenceH;
              }
            }
            else {
              double Q = m_UnitCell.Q(ReferenceH);
              if (Q != 0 && Q <= m_Qhigh && !set_sys_abs_test(ReferenceH)) {
                return ReferenceH;
              }
            }
          }
          else {
            TrVec HR(ReferenceH * m_ASU.CBOp().M().Rpart(), RBF);
            HR = HR.cancel();
            if (HR.BF() == 1) {
              Miller::Index H(HR.vec());
              if (m_Qhigh < 0.) {
                if (!H.is000() && !set_sys_abs_test(H)) {
                  return H;
                }
              }
              else {
                double Q = m_UnitCell.Q(H);
                if (Q != 0 && Q <= m_Qhigh && !set_sys_abs_test(H)) {
                  return H;
                }
              }
            }
          }
        }
      }
      return Miller::Index(0, 0, 0);
    }

    Miller::Index MillerIndexGenerator::next()
    {
      if (m_FriedelFlag) return next_under_friedel_symmetry();
      if (m_next_is_minus_previous) {
        m_next_is_minus_previous = false;
        return -m_previous;
      }
      m_previous = next_under_friedel_symmetry();
      m_next_is_minus_previous = !m_phase_info.isCentric();
      return m_previous;
    }

  } // namespace sgtbx

  namespace Miller {

    AsymIndex::AsymIndex(
      const sgtbx::SpaceGroup& SgOps,
      const sgtbx::ReciprocalSpaceASU& ASU,
      const Index& H)
    {
      m_TBF = SgOps.TBF();
      m_FriedelFlag = false;
      for(int iInv=0;iInv<SgOps.fInv();iInv++) {
        for(int iSMx=0;iSMx<SgOps.nSMx();iSMx++) {
          sgtbx::RTMx M = SgOps(0, iInv, iSMx);
          m_HR = H * M.Rpart();
          if (ASU.isInASU(m_HR)) {
            m_HT = H * M.Tpart();
            return;
          }
        }
      }
      cctbx_assert(!SgOps.isCentric());
      for(int iSMx=0;iSMx<SgOps.nSMx();iSMx++) {
        sgtbx::RTMx M = SgOps(0, 0, iSMx);
        m_HR = H * M.Rpart();
        if (ASU.isInASU(-m_HR)) {
          m_HT = H * M.Tpart();
          m_FriedelFlag = true;
          return;
        }
      }
      throw cctbx_internal_error();
    }

    AsymIndex::AsymIndex(
      const sgtbx::SymEquivMillerIndices& SEMI)
    {
      m_TBF = SEMI[0].TBF();
      int iSelected = 0;
      Index SelectedH = SEMI[0].HR();
      m_FriedelFlag = false;
      for(int iList=0;iList<SEMI.N();iList++) {
        const SymEquivIndex& SEI = SEMI[iList];
        Index TrialH = SEI.HR();
        for(int iMate = 0; iMate < SEMI.fMates(true); iMate++) {
          if (iMate) TrialH = -TrialH;
          if (TrialH < SelectedH) {
            iSelected = iList;
            SelectedH = TrialH;
            m_FriedelFlag = (iMate != 0);
          }
        }
      }
      m_HR = SEMI[iSelected].HR();
      m_HT = SEMI[iSelected].HT();
    }

    AsymIndex::AsymIndex(
      const sgtbx::SpaceGroup& SgOps,
      const Index& H) {
      *this = AsymIndex(SgOps.getEquivMillerIndices(H));
    }

  } // namespace Miller
} // namespace cctbx
