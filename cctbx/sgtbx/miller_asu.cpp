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


  } // namespace sgtbx

  namespace miller {

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
            m_HT = sgtbx::HT_mod_1(H, M.Tpart());
            return;
          }
        }
      }
      cctbx_assert(!SgOps.isCentric());
      for(int iSMx=0;iSMx<SgOps.nSMx();iSMx++) {
        sgtbx::RTMx M = SgOps(0, 0, iSMx);
        m_HR = H * M.Rpart();
        if (ASU.isInASU(-m_HR)) {
          m_HT = sgtbx::HT_mod_1(H, M.Tpart());
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

  } // namespace miller
} // namespace cctbx
