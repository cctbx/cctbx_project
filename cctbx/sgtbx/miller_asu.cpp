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

    namespace detail {

      class ReferenceReciprocalSpaceASU_1b
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_1b;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[2]>0 || (H[2]==0 && (H[0]>0 || (H[0]==0 && H[1]>=0))));
          }
          virtual const char* representation() const {
            return "l>0 or (l==0 and (h>0 or (h==0 and k>=0)))";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(-1, -1, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_2_m
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_2_m;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[1]>=0 && (H[2]>0 || (H[2]==0 && H[0]>=0)));
          }
          virtual const char* representation() const {
            return "k>=0 and (l>0 or (l=0 and h>=0))";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(-1, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_mmm
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_mmm;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=0 && H[1]>=0 && H[2]>=0);
          }
          virtual const char* representation() const {
            return "h>=0 and k>=0 and l>=0";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_4_m
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_4_m;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[2]>=0 && ((H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0)));
          }
          virtual const char* representation() const {
            return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_4_mmm
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_4_mmm;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=H[1] && H[1]>=0 && H[2]>=0);
          }
          virtual const char* representation() const {
            return "h>=k and k>=0 and l>=0";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_3b
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_3b;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0 && H[2]>=0);
          }
          virtual const char* representation() const {
            return "(h>=0 and k>0) or (h=0 and k=0 and l>=0)";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, -1);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_3b1m
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_3b1m;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=H[1] && H[1]>=0 && (H[1]>0 || H[2]>=0));
          }
          virtual const char* representation() const {
            return "h>=k and k>=0 and (k>0 or l>=0)";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, -1);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_3bm1
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_3bm1;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=H[1] && H[1]>=0 && (H[0]>H[1] || H[2]>=0));
          }
          virtual const char* representation() const {
            return "h>=k and k>=0 and (h>k or l>=0)";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, -1);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_6_m
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_6_m;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[2]>=0 && ((H[0]>=0 && H[1]>0) || (H[0]==0 && H[1]==0)));
          }
          virtual const char* representation() const {
            return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_6_mmm
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_6_mmm;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=H[1] && H[1]>=0 && H[2]>=0);
          }
          virtual const char* representation() const {
            return "h>=k and k>=0 and l>=0";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_m3b
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_m3b;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[0]>=0 && (   (H[2]>=H[0] && H[1]>H[0])
                                || (H[2]==H[0] && H[1]==H[0])));
          }
          virtual const char* representation() const {
            return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };
      class ReferenceReciprocalSpaceASU_m3bm
        : public ReferenceReciprocalSpaceASU {
        public:
          virtual tables::MatrixGroup::Code LaueGroupCode() const {
            return tables::MatrixGroup::MGC_m3bm;
          }
          virtual bool isInASU(const Miller::Index& H) const {
            return (H[1]>=H[2] && H[2]>=H[0] && H[0]>=0);
          }
          virtual const char* representation() const {
            return "k>=l and l>=h and h>=0";
          }
          virtual const af::int3& getCutParameters() const {
            static const af::int3 result(0, 0, 0);
            return result;
          }
      };

      static const ReferenceReciprocalSpaceASU_1b    ReferenceASU_1b =
                   ReferenceReciprocalSpaceASU_1b();
      static const ReferenceReciprocalSpaceASU_2_m   ReferenceASU_2_m =
                   ReferenceReciprocalSpaceASU_2_m();
      static const ReferenceReciprocalSpaceASU_mmm   ReferenceASU_mmm =
                   ReferenceReciprocalSpaceASU_mmm();
      static const ReferenceReciprocalSpaceASU_4_m   ReferenceASU_4_m =
                   ReferenceReciprocalSpaceASU_4_m();
      static const ReferenceReciprocalSpaceASU_4_mmm ReferenceASU_4_mmm =
                   ReferenceReciprocalSpaceASU_4_mmm();
      static const ReferenceReciprocalSpaceASU_3b    ReferenceASU_3b =
                   ReferenceReciprocalSpaceASU_3b();
      static const ReferenceReciprocalSpaceASU_3b1m  ReferenceASU_3b1m =
                   ReferenceReciprocalSpaceASU_3b1m();
      static const ReferenceReciprocalSpaceASU_3bm1  ReferenceASU_3bm1 =
                   ReferenceReciprocalSpaceASU_3bm1();
      static const ReferenceReciprocalSpaceASU_6_m   ReferenceASU_6_m =
                   ReferenceReciprocalSpaceASU_6_m();
      static const ReferenceReciprocalSpaceASU_6_mmm ReferenceASU_6_mmm =
                   ReferenceReciprocalSpaceASU_6_mmm();
      static const ReferenceReciprocalSpaceASU_m3b   ReferenceASU_m3b =
                   ReferenceReciprocalSpaceASU_m3b();
      static const ReferenceReciprocalSpaceASU_m3bm  ReferenceASU_m3bm =
                   ReferenceReciprocalSpaceASU_m3bm();

      static const ReferenceReciprocalSpaceASU*
      TableReferenceReciprocalSpaceASU[] = {
        &ReferenceASU_1b,
        &ReferenceASU_2_m,
        &ReferenceASU_mmm,
        &ReferenceASU_4_m, &ReferenceASU_4_mmm,
        &ReferenceASU_3b, &ReferenceASU_3b1m, &ReferenceASU_3bm1,
        &ReferenceASU_6_m, &ReferenceASU_6_mmm,
        &ReferenceASU_m3b, &ReferenceASU_m3bm,
        0,
      };

    } // namespace detail

    ReciprocalSpaceASU::ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo)
      : m_CBOp(), m_isReferenceASU(true), m_ReferenceASU()
    {
      m_CBOp = SgInfo.CBOp();
      m_isReferenceASU = m_CBOp.M().Rpart().isUnit();
      using namespace tables::MatrixGroup;
      Code
      MGC = tables::ReferenceSettings::MatrixGroupCodes[SgInfo.SgNumber()];
      Code LG_MGC = MGC.LaueGroupType();
      if (LG_MGC == MGC_3bm) {
       if (   MGC == MGC_312
           || MGC == MGC_31m
           || MGC == MGC_3b1m) LG_MGC = MGC_3b1m;
       else                    LG_MGC = MGC_3bm1;
      }
      for(std::size_t i=0;; i++) {
        m_ReferenceASU = detail::TableReferenceReciprocalSpaceASU[i];
        if (m_ReferenceASU == 0) throw cctbx_internal_error();
        if (m_ReferenceASU->LaueGroupCode() == LG_MGC) break;
      }
    }

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
        m_ASU(ReciprocalSpaceASU(SgInfo))
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
        m_ASU(ReciprocalSpaceASU(SgInfo)),
        m_Qhigh(-1.)
    {
      InitializeLoop(Miller::Index(af::abs(MaxIndex)));
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
              if (!ReferenceH.is000() && !m_SgOps.isSysAbsent(ReferenceH)) {
                return ReferenceH;
              }
            }
            else {
              double Q = m_UnitCell.Q(ReferenceH);
              if (Q != 0 && Q <= m_Qhigh && !m_SgOps.isSysAbsent(ReferenceH)) {
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
                if (!H.is000() && !m_SgOps.isSysAbsent(H)) {
                  return H;
                }
              }
              else {
                double Q = m_UnitCell.Q(H);
                if (Q != 0 && Q <= m_Qhigh && !m_SgOps.isSysAbsent(H)) {
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
      m_next_is_minus_previous = !m_SgOps.isCentric(m_previous);
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
        if (ASU.isInASU(m_HR.FriedelMate())) {
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
          if (iMate) TrialH = TrialH.FriedelMate();
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
