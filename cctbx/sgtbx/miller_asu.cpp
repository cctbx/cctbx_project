/* CCP4_ReciprocalSpaceASU is by Kevin Cowtan and is placed in the
   public domain to facilitate compatibility and
   interoperability. This code may be licensed under the cctbx
   license, see files COPYRIGHT.txt and cctbx/LICENSE.txt for further
   details.
*/


#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/sgtbx/reference.h>


namespace sgtbx
{

  void CCP4_ReciprocalSpaceASU::init(const SgOps& sgops)
  {
    // Calculate the reduced (Patterson) spacegroup
    SgOps sg_reduced("P -1");                     // create as -P1
    // now construct reduced sg using only primitive non-translation ops
    for ( int s = 0; s < sgops.nSMx(); s++ )
      sg_reduced.expandSMx( RTMx( sgops(0,0,s).Rpart() ) );
    SpaceGroupType type = sg_reduced.getSpaceGroupType(true);

    // Get the spacegroup number
    rsg = type.SgNumber();
    // Get the change-of-basis to get to a tabulated setting
    cbm = type.CBOp().InvM().Rpart();

    // now pick the appropriate ASU function
    if ( cbm.isUnit() ) {
      // Spacegroup can reduce to one of 14 standard settings
      if      ( rsg == 2   ) asufn = &CCP4_ReciprocalSpaceASU::ASU_1b   ;
      else if ( rsg == 10  ) asufn = &CCP4_ReciprocalSpaceASU::ASU_2_m  ;
      else if ( rsg == 47  ) asufn = &CCP4_ReciprocalSpaceASU::ASU_mmm  ;
      else if ( rsg == 83  ) asufn = &CCP4_ReciprocalSpaceASU::ASU_4_m  ;
      else if ( rsg == 123 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_4_mmm;
      else if ( rsg == 147 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_3b   ;
      else if ( rsg == 148 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_3b   ;
      else if ( rsg == 162 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_3bm  ;
      else if ( rsg == 164 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_3bmx ;
      else if ( rsg == 166 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_3bmx ;
      else if ( rsg == 175 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_6_m  ;
      else if ( rsg == 191 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_6_mmm;
      else if ( rsg == 200 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_m3b  ;
      else if ( rsg == 221 ) asufn = &CCP4_ReciprocalSpaceASU::ASU_m3bm ;
      else throw cctbx_internal_error();
    } else {
      // Non-standard settings involve a change-of-basis op first
      if      ( rsg == 2   ) asufn = &CCP4_ReciprocalSpaceASU::cASU_1b   ;
      else if ( rsg == 10  ) asufn = &CCP4_ReciprocalSpaceASU::cASU_2_m  ;
      else if ( rsg == 47  ) asufn = &CCP4_ReciprocalSpaceASU::cASU_mmm  ;
      else if ( rsg == 83  ) asufn = &CCP4_ReciprocalSpaceASU::cASU_4_m  ;
      else if ( rsg == 123 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_4_mmm;
      else if ( rsg == 147 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_3b   ;
      else if ( rsg == 148 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_3b   ;
      else if ( rsg == 162 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_3bm  ;
      else if ( rsg == 164 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_3bmx ;
      else if ( rsg == 166 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_3bmx ;
      else if ( rsg == 175 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_6_m  ;
      else if ( rsg == 191 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_6_mmm;
      else if ( rsg == 200 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_m3b  ;
      else if ( rsg == 221 ) asufn = &CCP4_ReciprocalSpaceASU::cASU_m3bm ;
      else throw cctbx_internal_error();
    }
  }

  tables::MatrixGroup::Code CCP4_ReciprocalSpaceASU::getLaueGroupType() const
  {
    if      ( rsg == 2   ) return tables::MatrixGroup::MGC_1b;
    else if ( rsg == 10  ) return tables::MatrixGroup::MGC_2_m;
    else if ( rsg == 47  ) return tables::MatrixGroup::MGC_mmm;
    else if ( rsg == 83  ) return tables::MatrixGroup::MGC_4_m;
    else if ( rsg == 123 ) return tables::MatrixGroup::MGC_4_mmm;
    else if ( rsg == 147 ) return tables::MatrixGroup::MGC_3b;
    else if ( rsg == 148 ) return tables::MatrixGroup::MGC_3b;
    else if ( rsg == 162 ) return tables::MatrixGroup::MGC_3bm;
    else if ( rsg == 164 ) return tables::MatrixGroup::MGC_3bm;
    else if ( rsg == 166 ) return tables::MatrixGroup::MGC_3bm;
    else if ( rsg == 175 ) return tables::MatrixGroup::MGC_6_m;
    else if ( rsg == 191 ) return tables::MatrixGroup::MGC_6_mmm;
    else if ( rsg == 200 ) return tables::MatrixGroup::MGC_m3b;
    else if ( rsg == 221 ) return tables::MatrixGroup::MGC_m3bm;
    throw cctbx_internal_error();
  }

  std::string CCP4_ReciprocalSpaceASU::getConditions() const
  {
    if      ( rsg == 2   ) return "l>0 or (l==0 and (h>0 or (h==0 and k>=0)))";
    else if ( rsg == 10  ) return "k>=0 and (l>0 or (l=0 and h>=0))";
    else if ( rsg == 47  ) return "h>=0 and k>=0 and l>=0";
    else if ( rsg == 83  ) return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
    else if ( rsg == 123 ) return "h>=k and k>=0 and l>=0";
    else if ( rsg == 147 ) return "(h>=0 and k>0) or (h=0 and k=0 and l >= 0)";
    else if ( rsg == 148 ) return "(h>=0 and k>0) or (h=0 and k=0 and l >= 0)";
    else if ( rsg == 162 ) return "h>=k and k>=0 and (k>0 or l>=0)";
    else if ( rsg == 164 ) return "h>=k and k>=0 and (h>k or l>=0)";
    else if ( rsg == 166 ) return "h>=k and k>=0 and (h>k or l>=0)";
    else if ( rsg == 175 ) return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
    else if ( rsg == 191 ) return "h>=k and k>=0 and l>=0";
    else if ( rsg == 200 ) return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
    else if ( rsg == 221 ) return "k>=l and l>=h and h>=0";
    throw cctbx_internal_error();
  }

  RotMx CCP4_ReciprocalSpaceASU::getRotToStandardSetting() const
  {
    return cbm;
  }

  Miller::Vec3 CCP4_ReciprocalSpaceASU::getCutParameters() const
  {
    Miller::Vec3 v000 = { 0, 0, 0};
    Miller::Vec3 v001 = { 0, 0,-1};
    Miller::Vec3 v100 = {-1, 0, 0};
    Miller::Vec3 v110 = {-1,-1, 0};
    Miller::Vec3 v111 = {-1,-1,-1};

    if (!cbm.isUnit() ) return v111;
    if      ( rsg == 2   ) return v110;
    else if ( rsg == 10  ) return v100;
    else if ( rsg == 47  ) return v000;
    else if ( rsg == 83  ) return v000;
    else if ( rsg == 123 ) return v000;
    else if ( rsg == 147 ) return v001;
    else if ( rsg == 148 ) return v001;
    else if ( rsg == 162 ) return v001;
    else if ( rsg == 164 ) return v001;
    else if ( rsg == 166 ) return v001;
    else if ( rsg == 175 ) return v000;
    else if ( rsg == 191 ) return v000;
    else if ( rsg == 200 ) return v000;
    else if ( rsg == 221 ) return v000;
    throw cctbx_internal_error();
  }

}

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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {-1, -1, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {-1, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, -1};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, -1};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, -1};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
          return
        (H[0]>=0 && ((H[2]>=H[0] && H[1]>H[0]) || (H[2]==H[0] && H[1]==H[0])));
        }
        virtual const char* representation() const {
          return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
        }
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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
        virtual const Miller::Vec3& getCutParameters() const {
          static const Miller::Vec3 result = {0, 0, 0};
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

  ReciprocalSpaceASU::ReciprocalSpaceASU(const SpaceGroupType& SgType)
    : m_CBOp(), m_isReferenceASU(true), m_ReferenceASU()
  {
    m_CBOp = SgType.CBOp();
    m_isReferenceASU = m_CBOp.M().Rpart().isUnit();
    using namespace tables::MatrixGroup;
    Code MGC = tables::ReferenceSettings::MatrixGroupCodes[SgType.SgNumber()];
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

  MillerIndexGenerator::MillerIndexGenerator(const uctbx::UnitCell& uc,
                                             const SgOps& sgo,
                                             double Resolution_d_min)
    : m_UnitCell(uc), m_SgOps(sgo)
  {
    if (Resolution_d_min <= 0.) {
      throw error("Resolution limit must be greater than zero.");
    }
    m_Qhigh = 1. / (Resolution_d_min * Resolution_d_min);

    SpaceGroupType SgType = m_SgOps.getSpaceGroupType();
    uctbx::UnitCell
    ReferenceUnitCell = m_UnitCell.ChangeBasis(SgType.CBOp().InvM().Rpart());
    SgOps ReferenceSgOps = SgOps(SpaceGroupSymbols(SgType.SgNumber()).Hall());
    m_ASU = ReciprocalSpaceASU(SgType);
    Miller::Vec3 CutP = m_ASU.ReferenceASU()->getCutParameters();
    Miller::Index
    ReferenceHmax = ReferenceUnitCell.MaxMillerIndices(Resolution_d_min);
    Miller::Index ReferenceHmin;
    for(std::size_t i=0;i<3;i++) ReferenceHmin[i] = ReferenceHmax[i] * CutP[i];
    m_loop = NestedLoop<Miller::Index>(ReferenceHmin, ReferenceHmax);
  }

  Miller::Index MillerIndexGenerator::next()
  {
    const int RBF = m_ASU.CBOp().M().RBF();
    for (; m_loop.over() == 0;) {
      Miller::Index ReferenceH = m_loop();
      m_loop.incr();
      if (m_ASU.ReferenceASU()->isInASU(ReferenceH)) {
        if (m_ASU.isReferenceASU()) {
          double Q = m_UnitCell.Q(ReferenceH);
          if (Q != 0 && Q <= m_Qhigh && !m_SgOps.isSysAbsent(ReferenceH)) {
            return ReferenceH;
          }
        }
        else {
          TrVec HR(ReferenceH * m_ASU.CBOp().M().Rpart(), RBF);
          HR = HR.cancel();
          if (HR.BF() == 1) {
            Miller::Index H(HR.elems);
            double Q = m_UnitCell.Q(H);
            if (Q != 0 && Q <= m_Qhigh && !m_SgOps.isSysAbsent(H)) {
              return H;
            }
          }
        }
      }
    }
    return Miller::Index(0, 0, 0);
  }

} // namespace sgtbx

namespace cctbx {
  namespace Miller {

    SymUniqueIndex::SymUniqueIndex(const sgtbx::SgOps& sgo,
                                   const sgtbx::ReciprocalSpaceASU& ASU,
                                   const Index& H)
    {
      m_TBF = sgo.TBF();
      for(int iInv=0;iInv<sgo.fInv();iInv++) {
        for(int iSMx=0;iSMx<sgo.nSMx();iSMx++) {
          sgtbx::RTMx M = sgo(0, iInv, iSMx);
          m_HR = H * M.Rpart();
          if (ASU.isInASU(m_HR)) {
            m_HT = H * M.Tpart();
            m_iMate = 0;
            m_H = m_HR;
            return;
          }
        }
      }
      cctbx_assert(!sgo.isCentric());
      for(int iSMx=0;iSMx<sgo.nSMx();iSMx++) {
        sgtbx::RTMx M = sgo(0, 0, iSMx);
        m_HR = H * M.Rpart();
        Miller::Index HRF = m_HR.FriedelMate();
        if (ASU.isInASU(HRF)) {
          m_HT = H * M.Tpart();
          m_iMate = 1;
          m_H = HRF;
          return;
        }
      }
      throw cctbx_internal_error();
    }

  } // namespace Miller
} // namespace cctbx
