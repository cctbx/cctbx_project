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

    class StdReciprocalSpaceASU_1b : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_2_m : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_mmm : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_4_m : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_4_mmm : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_3b : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_3b1m : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_3bm1 : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_6_m : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_6_mmm : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_m3b : public StdReciprocalSpaceASU {
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
    class StdReciprocalSpaceASU_m3bm: public StdReciprocalSpaceASU {
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

    static const StdReciprocalSpaceASU_1b    StdASU_1b =
                 StdReciprocalSpaceASU_1b();
    static const StdReciprocalSpaceASU_2_m   StdASU_2_m =
                 StdReciprocalSpaceASU_2_m();
    static const StdReciprocalSpaceASU_mmm   StdASU_mmm =
                 StdReciprocalSpaceASU_mmm();
    static const StdReciprocalSpaceASU_4_m   StdASU_4_m =
                 StdReciprocalSpaceASU_4_m();
    static const StdReciprocalSpaceASU_4_mmm StdASU_4_mmm =
                 StdReciprocalSpaceASU_4_mmm();
    static const StdReciprocalSpaceASU_3b    StdASU_3b =
                 StdReciprocalSpaceASU_3b();
    static const StdReciprocalSpaceASU_3b1m  StdASU_3b1m =
                 StdReciprocalSpaceASU_3b1m();
    static const StdReciprocalSpaceASU_3bm1  StdASU_3bm1 =
                 StdReciprocalSpaceASU_3bm1();
    static const StdReciprocalSpaceASU_6_m   StdASU_6_m =
                 StdReciprocalSpaceASU_6_m();
    static const StdReciprocalSpaceASU_6_mmm StdASU_6_mmm =
                 StdReciprocalSpaceASU_6_mmm();
    static const StdReciprocalSpaceASU_m3b   StdASU_m3b =
                 StdReciprocalSpaceASU_m3b();
    static const StdReciprocalSpaceASU_m3bm  StdASU_m3bm =
                 StdReciprocalSpaceASU_m3bm();

    static const StdReciprocalSpaceASU* TableStdReciprocalSpaceASU[] = {
      &StdASU_1b,
      &StdASU_2_m,
      &StdASU_mmm,
      &StdASU_4_m, &StdASU_4_mmm,
      &StdASU_3b, &StdASU_3b1m, &StdASU_3bm1,
      &StdASU_6_m, &StdASU_6_mmm,
      &StdASU_m3b, &StdASU_m3bm,
      0,
    };

  } // namespace detail

  ReciprocalSpaceASU::ReciprocalSpaceASU(const SpaceGroupType& SgType)
    : m_CBOp(), m_isStdASU(true), m_StdASU()
  {
    m_CBOp = SgType.CBOp();
    m_isStdASU = m_CBOp.M().isUnit();
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
      m_StdASU = detail::TableStdReciprocalSpaceASU[i];
      if (m_StdASU == 0) throw cctbx_internal_error();
      if (m_StdASU->LaueGroupCode() == LG_MGC) break;
    }
  }

  class MillerIndexGenerator
  {
  };

} // namespace sgtbx
