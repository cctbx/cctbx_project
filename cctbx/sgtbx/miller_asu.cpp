/* CCP4_ReciprocalSpaceASU is by Kevin Cowtan and is placed in the
   public domain to facilitate compatibility and
   interoperability. This code may be licensed under the cctbx
   license, see files COPYRIGHT.txt and cctbx/LICENSE.txt for further
   details.
*/


#include <cctbx/sgtbx/miller_asu.h>


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
