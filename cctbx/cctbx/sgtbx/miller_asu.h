/* CCP4_ReciprocalSpaceASU is by Kevin Cowtan and is placed in the
   public domain to facilitate compatibility and
   interoperability. This code may be licensed under the cctbx
   license, see files COPYRIGHT.txt and cctbx/LICENSE.txt for further
   details.
*/

#ifndef CCTBX_SGTBX_MILLER_ASU_H
#define CCTBX_SGTBX_MILLER_ASU_H


#include <cctbx/miller.h>
#include <cctbx/sgtbx/groups.h>


namespace sgtbx
{

  //! Class for accessing default CCP4 reciprocal space ASUs
  /*!

    This class is initialised with a set of SgOps, and provides a
    single computationally efficient method for determining whether a
    Miller index is inside or outside the default CCP4 assymetric
    unit.

    The sgops are used to generate a reduced spacegroup, being the
    spacegroup of the diffraction pattern, or oriented Laue group.
    From this the tabulated spacegroup number and change-of-basis are
    determined. The 230 tabulated spacegroups reduce to 14 different
    spacegroups. Therefore ASU functions are implemented for each of
    these 14 groups. (11 Laue groups, plus 312 and two rhomohedral
    settings.

  */
  class CCP4_ReciprocalSpaceASU
  {
    //! Function pointer to the appropriate ASU function
    bool (CCP4_ReciprocalSpaceASU::*asufn)(const Miller::Index& H) const;
    //! Change of basis op to get back to the default setting
    RotMx cbm;
    //! reduced spacegroup (oriented Laue group)
    int rsg;

    //! Error function for uninitialised pointer
    bool ASU_uninit(const Miller::Index& r) const
      { throw error("CCP4_ReciprocalSpaceASU: uninitialised"); }

    //! internal ASU function for Laue group -1
    inline bool ASU_1b   (const Miller::Index& r) const
      { return (r[2]>0 || (r[2]==0 && (r[0]>0 || (r[0]==0 && r[1]>=0)))); }
    inline bool ASU_2_m  (const Miller::Index& r) const
      { return (r[1]>=0 && (r[2]>0 || (r[2]==0 && r[0]>=0))); }
    inline bool ASU_mmm  (const Miller::Index& r) const
      { return (r[0]>=0 && r[1]>=0 && r[2]>=0); }
    inline bool ASU_4_m  (const Miller::Index& r) const
      { return (r[2]>=0 && ((r[0]>=0 && r[1]>0) || (r[0]==0 && r[1]==0))); }
    inline bool ASU_4_mmm(const Miller::Index& r) const
      { return (r[0]>=r[1] && r[1]>=0 && r[2]>=0); }
    inline bool ASU_3b   (const Miller::Index& r) const
      { return (r[0]>=0 && r[1]>0) || (r[0]==0 && r[1]==0 && r[2] >= 0); }
    inline bool ASU_3bm  (const Miller::Index& r) const
      { return (r[0]>=r[1] && r[1]>=0 && (r[1]>0 || r[2]>=0)); }
    inline bool ASU_3bmx (const Miller::Index& r) const
      { return (r[0]>=r[1] && r[1]>=0 && (r[0]>r[1] || r[2]>=0)); }
    inline bool ASU_6_m  (const Miller::Index& r) const
      { return (r[2]>=0 && ((r[0]>=0 && r[1]>0) || (r[0]==0 && r[1]==0))); }
    inline bool ASU_6_mmm(const Miller::Index& r) const
      { return (r[0]>=r[1] && r[1]>=0 && r[2]>=0); }
    inline bool ASU_m3b  (const Miller::Index& r) const
      { return (r[0]>=0 && ((r[2]>=r[0] && r[1]>r[0]) || (r[2]==r[0] && r[1]==r[0]))); }
    inline bool ASU_m3bm (const Miller::Index& r) const
      { return (r[1]>=r[2] && r[2]>=r[0] && r[0]>=0); }

    //! internal ASU function for Laue group -1 with change-of-basis
    inline bool cASU_1b   (const Miller::Index& r) const
      { return ASU_1b   ( r * cbm ); }
    inline bool cASU_2_m  (const Miller::Index& r) const
      { return ASU_2_m  ( r * cbm ); }
    inline bool cASU_mmm  (const Miller::Index& r) const
      { return ASU_mmm  ( r * cbm ); }
    inline bool cASU_4_m  (const Miller::Index& r) const
      { return ASU_4_m  ( r * cbm ); }
    inline bool cASU_4_mmm(const Miller::Index& r) const
      { return ASU_4_mmm( r * cbm ); }
    inline bool cASU_3b   (const Miller::Index& r) const
      { return ASU_3b   ( r * cbm ); }
    inline bool cASU_3bm  (const Miller::Index& r) const
      { return ASU_3bm  ( r * cbm ); }
    inline bool cASU_3bmx (const Miller::Index& r) const
      { return ASU_3bmx ( r * cbm ); }
    inline bool cASU_6_m  (const Miller::Index& r) const
      { return ASU_6_m  ( r * cbm ); }
    inline bool cASU_6_mmm(const Miller::Index& r) const
      { return ASU_6_mmm( r * cbm ); }
    inline bool cASU_m3b  (const Miller::Index& r) const
      { return ASU_m3b  ( r * cbm ); }
    inline bool cASU_m3bm (const Miller::Index& r) const
      { return ASU_m3bm ( r * cbm ); }

   public:
    //! Null constructor
    CCP4_ReciprocalSpaceASU() { asufn = &CCP4_ReciprocalSpaceASU::ASU_uninit; }
    //! Constructor: takes SgOps and deduces Laue group and change-of-basis
    CCP4_ReciprocalSpaceASU(const SgOps& sgops) { init( sgops ); }
    //! Initialiser: takes SgOps and deduces Laue group and change-of-basis
    void init(const SgOps& sgops);
    //! Returns true if the reflection is in the default CCP4 reciprocal ASU
    /*! Note: only one of a pair of Friedel opposites will yield
      `true'. This method is highly optimised for standard settings,
      but will be slower for non-standard settings where a
      change-of-basis is required. */
    inline bool isInASU(const Miller::Index& r) const
      { return (this->*asufn)(r); }

    //! This returns one of exactly 11 codes.
    tables::MatrixGroup::Code getLaueGroupType() const;
    //! This returns something like "l>=0 and (k>0 or (k==0 and h>=0))"
    /*!  The string representation of the tabulated settings is true
      only for the standard settings. Otherwise it must be used in
      conjunction with getRotToStandardSetting(). */
    std::string getConditions() const;
    //! Return a RotMx rotation matrix representing the change-of-basis
    /*! getRotToStandardSetting().isUnit() will tell you if you have a
     standard setting or not. */
    RotMx getRotToStandardSetting() const;
    //! Similar to SgOps::getCutParameters.
    Miller::Vec3 getCutParameters() const;
  };

}


#endif // CCTBX_SGTBX_MILLER_ASU_H
