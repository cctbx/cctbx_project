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

} // namespace sgtbx

namespace sgtbx {

  /*! \brief Contiguous reciprocal space asymmetric units for
      the 230 reference settings.
   */
  /*! 12 contiguous reciprocal space asymmetric units (11 Laue
      classes, two settings for Laue class -3m) are
      tabulated. The tabulated asymmetric units are
      compatible with the asymmetric units of the CCP4
      package.
      <p>
      This implementation is based on work by Kevin Cowtan.
      <p>
      See also: class ReciprocalSpaceASU, class MillerIndexGenerator
   */
  class ReferenceReciprocalSpaceASU {
    public:
      //! Returns one of exactly 12 Laue group codes.
      /*! The labels of the possible return codes are:<br>
            -1, 2/m, mmm, 4/m, 4/mmm, -3, -31m, -3m1, 6/m, 6/mmm, m-3, m-3m
          <p>
          For Laue class -3m there are two possible orientations of the
          mirror plane with respect to the periodic lattice.
       */
      virtual tables::MatrixGroup::Code LaueGroupCode() const {
        throw cctbx_internal_error();
      }
      //! Test if given Miller index is in the tabulated asymmetric unit.
      virtual bool isInASU(const Miller::Index& H) const {
        throw cctbx_internal_error();
      }
      //! String representation of the tabluated asymmetric unit.
      /*! Example: "h>=k and k>=0 and (k>0 or l>=0)"
       */
      virtual const char* representation() const {
        throw cctbx_internal_error();
      }
      //! "Cut parameters" for building Miller indices.
      /*! When building (or generating) a large list of Miller indices,
          it is useful to restrict the loop over all possible indices
          to 1/2, 1/4, or 1/8 of reciprocal space, if possible.
          <p>
          The cut parameters are used in the next() method
          of the class MillerIndexGenerator. In general it
          should be much more convenient to use that higher-level
          class rather than hand-crafting a loop for building
          Miller indices.
          <p>
          Each element of CutP is either -1 or 0. A value
          of 0 indicates that the corresponding negative half-space
          can be omitted in the loop over possible indices.
          <p>
          Friedel symmetry is implied. If the Friedel mates
          are needed explicitly, they have to be added in a
          separate step. Note that the Friedel mate appears
          explicitly only for acentric reflections (use e.g.
          !SgOps::isCentric(H) to determine which reflections
          are acentric).
       */
      virtual const Miller::Vec3& getCutParameters() const {
        throw cctbx_internal_error();
      }
  };

  //! Access to general contiguous reciprocal space asymmetric units.
  /*! class ReferenceReciprocalSpaceASU implements 12 contiguous
      reciprocal space asymmetric units that cover the 230
      reference settings. The algorithm for the determination of the
      space group type (SgOps::getSpaceGroupType) is used to derive a
      change-of-basis matrix for the transformation of the tabulated
      asymmetric units. In this way a contiguous asymmetric unit is
      available for any arbitrary setting.
   */
  class ReciprocalSpaceASU
  {
    public:
      //! Default constructor.
      /*! Default constructed instances will throw exceptions if
          some of the member functions are used.
       */
      inline ReciprocalSpaceASU()
        : m_CBOp(), m_isReferenceASU(true), m_ReferenceASU() {}
      //! Initialization.
      /*! Based on the space group number (SpaceGroupType::SgNumber),
          the Laue group class is derived which is in turn used
          to select the corresponding tabulated
          ReferenceReciprocalSpaceASU.
       */
      ReciprocalSpaceASU(const SpaceGroupType& SgType);
      //! Access to the selected tabulated ReferenceReciprocalSpaceASU.
      inline const ReferenceReciprocalSpaceASU* ReferenceASU() const {
        return m_ReferenceASU;
      }
      //! Access to the change-of-basis operator.
      /*! This operator is a copy of SgType.CBOp() as passed to
          the constructor.
       */
      inline const ChOfBasisOp& CBOp() const { return m_CBOp; }
      //! Test if the given asymmetric unit is one of the tabulated units.
      /*! This test is equivalent to the test CBOp().M().Rpart().isUnit().
          That is, it is tested if the rotation part of the
          change-of-basis matrix is the unit matrix. (If this
          is the case, some optimizations are possible.)
       */
      inline bool isReferenceASU() const { return m_isReferenceASU; }
      //! Test if the given Miller index is in the asymmetric unit.
      /*! The change-of-basis matrix is used to transform
          the Miller index (H * CBOp().InvM()). It is then
          tested if the result is in the tabulated reference
          asymmetric unit.
       */
      inline bool isInASU(const Miller::Index& H) const {
        if (m_isReferenceASU) return m_ReferenceASU->isInASU(H);
        return m_ReferenceASU->isInASU(H * m_CBOp.InvM().Rpart());
      }
    private:
      ChOfBasisOp m_CBOp;
      bool m_isReferenceASU;
      const ReferenceReciprocalSpaceASU* m_ReferenceASU;
  };

  /*! \brief Efficient, easy-to-use algorithm for building
      Miller indices up to a given high-resolution limit.
   */
  /*! Example (Python syntax):<pre>
        UnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
        MIG = sgtbx.MillerIndexGenerator(UnitCell, SgOps, 3.0)
        for H in MIG: print H
      </pre>
      This class is implemented as an iterator. That is,
      the generation of Miller indices does not consume any
      significant amounts of memory. The key to efficiency
      is class ReferenceReciprocalSpaceASU.
   */
  class MillerIndexGenerator
  {
    public:
      //! Default constructor.
      /*! Default constructed instances will throw exceptions if
          some of the member functions are used.
       */
      MillerIndexGenerator() {}
      //! Initialization.
      /*! Miller indices up to and including Resolution_d_min will
          be generated.
       */
      MillerIndexGenerator(const uctbx::UnitCell& uc,
                           const SgOps& sgo,
                           double Resolution_d_min);
      //! Access to the reciprocal space asymmetric unit.
      /*! The Miller indices that are generated by this class (member
          function next()) are inside this asymmetric unit.
       */
      inline const ReciprocalSpaceASU& ASU() const { return m_ASU; }
      //! Iterator over Miller indices.
      /*! Each call to this member function will return the next
          Miller index in the sequence. The indices are inside
          ReciprocalSpaceASU().
       */
      Miller::Index next();
    private:
      uctbx::UnitCell m_UnitCell;
      SgOps m_SgOps;
      double m_Qhigh;
      ReciprocalSpaceASU m_ASU;
      NestedLoop<Miller::Index> m_loop;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_MILLER_ASU_H
