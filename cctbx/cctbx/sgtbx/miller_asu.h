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
      /*! Default-constructed instances will throw exceptions if
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
      /*! Default-constructed instances will throw exceptions if
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

namespace cctbx {
  namespace Miller {

    //! Symmetry-unique ("asymmetric") Miller index class.
    /*! The selection of the symmetry-unique ("asymmetric")
        Miller index is based on class sgtbx::ReciprocalSpaceASU.
        <p>
        The SymUniqueIndex class supports the following data layouts
        for lists of symmetry-unique ("asymmetric") Miller indices:
        <p>
        1. Assuming Friedel symmetry (no anomalous signal), only one
        value is stored for a Friedel pair:<pre>
         h  k  l  F</pre>
        The values associated with h,k,l and -h,-k,-l are assumed to be
        equal, and the phases are related by the equation phi(h,k,l) =
        -phi(-h,-k,-l).<br>
        In this case, SymUniqueIndex.H() is intended for use as the
        symmetry-unique index in a list.
        <p>
        2. No Friedel symmetry (i.e. in presence of an anomalous signal),
        two columns of data:<pre>
         h  k  l  F+  F-</pre>
        Both Friedel mates are associated with the same index in a list.
        The Miller index for F+ is (h, k, l) and the implied Miller index
        for F- is (-h, -k, -l).<br>
        In this case, SymUniqueIndex.H() is intended for use as the
        symmetry-unique index in a list.
        SymUniqueIndex.iMate() is intended for selecting the data column.
        <p>
        3. No Friedel symmetry (i.e. in presence of an anomalous signal),
        one column of data:<pre>
         h  k  l  F
        -h -k -l  F</pre>
        There is a separate entry for each Friedel mate in a list.<br>
        In this case, SymUniqueIndex.HR() is intended for use as the
        symmetry-unique index in a list.
     */
    class SymUniqueIndex
    {
      public:
        //! Default constructor.
        /*! The member functions of default-constructed instances
            will return meaningless values.
         */
        inline SymUniqueIndex() : m_TBF(0) {}
        //! Initialization.
        /*! The symmetry-unique Miller index is determined in
            the constructor as the product of the Miller index
            H and the symmetry operation of SgOps which maps H
            into the reciprocal space asymmetric unit
            (ReciprocalSpaceASU).
         */
        SymUniqueIndex(const sgtbx::SgOps& sgo,
                       const sgtbx::ReciprocalSpaceASU& ASU,
                       const Index& H);
        /*! \brief The symmetry-unique index for data layouts 1 and 2.
            See class details.
         */
        inline const Index& H() const { return m_H; }
        //! Selection of F+ or F- data column. See class details.
        /*! The values returned are 0 for F+, and 1 for F-.
         */
        inline int iMate() const { return m_iMate; }
        //! The symmetry-unique index for data layout 3. See class details.
        /*! This index is the product of the input Miller index H that was
            passed to the constructor, and the rotation part of the
            symmetry operation that maps H into the asymmetric unit
            that was passed to the constructor.
         */
        inline const Index& HR() const { return m_HR; }
        //! Phase shift for HR() with respect to the input Miller index.
        /*! Low-level information for computing the phase of the
            symmetry-unique index, given the phase of the input Miller
            index. Note that high-level functions are also available
            (e.g. ShiftPhase()).<br>
            HT() is multiplied by the base factor TBF() in order to
            obtain an integer value.
         */
        inline int HT() const { return m_HT; }
        //! Translation base factor.
        /*! This is the factor by which HT() is multiplied.
         */
        inline int TBF() const { return m_TBF; }
        /*! \brief Phase for symmetry-unique index in radians, given
            phase for input Miller index.
         */
        /*! Formula used:<br>
            symunique_phi = phi - (2 * pi * HT()) / TBF();<br>
            if (FriedelSym && iMate()) symunique_phi = -symunique_phi;
         */
        template <class T>
        T Phase_rad(const T& phi, bool FriedelSym) const {
          using cctbx::constants::pi;
          T symunique_phi = phi - (2. * pi * HT()) / TBF();
          if (FriedelSym && iMate()) return -symunique_phi;
          return symunique_phi;
        }
        /*! \brief Phase for symmetry-unique index in degrees, given
            phase for input Miller index.
         */
        /*! Formula used:<br>
            symunique_phi = phi - (2 * 180 * HT()) / TBF();<br>
            if (FriedelSym && iMate()) symunique_phi = -symunique_phi;
         */
        template <class T>
        T Phase_deg(const T& phi, bool FriedelSym) const {
          T symunique_phi = phi - (2. * 180. * HT()) / TBF();
          if (FriedelSym && iMate()) return -symunique_phi;
          return symunique_phi;
        }
        /*! \brief Complex value for symmetry-unique index, given complex
            value for input index.
         */
        /*! Formula used:<br>
            symunique_F = F * exp(-2 * pi * j * HT() / TBF());<br>
            where j is the imaginary number.
         */
        template <class T>
        std::complex<T> ShiftPhase(const std::complex<T>& F,
                                   bool FriedelSym) const {
          using cctbx::constants::pi;
          T theta = (-2. * pi * HT()) / TBF();
          std::complex<T> symunique_F = F * std::polar(1., theta);
          if (FriedelSym && iMate()) return std::conj(symunique_F);
          return symunique_F;
        }

      private:
        int m_TBF;
        Index m_HR;
        int m_HT;
        Index m_H;
        bool m_iMate;
    };

  } // namespace Miller
} // namespace cctbx

#endif // CCTBX_SGTBX_MILLER_ASU_H
