// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#ifndef CCTBX_SGTBX_MILLER_ASU_H
#define CCTBX_SGTBX_MILLER_ASU_H

#include <cctbx/miller.h>
#include <cctbx/sgtbx/groups.h>

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
      This implementation is based on work by
      <a href="http://www.ysbl.york.ac.uk/~cowtan/"
      >Kevin Cowtan</a>.
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
          !SpaceGroup::isCentric(H) to determine which reflections
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
      space group type (class SpaceGroupInfo) is used to derive a
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
      /*! Based on the space group number (SpaceGroupInfo::SgNumber()),
          the Laue group class is derived which is in turn used
          to select the corresponding tabulated
          ReferenceReciprocalSpaceASU.
       */
      ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo);
      //! Access to the selected tabulated ReferenceReciprocalSpaceASU.
      inline const ReferenceReciprocalSpaceASU* ReferenceASU() const {
        return m_ReferenceASU;
      }
      //! Access to the change-of-basis operator.
      /*! This operator is a copy of SgInfo.CBOp() as passed to
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
      an asymmetric unit of Miller indices up to a given
      high-resolution limit or up to a given maximum Miller
      index.
   */
  /*! Example (Python syntax):<pre>
        # Given a resolution limit.
        UnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
        SgOps = sgtbx.SpaceGroup("P 41")
        MIG = sgtbx.MillerIndexGenerator(UnitCell, SgOps.Info(), 3.0)
        for H in MIG: print H
        #
        # Given a maximum Miller index.
        SgOps = sgtbx.SpaceGroup("P 31")
        MIG = sgtbx.MillerIndexGenerator(SgOps.Info(), (3,4,5))
        for H in MIG: print H
      </pre>
      This class is implemented as an iterator. Therefore
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
      //! Initialization with resolution limit.
      /*! Miller indices up to and including Resolution_d_min will
          be generated.
       */
      MillerIndexGenerator(const uctbx::UnitCell& uc,
                           const SpaceGroupInfo& SgInfo,
                           double Resolution_d_min);
      //! Initialization with maximum Miller index.
      /*! Miller indices in the range from -MaxIndex to +MaxIndex
          will be generated.
       */
      MillerIndexGenerator(const SpaceGroupInfo& SgInfo,
                           const Miller::Index& MaxIndex);
      //! Access to the reciprocal space asymmetric unit.
      /*! The Miller indices that are generated by this class (member
          function next()) are inside this asymmetric unit.
       */
      inline const ReciprocalSpaceASU& ASU() const { return m_ASU; }
      //! Iterator over Miller indices.
      /*! Each call to this member function will return the next
          Miller index in the sequence. The indices are inside
          ReciprocalSpaceASU(). Systematically absent reflections
          are automatically filtered out.
          <p>
          The Miller index 0,0,0 indicates the end of the iteration.
       */
      Miller::Index next();
    private:
      void InitializeLoop(const Miller::Index& ReferenceHmax);
      uctbx::UnitCell m_UnitCell;
      SpaceGroup m_SgOps;
      ReciprocalSpaceASU m_ASU;
      double m_Qhigh;
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
            H and the symmetry operation of SpaceGroup which maps H
            into the reciprocal space asymmetric unit
            (ReciprocalSpaceASU).
         */
        SymUniqueIndex(const sgtbx::SpaceGroup& SgOps,
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
