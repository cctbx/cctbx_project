// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 27: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#ifndef CCTBX_SGTBX_MILLER_ASU_H
#define CCTBX_SGTBX_MILLER_ASU_H

#include <cctbx/array_family/loops.h>
#include <cctbx/sgtbx/groups.h>

namespace cctbx {
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
        virtual const af::int3& getCutParameters() const {
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
        ReciprocalSpaceASU()
          : m_CBOp(), m_isReferenceASU(true), m_ReferenceASU() {}
        //! Initialization.
        /*! Based on the space group number (SpaceGroupInfo::SgNumber()),
            the Laue group class is derived which is in turn used
            to select the corresponding tabulated
            ReferenceReciprocalSpaceASU.
         */
        ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo);
        //! Access to the selected tabulated ReferenceReciprocalSpaceASU.
        const ReferenceReciprocalSpaceASU* ReferenceASU() const {
          return m_ReferenceASU;
        }
        //! Access to the change-of-basis operator.
        /*! This operator is a copy of SgInfo.CBOp() as passed to
            the constructor.
         */
        const ChOfBasisOp& CBOp() const { return m_CBOp; }
        //! Test if the given asymmetric unit is one of the tabulated units.
        /*! This test is equivalent to the test CBOp().M().Rpart().isUnit().
            That is, it is tested if the rotation part of the
            change-of-basis matrix is the unit matrix. (If this
            is the case, some optimizations are possible.)
         */
        bool isReferenceASU() const { return m_isReferenceASU; }
        //! Test if the given Miller index is in the asymmetric unit.
        /*! The change-of-basis matrix is used to transform
            the Miller index (H * CBOp().InvM()). It is then
            tested if the result is in the tabulated reference
            asymmetric unit.
         */
        bool isInASU(const Miller::Index& H) const {
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
        const ReciprocalSpaceASU& ASU() const { return m_ASU; }
        //! Iterator over Miller indices.
        /*! Each call to this member function will return the next
            Miller index in the sequence. The indices are inside
            ReciprocalSpaceASU(). Systematically absent reflections
            are automatically filtered out.
            <p>
            The Miller index 0,0,0 indicates the end of the iteration.
         */
        Miller::Index next();
        //! Add all Miller indices to an array.
        /*! The next() method is called in a loop until the list
            of Miller indices is exhausted. The results are added
            to ArrayOfH.
            <p>
            Requirements for MillerArrayType:
            <ul>
            <li>Must contain objects of type Miller::Index.
            <li>Must support the push_back() method.
            </ul>
            <p>
            Example:<pre>
            cctbx::uctbx::UnitCell UC(...);
            cctbx::sgtbx::SpaceGroup SgOps(...);
            cctbx::af::shared<Miller::Index> ArrayOfH;
            MillerIndexGenerator(UC, SgOps, 1.0).AddtoArray(ArrayOfH);
            </pre>
         */
        template <class MillerArrayType>
        void
        AddToArray(MillerArrayType& ArrayOfH)
        {
          for (;;) {
            Miller::Index H = next();
            if (H.is000()) break;
            ArrayOfH.push_back(H);
          }
        }
      private:
        void InitializeLoop(const Miller::Index& ReferenceHmax);
        uctbx::UnitCell m_UnitCell;
        SpaceGroup m_SgOps;
        ReciprocalSpaceASU m_ASU;
        double m_Qhigh;
        af::nested_loop<Miller::Index> m_loop;
    };

  } // namespace sgtbx

  namespace Miller {

    /*! \brief Support for common layouts of tables of
        asymmetric Miller indices.
     */
    /*! See AsymIndex for details.
     */
    class IndexTableLayoutAdaptor : private SymEquivIndex
    {
      public:
        //! Default constructor. Some data members are not initialized!
        IndexTableLayoutAdaptor() {}
        //! The entry in the table of asymmetric Miller indices.
        /*! See AsymIndex.
         */
        Index H() const {
          return m_H;
        }
        /*! \brief %Index (0 or 1) of the column for the
            AsymIndex::FplusFminusLayout().
         */
        int iColumn() const {
          return m_iColumn;
        }
        //! Phase in radians, given phase for input Miller index.
        template <class FloatType>
        FloatType
        Phase_rad(const FloatType& phi) const {
          if (m_ConjF) return -SymEquivIndex::Phase_rad(phi);
          return SymEquivIndex::Phase_rad(phi);
        }
        //! Phase in degrees, given phase for input Miller index.
        template <class FloatType>
        FloatType
        Phase_deg(const FloatType& phi) const {
          if (m_ConjF) return -SymEquivIndex::Phase_deg(phi);
          return SymEquivIndex::Phase_deg(phi);
        }
        /*! \brief Complex structure factor, given structure factor
            for input Miller index.
         */
        template <class FloatType>
        std::complex<FloatType>
        ShiftPhase(const std::complex<FloatType>& F) const {
          if (m_ConjF) return std::conj(SymEquivIndex::ShiftPhase(F));
          return SymEquivIndex::ShiftPhase(F);
        }
      private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
        friend class AsymIndex;
#endif // DOXYGEN_SHOULD_SKIP_THIS
        IndexTableLayoutAdaptor(const SymEquivIndex& SEI,
                                bool ConjH, bool FminusColumn, bool ConjF)
          : SymEquivIndex(SEI), m_iColumn(0), m_ConjF(ConjF)
        {
          m_H = m_HR;
          if (ConjH) m_H = m_H.FriedelMate();
          if (FminusColumn) m_iColumn = 1;
        }
      private:
        Index m_H;
        int m_iColumn;
        bool m_ConjF;
    };

    //! Selection of an asymmetric Miller index.
    /*! The selection of a particular asymmetric
        Miller index amongst symmetrically equivalent Miller indices
        is based on class sgtbx::ReciprocalSpaceASU,
        or alternatively on an <i>ad hoc</i> scheme for the
        determination of the "prettiest" index for human-
        readable listings.
        <p>
        Friedel's law is always considered in the
        determination of the asymmetric index.
        FriedelFlag() indicates if Friedel's law was
        actually applied.
        <p>
        Support for the following data
        layouts is provided using the class
        IndexTableLayoutAdaptor:
        <p>
        <ul>
        <li>AnomalousLayout()
        <li>HermitianLayout()
        <li>FplusFminusLayout()
        </ul>
     */
    class AsymIndex : public SymEquivIndex
    {
      public:
        //! Default constructor. Some data members are not initialized!
        AsymIndex() {}
        //! Selection of an asymmetric Miller index.
        /*! The asymmetric index is determined as the product of the
            Index H and the symmetry operation of sgtbx::SpaceGroup
            which maps H into the reciprocal space asymmetric unit
            (sgtbx::ReciprocalSpaceASU).
         */
        AsymIndex(const sgtbx::SpaceGroup& SgOps,
                  const sgtbx::ReciprocalSpaceASU& ASU,
                  const Index& H);
        /*! \brief Selection of a "pretty" Index for human-readable
            listings, given symmetry operations and an input Index.
         */
        /*! The selection is based on Miller::Index::operator<().
            The asymmetric unit from which the indices are
            selected is not necessarily contiguous.
         */
        AsymIndex(const sgtbx::SpaceGroup& SgOps,
                  const Index& H);
        /*! \brief Selection of a "pretty" Index for human-readable
            listings, given a list of symmetrically equivalent %Miller
            indices.
         */
        /*! The selection is based on Miller::Index::operator<().
            The asymmetric unit from which the indices are
            selected is not necessarily contiguous.
         */
        AsymIndex(const sgtbx::SymEquivMillerIndices& SEMI);
        //! Flag for application of Friedel's law.
        /*! This flag indicates if Friedel's law was applied
            in the selection of the index. If true, the
            selected index is the Friedel mate of one of
            the indices that are symmetrically equivalent
            to the input index. (For centric reflections,
            FriedelFlag() is always false.)
            <p>
            Note that support for common layouts
            of tables of asymmetric Miller indices is
            included in this class (AnomalousLayout(),
            HermitianLayout(), FplusFminusLayout()).
            For these layouts it should not be necessary to
            inspect FriedelFlag() directly.
         */
        bool FriedelFlag() const {
          return m_FriedelFlag;
        }
        //! Adaptor for table of asymmetric Miller indices.
        /*! No Friedel symmetry (i.e. in the presence of an
            anomalous signal), one column of data:<pre>
              h  k  l  F
             -h -k -l  F</pre>
            There is a separate entry for each Friedel mate
            in a table.
         */
        IndexTableLayoutAdaptor AnomalousLayout() const {
          return IndexTableLayoutAdaptor(*this,
            false, false, false);
        }
        //! Adaptor for table of asymmetric Miller indices.
        /*! Assuming Friedel symmetry (no anomalous signal),
            and only one value is stored for a Friedel pair:<pre>
             h  k  l  F</pre>
            The values associated with h,k,l and -h,-k,-l
            are assumed to be equal, and the phases are
            related by the equation
            phi(h,k,l) = -phi(-h,-k,-l).
         */
        IndexTableLayoutAdaptor HermitianLayout() const {
          return IndexTableLayoutAdaptor(*this,
            m_FriedelFlag, false, m_FriedelFlag);
        }
        //! Adaptor for table of asymmetric Miller indices.
        /*! No Friedel symmetry (i.e. the in presence of an anomalous
            signal), two columns of data:<pre>
              h  k  l  F+  F-</pre>
            Both Friedel mates are associated with the same
            index in a table.  The Miller index for F+ is
            (h, k, l) and the implied Miller index for F-
            is (-h, -k, -l).
         */
        IndexTableLayoutAdaptor FplusFminusLayout() const {
          return IndexTableLayoutAdaptor(*this,
            m_FriedelFlag, m_FriedelFlag, false);
        }
      private:
        bool m_FriedelFlag;
    };

  } // namespace Miller
} // namespace cctbx

#endif // CCTBX_SGTBX_MILLER_ASU_H
