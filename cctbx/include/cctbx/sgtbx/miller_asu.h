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

#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_ref_asu.h>

namespace cctbx {
  namespace sgtbx {

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
        bool isInASU(const miller::Index& H) const {
          if (m_isReferenceASU) return m_ReferenceASU->isInASU(H);
          return m_ReferenceASU->isInASU(H * m_CBOp.InvM().Rpart());
        }
        int asu_sign(const miller::Index& h,
                     const miller::Index& minus_h) const {
          if (m_isReferenceASU) return m_ReferenceASU->asu_sign(h, minus_h);
          miller::Index h_ref = h * m_CBOp.InvM().Rpart();
          return m_ReferenceASU->asu_sign(h_ref);
        }
        int asu_sign(const miller::Index& h) const {
          if (m_isReferenceASU) return m_ReferenceASU->asu_sign(h);
          miller::Index h_ref = h * m_CBOp.InvM().Rpart();
          return m_ReferenceASU->asu_sign(h_ref);
        }
      private:
        ChOfBasisOp m_CBOp;
        bool m_isReferenceASU;
        const ReferenceReciprocalSpaceASU* m_ReferenceASU;
    };

  } // namespace sgtbx

  namespace miller {

    /*! \brief Support for common layouts of tables of
        asymmetric Miller indices.
     */
    /*! See AsymIndex for details.
     */
    class IndexTableLayoutAdaptor : public SymEquivIndex
    {
      public:
        //! Default constructor. Some data members are not initialized!
        IndexTableLayoutAdaptor() {}
        //! The entry in the table of asymmetric Miller indices.
        /*! See AsymIndex.
         */
        Index H() const { return m_H; }
        /*! \brief %Index (0 or 1) of the column for the
            AsymIndex::TwoColumn() layout.
         */
        int iColumn() const { return m_iColumn; }
      private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
        friend class AsymIndex;
#endif // DOXYGEN_SHOULD_SKIP_THIS
        IndexTableLayoutAdaptor(const SymEquivIndex& SEI,
                                bool ConjH, bool FminusColumn, bool ConjF)
          : SymEquivIndex(SEI), m_iColumn(0)
        {
          m_H = m_HR;
          if (ConjH) m_H = -m_H;
          if (FminusColumn) m_iColumn = 1;
          m_FriedelFlag = ConjF;
        }

        Index m_H;
        int m_iColumn;
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
        Support one- and two-column data
        layouts is provided using the class
        IndexTableLayoutAdaptor:
        <p>
        <ul>
        <li>one_column()
        <li>two_column()
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
        /*! The selection is based on miller::Index::operator<().
            The asymmetric unit from which the indices are
            selected is not necessarily contiguous.
         */
        AsymIndex(const sgtbx::SpaceGroup& SgOps,
                  const Index& H);
        /*! \brief Selection of a "pretty" Index for human-readable
            listings, given a list of symmetrically equivalent Miller
            indices.
         */
        /*! The selection is based on miller::Index::operator<().
            The asymmetric unit from which the indices are
            selected is not necessarily contiguous.
         */
        AsymIndex(const sgtbx::SymEquivMillerIndices& SEMI);
        //! Adaptor for one-column table of asymmetric Miller indices.
        /*! FriedelFlag == true (no anomalous signal),
            and only one value is stored for a Friedel pair:<pre>
             h  k  l  F</pre>
            The values associated with h,k,l and -h,-k,-l
            are assumed to be equal, and the phases are
            related by the equation
            phi(h,k,l) = -phi(-h,-k,-l).
            <p>
            FriedelFlag == false (i.e. in the presence of an
            anomalous signal):<pre>
              h  k  l  F
             -h -k -l  F</pre>
            There is a separate entry for each Friedel mate
            in a table.
         */
        IndexTableLayoutAdaptor one_column(bool friedel_flag) const
        {
          if (friedel_flag) {
            return IndexTableLayoutAdaptor(*this,
              m_FriedelFlag, false, m_FriedelFlag);
          }
          return IndexTableLayoutAdaptor(*this, false, false, false);
        }
        //! Adaptor for two-column table of asymmetric Miller indices.
        /*! FriedelFlag == true (no anomalous signal): same as
            OneColumn(). Only one column is used. Provided
            for completeness.
            <p>
            FriedelFlag == false (i.e. the in presence of an anomalous
            signal):<pre>
              h  k  l  F+  F-</pre>
            Both Friedel mates are associated with the same
            index in a table.  The Miller index for F+ is
            (h, k, l) and the implied Miller index for F-
            is (-h, -k, -l).
         */
        IndexTableLayoutAdaptor two_column(bool friedel_flag) const
        {
          if (friedel_flag) {
            return IndexTableLayoutAdaptor(*this,
              m_FriedelFlag, false, m_FriedelFlag);
          }
          return IndexTableLayoutAdaptor(*this,
            m_FriedelFlag, m_FriedelFlag, false);
        }
    };

  } // namespace miller
} // namespace cctbx

#endif // CCTBX_SGTBX_MILLER_ASU_H
