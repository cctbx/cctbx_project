// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/sgtbx/miller_asu.h (rwgk)
 */

#ifndef CCTBX_MILLER_ASU_H
#define CCTBX_MILLER_ASU_H

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/sgtbx/miller_asu.h>

namespace cctbx { namespace miller {

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
      AsymIndex(SymEquivIndices const& SEMI);
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

  template <typename DataType>
  class map_to_asu
  {
    public:
      map_to_asu() {}

      map_to_asu(
        const sgtbx::SpaceGroupInfo& sginfo,
        bool friedel_flag,
        af::shared<miller::Index> miller_indices,
        af::shared<DataType> data_array,
        bool in_place = false)
        : friedel_flag_(friedel_flag),
          asu_(sginfo)
      {
        cctbx_assert(miller_indices.size() == data_array.size());
        const sgtbx::SpaceGroup& sgops = sginfo.SgOps();
        for(std::size_t i=0;i<miller_indices.size();i++) {
          AsymIndex ai(sgops, asu_, miller_indices[i]);
          IndexTableLayoutAdaptor ila = ai.one_column(friedel_flag);
          if (in_place) {
            miller_indices[i] = ila.H();
            data_array[i] = ila.complex_eq(data_array[i]);
          }
          else {
            asym_indices_.push_back(ila.H());
            asym_data_array_.push_back(ila.complex_eq(data_array[i]));
          }
        }
        if (in_place) {
          asym_indices_ = miller_indices;
          asym_data_array_ = data_array;
        }
      }

      bool friedel_flag() const { return friedel_flag_; }

      sgtbx::ReciprocalSpaceASU const& asu() const { return asu_; }

      af::shared<miller::Index> asym_indices() const {
        return asym_indices_;
      }

      af::shared<DataType> asym_data_array() const {
        return asym_data_array_;
      }
    protected:
      bool friedel_flag_;
      sgtbx::ReciprocalSpaceASU asu_;
      af::shared<miller::Index> asym_indices_;
      af::shared<DataType> asym_data_array_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_ASU_H
