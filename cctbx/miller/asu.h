#ifndef CCTBX_MILLER_ASU_H
#define CCTBX_MILLER_ASU_H

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>

namespace cctbx { namespace miller {

  /*! \brief Support for common layouts of tables of
      asymmetric Miller indices.
   */
  /*! See asym_index for details.
   */
  class index_table_layout_adaptor : public sym_equiv_index
  {
    public:
      //! Default constructor. Some data members are not initialized!
      index_table_layout_adaptor() {}

      //! The entry in the table of asymmetric Miller indices.
      /*! See asym_index.
       */
      index<>
      h() const { return h_; }

      /*! \brief Index (0 or 1) of the column for the
          asym_index::two_column() layout.
       */
      std::size_t
      i_column() const { return i_column_; }

    private:
      friend class asym_index;

      index_table_layout_adaptor(
        sym_equiv_index const& eq,
        bool conj_h,
        bool f_minus_column,
        bool conj_f)
      :
        sym_equiv_index(eq), i_column_(0)
      {
        h_ = hr_;
        if (conj_h) h_ = -h_;
        if (f_minus_column) i_column_ = 1;
        friedel_flag_ = conj_f;
      }

      index<> h_;
      std::size_t i_column_;
  };

  //! Selection of an asymmetric Miller index.
  /*! The selection of a particular asymmetric
      Miller index amongst symmetrically equivalent Miller indices
      is based on class sgtbx::reciprocal_space::asu,
      or alternatively on an <i>ad hoc</i> scheme for the
      determination of the "prettiest" index for human-
      readable listings.
      <p>
      Friedel's law is always considered in the
      determination of the asymmetric index.
      friedel_flag() indicates if Friedel's law was
      actually applied.
      <p>
      Support for one- and two-column data
      layouts is provided using the class
      index_table_layout_adaptor:
      <p>
      <ul>
      <li>one_column()
      <li>two_column()
      </ul>
   */
  class asym_index : public sym_equiv_index
  {
    public:
      //! Default constructor. Some data members are not initialized!
      asym_index() {}

      //! Selection of an asymmetric Miller index.
      /*! The asymmetric index is determined as the product of the
          index h and the symmetry operation of sgtbx::space_group
          which maps h into the reciprocal space asymmetric unit
          (sgtbx::reciprocal_space::asu).
       */
      asym_index(sgtbx::space_group const& space_group,
                 sgtbx::reciprocal_space::asu const& asu,
                 index<> const& h);

      /*! \brief Selection of a "pretty" index for human-readable
          listings, given symmetry operations and an input index.
       */
      /*! The selection is based on miller::index<>::operator<().
          The asymmetric unit from which the indices are
          selected is not necessarily contiguous.
       */
      asym_index(sgtbx::space_group const& space_group,
                 index<> const& h);

      /*! \brief Selection of a "pretty" index for human-readable
          listings, given a list of symmetrically equivalent Miller
          indices.
       */
      /*! The selection is based on miller::index<>::operator<().
          The asymmetric unit from which the indices are
          selected is not necessarily contiguous.
       */
      asym_index(sym_equiv_indices const& h_eq);

      //! Adaptor for one-column table of asymmetric Miller indices.
      /*! anomalous_flag == false
          and only one value is stored for a Friedel pair:<pre>
           h  k  l  F</pre>
          The values associated with h,k,l and -h,-k,-l
          are assumed to be equal, and the phases are
          related by the equation
          phi(h,k,l) = -phi(-h,-k,-l).
          <p>
          anomalous_flag == true (i.e. in the presence of an
          anomalous signal):<pre>
            h  k  l  F
           -h -k -l  F</pre>
          There is a separate entry for each Friedel mate
          in a table.
       */
      index_table_layout_adaptor
      one_column(bool anomalous_flag) const
      {
        if (!anomalous_flag) {
          return index_table_layout_adaptor(*this,
            friedel_flag_, false, friedel_flag_);
        }
        return index_table_layout_adaptor(*this, false, false, false);
      }

      //! Adaptor for two-column table of asymmetric Miller indices.
      /*! anomalous_flag == false: same as
          one_column(). Only one column is used. Provided
          for completeness.
          <p>
          anomalous_flag == true:<pre>
            h  k  l  F+  F-</pre>
          Both Friedel mates are associated with the same
          index in a table.  The Miller index for F+ is
          (h, k, l) and the implied Miller index for F-
          is (-h, -k, -l).
       */
      index_table_layout_adaptor
      two_column(bool anomalous_flag) const
      {
        if (!anomalous_flag) {
          return index_table_layout_adaptor(*this,
            friedel_flag_, false, friedel_flag_);
        }
        return index_table_layout_adaptor(*this,
          friedel_flag_, friedel_flag_, false);
      }

      int
      isym() const { return isym_; }

    private:
      int isym_;

  };

  namespace data_classes {

    struct scalar_type {};
    struct complex_type {};
    struct hendrickson_lattman_type {};
    struct phase_type {};

  }

  template <typename DataClass>
  struct map_to_asu_policy;

  template <>
  struct map_to_asu_policy<data_classes::scalar_type>
  {
    template <typename ValueType>
    static
    void
    eq(index_table_layout_adaptor const& /*ila*/, ValueType& /*value*/, bool)
    {
    }
  };

  template <>
  struct map_to_asu_policy<data_classes::complex_type>
  {
    template <typename ValueType>
    static
    void
    eq(index_table_layout_adaptor const& ila, ValueType& value, bool)
    {
      value = ila.complex_eq(value);
    }
  };

  template <>
  struct map_to_asu_policy<data_classes::hendrickson_lattman_type>
  {
    template <typename ValueType>
    static
    void
    eq(index_table_layout_adaptor const& ila, ValueType& value, bool)
    {
      value = ila.hendrickson_lattman_eq(value);
    }
  };

  template <>
  struct map_to_asu_policy<data_classes::phase_type>
  {
    template <typename ValueType>
    static
    void
    eq(index_table_layout_adaptor const& ila, ValueType& value, bool deg)
    {
      value = ila.phase_eq(value, deg);
    }
  };

  template <> struct map_to_asu_policy<double>
  : map_to_asu_policy<data_classes::scalar_type> {};

  template <> struct map_to_asu_policy<std::complex<double> >
  : map_to_asu_policy<data_classes::complex_type> {};

  template <> struct map_to_asu_policy<hendrickson_lattman<> >
  : map_to_asu_policy<data_classes::hendrickson_lattman_type> {};

  namespace detail {

    template <typename ValueType,
              typename PolicySelectType>
    void
    map_to_asu(
      sgtbx::space_group_type const& sg_type,
      bool anomalous_flag,
      af::ref<index<> > const& miller_indices,
      af::ref<ValueType> const& data,
      bool deg)
    {
      CCTBX_ASSERT(miller_indices.size() == data.size());
      sgtbx::reciprocal_space::asu asu(sg_type);
      sgtbx::space_group const& sg = sg_type.group();
      for(std::size_t i=0;i<miller_indices.size();i++) {
        asym_index ai(sg, asu, miller_indices[i]);
        index_table_layout_adaptor ila = ai.one_column(anomalous_flag);
        miller_indices[i] = ila.h();
        map_to_asu_policy<PolicySelectType>::eq(ila, data[i], deg);
      }
    }

  } // namespace detail

  void
  map_to_asu(
    sgtbx::space_group_type const& sg_type,
    bool anomalous_flag,
    af::ref<index<> > const& miller_indices);

  void
  map_to_asu_isym(
    sgtbx::space_group_type const& sg_type,
    bool anomalous_flag,
    af::ref<index<> > const& miller_indices,
    af::ref<int> const& isym);

  template <typename ValueType>
  void
  map_to_asu(
    sgtbx::space_group_type const& sg_type,
    bool anomalous_flag,
    af::ref<index<> > const& miller_indices,
    af::ref<ValueType> const& data)
  {
    detail::map_to_asu<ValueType, ValueType>(
      sg_type, anomalous_flag, miller_indices, data, false);
  }

  template <typename ValueType>
  void
  map_to_asu(
    sgtbx::space_group_type const& sg_type,
    bool anomalous_flag,
    af::ref<index<> > const& miller_indices,
    af::ref<ValueType> const& data,
    bool deg)
  {
    detail::map_to_asu<ValueType, data_classes::phase_type>(
      sg_type, anomalous_flag, miller_indices, data, deg);
  }

  bool
  is_unique_set_under_symmetry(
    sgtbx::space_group_type const& space_group_type,
    bool anomalous_flag,
    af::const_ref<index<> > const& miller_indices);

  af::shared<std::size_t>
  unique_under_symmetry_selection(
    sgtbx::space_group_type const& space_group_type,
    bool anomalous_flag,
    af::const_ref<index<> > const& miller_indices);

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_ASU_H
