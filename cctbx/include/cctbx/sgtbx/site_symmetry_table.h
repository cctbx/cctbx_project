#ifndef CCTBX_SGTBX_SITE_SYMMETRY_TABLE_H
#define CCTBX_SGTBX_SITE_SYMMETRY_TABLE_H

#include <cctbx/sgtbx/site_symmetry.h>

namespace cctbx { namespace sgtbx {

  class site_symmetry_table
  {
    public:
      //! Default constructor.
      site_symmetry_table()
      :
        indices_const_ref_(indices_.const_ref()),
        table_const_ref_(table_.const_ref())
      {}

      //! Add new site_symmetry_ops to internal table.
      /*! The internal table is searched for identical site_symmetry_ops.
          The indices() for duplicate site_symmetry_ops will point to the
          same table() entry.
       */
      void
      process(site_symmetry_ops const& site_symmetry_ops_)
      {
        CCTBX_ASSERT(indices_const_ref_.end() == indices_.end());
        CCTBX_ASSERT(table_const_ref_.end() == table_.end());
        if (table_const_ref_.size() == 0) {
          table_.push_back(site_symmetry_ops_.make_point_group_1());
          table_const_ref_ = table_.const_ref();
        }
        if (site_symmetry_ops_.is_point_group_1()) {
          indices_.push_back(0);
        }
        else {
          std::size_t i;
          for(i=0;i<table_const_ref_.size();i++) {
            if (   table_const_ref_[i].special_op()
                != site_symmetry_ops_.special_op()) {
              continue;
            }
            af::const_ref<rt_mx> tm=table_const_ref_[i].matrices().const_ref();
            af::const_ref<rt_mx> nm=site_symmetry_ops_.matrices().const_ref();
            if (tm.size() != nm.size()) continue;
            std::size_t j;
            for(j=0;j<tm.size();j++) {
              if (tm[j] != nm[j]) break;
            }
            if (j == tm.size()) break;
          }
          if (i == table_const_ref_.size()) {
            table_.push_back(site_symmetry_ops_);
            table_const_ref_ = table_.const_ref();
          }
          special_position_indices_.push_back(indices_.size());
          indices_.push_back(i);
        }
        indices_const_ref_ = indices_.const_ref();
      }

      //! Shorthand for: indices()[i_seq] != 0.
      /*! An exception is thrown if i_seq is out of bounds.
       */
      bool
      is_special_position(std::size_t i_seq) const
      {
        CCTBX_ASSERT(i_seq < indices_const_ref_.size());
        return indices_const_ref_[i_seq] != 0;
      }

      //! Shorthand for: table()[indices()[i_seq]]
      /*! An exception is thrown if i_seq is out of bounds.
       */
      site_symmetry_ops const&
      get(std::size_t i_seq) const
      {
        CCTBX_ASSERT(i_seq < indices_const_ref_.size());
        return table_const_ref_[indices_const_ref_[i_seq]];
      }

      //! Number of sites in special positions.
      std::size_t
      n_special_positions() const { return special_position_indices_.size(); }

      //! Indices of sites in special positions.
      /*! The indices refer to the order in which the site symmetries
          were passed to process().
       */
      af::shared<std::size_t> const&
      special_position_indices() const { return special_position_indices_; }

      //! Indices to internal table entries.
      af::shared<std::size_t> const&
      indices() const { return indices_; }

      //! Reference to indices to internal table entries; use for efficiency.
      /*! Not available in Python.
       */
      af::const_ref<std::size_t> const&
      indices_const_ref() const { return indices_const_ref_; }

      //! Number of internal table entries.
      std::size_t
      n_unique() const { return table_const_ref_.size(); }

      //! Internal table of unique site_symmetry_ops.
      /*! The table is organized such that table()[0].is_point_group_1()
          is always true after the first call of process().
       */
      af::shared<site_symmetry_ops> const&
      table() const { return table_; }

      /*! \brief Reference to internal table of unique site_symmetry_ops;
          use for efficiency.
       */
      /*! Not available in Python.
       */
      af::const_ref<site_symmetry_ops> const&
      table_const_ref() const { return table_const_ref_; }

      //! Pre-allocates memory for indices(); for efficiency.
      void
      reserve(std::size_t n_sites_final)
      {
        indices_.reserve(n_sites_final);
        indices_const_ref_ = indices_.const_ref();
      }

      //! Creates independent copy.
      site_symmetry_table
      deep_copy() const
      {
        site_symmetry_table result;
        result.indices_ = indices_.deep_copy();
        result.indices_const_ref_ = result.indices_.const_ref();
        result.table_ = table_.deep_copy();
        result.table_const_ref_ = result.table_.const_ref();
        result.special_position_indices_=special_position_indices_.deep_copy();
        return result;
      }

      //! Apply change-of-basis operator.
      site_symmetry_table
      change_basis(change_of_basis_op const& cb_op) const
      {
        site_symmetry_table result;
        result.indices_ = indices_.deep_copy();
        result.indices_const_ref_ = result.indices_.const_ref();
        result.table_.reserve(table_const_ref_.size());
        for(std::size_t i_seq=0;i_seq<table_const_ref_.size();i_seq++) {
          result.table_.push_back(table_const_ref_[i_seq].change_basis(cb_op));
        }
        result.table_const_ref_ = result.table_.const_ref();
        result.special_position_indices_=special_position_indices_.deep_copy();
        return result;
      }

      //! Apply selection.
      site_symmetry_table
      select(af::const_ref<std::size_t> const& indices) const
      {
        site_symmetry_table result;
        result.reserve(indices.size());
        for(std::size_t i=0;i<indices.size();i++) {
          result.process(get(indices[i]));
        }
        return result;
      }

      //! Support for Python's pickle facility. Do not use for other purposes.
      site_symmetry_table(
        af::shared<std::size_t> const& indices,
        af::shared<site_symmetry_ops> const& table,
        af::shared<std::size_t> const& special_position_indices)
      :
        indices_(indices),
        table_(table),
        special_position_indices_(special_position_indices)
      {
        indices_const_ref_ = indices_.const_ref();
        table_const_ref_ = table_.const_ref();
      }

    protected:
      af::shared<std::size_t> indices_;
      af::const_ref<std::size_t> indices_const_ref_;
      af::shared<site_symmetry_ops> table_;
      af::const_ref<site_symmetry_ops> table_const_ref_;
      af::shared<std::size_t> special_position_indices_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SITE_SYMMETRY_TABLE_H
