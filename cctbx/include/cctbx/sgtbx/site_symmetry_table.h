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
          indices_.push_back(i);
        }
        indices_const_ref_ = indices_.const_ref();
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

          Not available in Python.
       */
      af::const_ref<site_symmetry_ops> const&
      table() const { return table_const_ref_; }

      //! Pre-allocates memory for indices(); for efficiency.
      void
      reserve(std::size_t n_sites_final)
      {
        indices_.reserve(n_sites_final);
        indices_const_ref_ = indices_.const_ref();
      }

    protected:
      af::shared<std::size_t> indices_;
      af::const_ref<std::size_t> indices_const_ref_;
      af::shared<site_symmetry_ops> table_;
      af::const_ref<site_symmetry_ops> table_const_ref_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SITE_SYMMETRY_TABLE_H
