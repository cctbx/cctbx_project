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

      //! Insert new site_symmetry_ops into internal table.
      /*! The internal table is searched for identical site_symmetry_ops.
          The indices() for duplicate site_symmetry_ops will point to the
          same table() entry.
       */
      void
      process(
        std::size_t insert_at_index,
        site_symmetry_ops const& site_symmetry_ops_);

      //! Add new site_symmetry_ops to internal table.
      /*! See also: other overload.
       */
      void
      process(site_symmetry_ops const& site_symmetry_ops_)
      {
        process(indices_const_ref_.size(), site_symmetry_ops_);
      }

      /*! \brief Compute and process site symmetries for an array
          of original_sites_frac.
       */
      /*! See also: constructor of class site_symmetry
       */
      void
      process(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<scitbx::vec3<double> > const& original_sites_frac,
        af::const_ref<bool> const& unconditional_general_position_flags=
          af::const_ref<bool>(0,0),
        double min_distance_sym_equiv=0.5,
        bool assert_min_distance_sym_equiv=true)
      {
        CCTBX_ASSERT(
             unconditional_general_position_flags.size() == 0
          || unconditional_general_position_flags.size()
               == original_sites_frac.size());
        bool const* ugpf = unconditional_general_position_flags.begin();
        for(std::size_t i_seq=0;i_seq<original_sites_frac.size();i_seq++) {
          process(site_symmetry(
            unit_cell,
            space_group,
            original_sites_frac[i_seq],
            (ugpf != 0 && ugpf[i_seq] != 0 ? 0 : min_distance_sym_equiv),
            assert_min_distance_sym_equiv));
        }
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

      //! Number of sites on special positions.
      std::size_t
      n_special_positions() const { return special_position_indices_.size(); }

      //! Indices of sites on special positions.
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

      //! Support for asu_mappings::discard_last().
      /*! Not available in Python.
       */
      void
      discard_last()
      {
        if (indices_const_ref_.back() != 0) {
          special_position_indices_.pop_back();
        }
        indices_.pop_back();
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
      select(af::const_ref<std::size_t> const& selection) const
      {
        site_symmetry_table result;
        result.reserve(selection.size());
        for(std::size_t i=0;i<selection.size();i++) {
          result.process(get(selection[i]));
        }
        return result;
      }

      //! Apply selection.
      site_symmetry_table
      select(af::const_ref<bool> const& selection) const
      {
        CCTBX_ASSERT(selection.size() == indices_.size());
        site_symmetry_table result;
        for(std::size_t i_seq=0;i_seq<selection.size();i_seq++) {
          if (selection[i_seq]) {
            result.process(table_const_ref_[indices_const_ref_[i_seq]]);
          }
        }
        return result;
      }

      //! Mainly to support Python's pickle facility.
      /*! The inputs are NOT checked for consistency.
       */
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

  void
  replace(
    std::size_t replace_at_index,
    site_symmetry_ops const& site_symmetry_ops_)
  {
    CCTBX_ASSERT(indices_const_ref_.end() == indices_.end());
    CCTBX_ASSERT(table_const_ref_.end() == table_.end());
    CCTBX_ASSERT(replace_at_index <= indices_const_ref_.size());

    std::size_t i_tab = 0;

    bool new_is_pg_1 = site_symmetry_ops_.is_point_group_1();
    std::size_t old_i_tab=indices_[replace_at_index];
    bool old_is_pg_1 = table_[old_i_tab].is_point_group_1();
    indices_[replace_at_index]=0;

    bool old_sym_used=false;
    if (!old_is_pg_1)
      for(std::size_t itr=1;itr<indices_const_ref_.size();itr++) {
        if (indices_[itr] == old_i_tab) {
          old_sym_used=true;
          break;
        }
    }

        // Make changes in table
    if (!new_is_pg_1) {
      for(i_tab=1;i_tab<table_const_ref_.size();i_tab++) {
        if (table_const_ref_[i_tab] == site_symmetry_ops_) {
          break;
        }
      }
      if (i_tab == table_const_ref_.size()) {
        if (old_sym_used)
          table_.push_back(site_symmetry_ops_);
        else
          table_[old_i_tab]=site_symmetry_ops_;
        table_const_ref_ = table_.const_ref();
      }
    }

        // Make changes in indices and special_position_indices
    if (!new_is_pg_1 && old_sym_used){
        indices_[replace_at_index]=i_tab;
    }

    if (new_is_pg_1){
        if (!old_is_pg_1){
        // Remove from table_ too ?
                indices_[replace_at_index]=0;
                std::size_t n = special_position_indices_.size();
                std::size_t i = n-1;
                std::size_t* si = special_position_indices_.begin();
                while (i != 0) {
                  i--;
                  if (si[i] < replace_at_index) break;
                  si[i] = si[i+1];
                }
                special_position_indices_.resize(n-1);
        }
    }
    else if (old_is_pg_1){
      std::size_t n = special_position_indices_.size();
      special_position_indices_.resize(n+1);
      std::size_t* si = special_position_indices_.begin();
      std::size_t i = n;
      while (i != 0) {
        i--;
        if (si[i] < replace_at_index) break;
        si[i+1] = si[i];
      }
      si[i] = replace_at_index;
    }
    indices_const_ref_ = indices_.const_ref();
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
