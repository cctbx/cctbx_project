/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H
#define SCITBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H

#include <scitbx/error.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/small_reductions.h>

namespace scitbx { namespace af {

  typedef small<long, 10> flex_grid_default_index_type;

  template <typename IndexType = flex_grid_default_index_type>
  class flex_grid
  {
    public:
      typedef IndexType index_type;
      typedef typename index_type::value_type index_value_type;

      flex_grid() {}

      flex_grid(index_value_type const& grid_0)
      : origin_(1, index_value_type(0)),
        grid_(1, grid_0)
      {}

      flex_grid(index_type const& grid)
      : origin_(grid.size(), index_value_type(0)),
        grid_(grid)
      {}

      template <typename OtherArrayType>
      flex_grid(array_adaptor<OtherArrayType> const& a_a)
      : origin_((a_a.pointee)->size(), index_value_type(0)),
        grid_(a_a)
      {}

      flex_grid(index_type const& origin,
                index_type const& last,
                bool open_range = true)
      : origin_(origin)
      {
        SCITBX_ASSERT(origin_.size() == last.size());
        set_grid_(origin.size(), origin.begin(), last.begin(), open_range);
      }

      flex_grid
      set_layout(index_value_type const& layout_0)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        return *this;
      }

      flex_grid
      set_layout(index_type const& layout)
      {
        layout_ = layout;
        return *this;
      }

      std::size_t nd() const { return grid_.size(); }

      std::size_t size_1d() const
      {
        return af::product(grid_);
      }

      index_type const& origin() const { return origin_; }

      index_type const& grid() const { return grid_; }

      index_type last(bool open_range = true) const
      {
        index_value_type incl = 1;
        if (open_range) incl = 0;
        index_type result = origin_;
        for(std::size_t i=0;i<result.size();i++) {
          result[i] += grid_[i] - incl;
        }
        return result;
      }

      index_type const& layout() const { return layout_; }

      bool is_0_based() const
      {
        return origin_.all_eq(0);
      }

      bool is_padded() const
      {
        if (layout_.size() == 0) return false;
        SCITBX_ASSERT(grid_.size() == layout_.size());
        SCITBX_ASSERT(grid_.all_ge(layout_));
        return !grid_.all_eq(layout_);
      }

      std::size_t operator()(index_type const& i) const
      {
        std::size_t n = nd();
        std::size_t result = 0;
        if (n) {
          for(std::size_t j=0;;) {
            result += i[j] - origin_[j];
            j++;
            if (j == n) break;
            result *= grid_[j];
          }
        }
        return result;
      }

      bool is_valid_index(index_type const& i) const
      {
        std::size_t n = nd();
        if (i.size() != n) return false;
        for(std::size_t j=0;j<n;j++) {
          if (i[j] < origin_[j] || i[j] >= (origin_[j] + grid_[j])) {
            return false;
          }
        }
        return true;
      }

      bool operator==(flex_grid<index_type> const& other) const
      {
        if (!origin_.all_eq(other.origin_)) return false;
        if (!grid_.all_eq(other.grid_)) return false;
        return layout_.all_eq(other.layout_);
      }

      bool operator!=(flex_grid<index_type> const& other) const
      {
        return !(*this == other);
      }

    protected:
      index_type origin_;
      index_type grid_;
      index_type layout_;

      void set_grid_(
        std::size_t sz,
        const index_value_type* origin,
        const index_value_type* last,
        bool open_range)
      {
        index_value_type incl = 1;
        if (open_range) incl = 0;
        for(std::size_t i=0;i<sz;i++) {
          grid_.push_back(last[i] - origin[i] + incl);
        }
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H
