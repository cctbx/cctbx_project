/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_FLEX_GRID_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_FLEX_GRID_H

#include <scitbx/error.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/small_reductions.h>
#include <scitbx/array_family/small_algebra.h>

namespace scitbx { namespace af {

  typedef small<long, 10> flex_grid_default_index_type;

  template <typename IndexType = flex_grid_default_index_type>
  class flex_grid
  {
    public:
      typedef IndexType index_type;
      typedef typename index_type::value_type index_value_type;

      flex_grid() {}

      flex_grid(index_type const& grid)
      : grid_(grid)
      {}

      template <typename OtherArrayType>
      flex_grid(array_adaptor<OtherArrayType> const& a_a)
      : grid_(a_a)
      {}

      flex_grid(index_value_type const& grid_0)
      :
        grid_(1, grid_0)
      {}

      flex_grid(index_value_type const& grid_0,
                index_value_type const& grid_1)
      :
        grid_(1, grid_0)
      {
        grid_.push_back(grid_1);
      }

      flex_grid(index_value_type const& grid_0,
                index_value_type const& grid_1,
                index_value_type const& grid_2)
      :
        grid_(1, grid_0)
      {
        grid_.push_back(grid_1);
        grid_.push_back(grid_2);
      }

      flex_grid(index_value_type const& grid_0,
                index_value_type const& grid_1,
                index_value_type const& grid_2,
                index_value_type const& grid_3)
      :
        grid_(1, grid_0)
      {
        grid_.push_back(grid_1);
        grid_.push_back(grid_2);
        grid_.push_back(grid_3);
      }

      flex_grid(index_value_type const& grid_0,
                index_value_type const& grid_1,
                index_value_type const& grid_2,
                index_value_type const& grid_3,
                index_value_type const& grid_4)
      :
        grid_(1, grid_0)
      {
        grid_.push_back(grid_1);
        grid_.push_back(grid_2);
        grid_.push_back(grid_3);
        grid_.push_back(grid_4);
      }

      flex_grid(index_value_type const& grid_0,
                index_value_type const& grid_1,
                index_value_type const& grid_2,
                index_value_type const& grid_3,
                index_value_type const& grid_4,
                index_value_type const& grid_5)
      :
        grid_(1, grid_0)
      {
        grid_.push_back(grid_1);
        grid_.push_back(grid_2);
        grid_.push_back(grid_3);
        grid_.push_back(grid_4);
        grid_.push_back(grid_5);
      }

      flex_grid(index_type const& origin,
                index_type const& last,
                bool open_range=true)
      :
        grid_(last),
        origin_(origin)
      {
        grid_ -= origin_;
        if (!open_range) grid_ += index_value_type(1);
      }

      flex_grid
      set_layout(index_type const& layout, bool open_range=true)
      {
        layout_ = layout;
        if (!open_range && layout_.size() > 0) layout_ += index_value_type(1);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0,
                 index_value_type const& layout_1)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        layout_.push_back(layout_1);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0,
                 index_value_type const& layout_1,
                 index_value_type const& layout_2)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        layout_.push_back(layout_1);
        layout_.push_back(layout_2);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0,
                 index_value_type const& layout_1,
                 index_value_type const& layout_2,
                 index_value_type const& layout_3)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        layout_.push_back(layout_1);
        layout_.push_back(layout_2);
        layout_.push_back(layout_3);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0,
                 index_value_type const& layout_1,
                 index_value_type const& layout_2,
                 index_value_type const& layout_3,
                 index_value_type const& layout_4)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        layout_.push_back(layout_1);
        layout_.push_back(layout_2);
        layout_.push_back(layout_3);
        layout_.push_back(layout_4);
        return *this;
      }

      flex_grid
      set_layout(index_value_type const& layout_0,
                 index_value_type const& layout_1,
                 index_value_type const& layout_2,
                 index_value_type const& layout_3,
                 index_value_type const& layout_4,
                 index_value_type const& layout_5)
      {
        layout_.clear();
        layout_.push_back(layout_0);
        layout_.push_back(layout_1);
        layout_.push_back(layout_2);
        layout_.push_back(layout_3);
        layout_.push_back(layout_4);
        layout_.push_back(layout_5);
        return *this;
      }

      std::size_t
      nd() const { return grid_.size(); }

      std::size_t
      size_1d() const
      {
        return af::product(grid_);
      }

      index_type const&
      grid() const { return grid_; }

      bool
      has_origin() const { return origin_.size() != 0; }

      index_type
      origin() const
      {
        if (has_origin()) return origin_;
        return index_type(grid_.size(), index_value_type(0));
      }

      index_type
      last(bool open_range=true) const
      {
        index_type result = origin();
        result += grid_;
        if (!open_range) result -= index_value_type(1);
        return result;
      }

      bool
      has_layout() const { return layout_.size() != 0; }

      index_type
      layout(bool open_range=true) const
      {
        if (has_layout()) {
          index_type result = layout_;
          if (!open_range) result -= index_value_type(1);
          return result;
        }
        return last(open_range);
      }

      std::size_t
      layout_size_1d() const
      {
        if (layout_.size() == 0) return size_1d();
        index_type n = layout();
        n -= origin();
        return af::product(n);
      }

      bool
      is_0_based() const
      {
        return !has_origin() || origin_.all_eq(0);
      }

      bool
      is_padded() const
      {
        if (!has_layout()) return false;
        SCITBX_ASSERT(grid_.size() == layout_.size());
        SCITBX_ASSERT(last().all_ge(layout_));
        return !last().all_eq(layout_);
      }

      flex_grid
      shift_origin() const
      {
        if (is_0_based()) return *this;
        if (!has_layout()) return flex_grid(grid_);
        SCITBX_ASSERT(layout_.size() == grid_.size());
        index_type result_layout = layout_; // step-wise to avoid
        result_layout -= origin();          // gcc 2.96 internal error
        return flex_grid(grid_).set_layout(result_layout);
      }

      bool
      is_valid_index(index_type const& i) const
      {
        std::size_t n = nd();
        if (i.size() != n) return false;
        if (has_origin()) {
          for(std::size_t j=0;j<n;j++) {
            if (i[j] < origin_[j] || i[j] >= (origin_[j] + grid_[j])) {
              return false;
            }
          }
        }
        else {
          for(std::size_t j=0;j<n;j++) {
            if (i[j] < 0 || i[j] >= grid_[j]) {
              return false;
            }
          }
        }
        return true;
      }

      std::size_t
      operator()(index_type const& i) const
      {
        std::size_t n = nd();
        std::size_t result = 0;
        if (n) {
          if (has_origin()) {
            for(std::size_t j=0;;) {
              result += i[j] - origin_[j];
              j++;
              if (j == n) break;
              result *= grid_[j];
            }
          }
          else {
            for(std::size_t j=0;;) {
              result += i[j];
              j++;
              if (j == n) break;
              result *= grid_[j];
            }
          }
        }
        return result;
      }

      bool
      operator==(flex_grid<index_type> const& other) const
      {
        if (!grid_.all_eq(other.grid_)) return false;
        if (!origin_.all_eq(other.origin_)) return false;
        return layout_.all_eq(other.layout_);
      }

      bool
      operator!=(flex_grid<index_type> const& other) const
      {
        return !(*this == other);
      }

    protected:
      index_type grid_;
      index_type origin_;
      index_type layout_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_FLEX_GRID_H
