// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H
#define CCTBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H

#include <cctbx/error.h>
#include <cctbx/array_family/small_plain.h>
#include <cctbx/array_family/reductions.h>

namespace cctbx { namespace af {

  typedef small_plain<long, 10> flex_grid_default_index_type;

  template <typename IndexType = flex_grid_default_index_type>
  class flex_grid
  {
    public:
      typedef IndexType index_type;
      typedef typename index_type::value_type index_value_type;

      flex_grid() {}

      explicit
      flex_grid(index_type const& grid)
      : origin_(grid.size(), index_value_type(0)),
        grid_(grid)
      {}

      flex_grid(index_type const& origin,
                index_type const& last,
                bool open_range = true)
      : origin_(origin)
      {
        cctbx_assert(origin_.size() == last.size());
        set_grid_(origin.size(), origin.begin(), last.begin(), open_range);
      }

      flex_grid
      set_layout(index_type const& layout)
      {
        layout_ = layout;
        return *this;
      }

      std::size_t nd() const { return grid_.size(); }

      std::size_t size1d() const
      {
        return af::product(grid_.const_ref());
      }

      index_type const& origin() const { return origin_; }

      index_type const& grid() const { return grid_; }

      index_type const& layout() const { return layout_; }

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

      bool is_valid_index(const index_type& i) const
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
        if (!index_cmp_(origin_, other.origin_)) return false;
        if (!index_cmp_(grid_, other.grid_)) return false;
        return index_cmp_(layout_, other.layout_);
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

      // XXX use generic cmp after reorganization
      static
      bool
      index_cmp_(index_type const& lhs, index_type const& rhs)
      {
        if (lhs.size() != rhs.size()) return false;
        for(std::size_t i=0;i<lhs.size();i++) {
          if (lhs[i] != rhs[i]) return false;
        }
        return true;
      }
  };

  inline
  flex_grid<>
  make_flex_grid_1d(flex_grid<>::index_value_type const& n)
  {
    flex_grid_default_index_type grid;
    grid.push_back(n);
    return flex_grid<>(grid);
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H
