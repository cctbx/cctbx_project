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
#include <cctbx/array_family/tiny_plain.h>
#include <cctbx/array_family/small_plain.h>
#include <cctbx/array_family/reductions.h>

namespace cctbx { namespace af {

  typedef small_plain<long, 10> flex_grid_default_index_type;

  template <typename IndexType = flex_grid_default_index_type>
  class flex_grid
  {
    public:
      typedef IndexType index_type;
      typedef typename IndexType::value_type value_type;

      flex_grid() {}

      explicit
      flex_grid(IndexType const& grid)
      : origin_(grid.size(), value_type(0)),
        grid_(grid)
      {}

      flex_grid(IndexType const& origin,
                IndexType const& last,
                bool open_range = true)
      : origin_(origin)
      {
        cctbx_assert(origin_.size() == last.size());
        set_grid_(origin.size(), origin.begin(), last.begin(), open_range);
      }

      template <typename value_type, std::size_t N>
      explicit
      flex_grid(
        tiny_plain<value_type, N> const& grid)
      : origin_(grid.size(), value_type(0))
      {
        grid_.insert(grid_.begin(), grid.begin(), grid.end());
      }

      template <typename value_type, std::size_t N>
      flex_grid(
        tiny_plain<value_type, N> const& origin,
        tiny_plain<value_type, N> const& last,
        bool open_range = true)
      {
        origin_.insert(origin_.begin(), origin.begin(), origin.end());
        set_grid_(N, origin.begin(), last.begin(), open_range);
      }

      std::size_t nd() const { return grid_.size(); }

      std::size_t size1d() const
      {
        return af::product(grid_.const_ref());
      }

      std::size_t operator()(IndexType const& i) const
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

      bool is_valid_index(const IndexType& i) const
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

    protected:
      index_type origin_;
      index_type grid_;

      void set_grid_(
        std::size_t sz,
        const value_type* origin,
        const value_type* last,
        bool open_range)
      {
        value_type incl = 1;
        if (open_range) incl = 0;
        for(std::size_t i=0;i<sz;i++) {
          grid_.push_back(last[i] - origin[i] + incl);
        }
      }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_FLEX_GRID_ACCESSOR_H
