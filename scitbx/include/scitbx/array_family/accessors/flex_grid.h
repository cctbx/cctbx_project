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

      flex_grid(index_type const& all)
      : all_(all)
      {}

      template <typename OtherArrayType>
      flex_grid(array_adaptor<OtherArrayType> const& a_a)
      : all_(a_a)
      {}

      flex_grid(index_value_type const& all_0)
      :
        all_(1, all_0)
      {}

      flex_grid(index_value_type const& all_0,
                index_value_type const& all_1)
      :
        all_(1, all_0)
      {
        all_.push_back(all_1);
      }

      flex_grid(index_value_type const& all_0,
                index_value_type const& all_1,
                index_value_type const& all_2)
      :
        all_(1, all_0)
      {
        all_.push_back(all_1);
        all_.push_back(all_2);
      }

      flex_grid(index_value_type const& all_0,
                index_value_type const& all_1,
                index_value_type const& all_2,
                index_value_type const& all_3)
      :
        all_(1, all_0)
      {
        all_.push_back(all_1);
        all_.push_back(all_2);
        all_.push_back(all_3);
      }

      flex_grid(index_value_type const& all_0,
                index_value_type const& all_1,
                index_value_type const& all_2,
                index_value_type const& all_3,
                index_value_type const& all_4)
      :
        all_(1, all_0)
      {
        all_.push_back(all_1);
        all_.push_back(all_2);
        all_.push_back(all_3);
        all_.push_back(all_4);
      }

      flex_grid(index_value_type const& all_0,
                index_value_type const& all_1,
                index_value_type const& all_2,
                index_value_type const& all_3,
                index_value_type const& all_4,
                index_value_type const& all_5)
      :
        all_(1, all_0)
      {
        all_.push_back(all_1);
        all_.push_back(all_2);
        all_.push_back(all_3);
        all_.push_back(all_4);
        all_.push_back(all_5);
      }

      flex_grid(index_type const& origin,
                index_type const& last,
                bool open_range=true)
      :
        all_(last),
        origin_(origin)
      {
        all_ -= origin_;
        if (!open_range) all_ += index_value_type(1);
      }

      flex_grid
      set_focus(index_type const& focus, bool open_range=true)
      {
        focus_ = focus;
        if (!open_range && focus_.size() > 0) focus_ += index_value_type(1);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2,
                index_value_type const& focus_3)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2,
                index_value_type const& focus_3,
                index_value_type const& focus_4)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        focus_.push_back(focus_4);
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2,
                index_value_type const& focus_3,
                index_value_type const& focus_4,
                index_value_type const& focus_5)
      {
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        focus_.push_back(focus_4);
        focus_.push_back(focus_5);
        return *this;
      }

      std::size_t
      nd() const { return all_.size(); }

      std::size_t
      size_1d() const
      {
        return af::product(all_);
      }

      index_type const&
      all() const { return all_; }

      bool
      has_origin() const { return origin_.size() != 0; }

      index_type
      origin() const
      {
        if (has_origin()) return origin_;
        return index_type(all_.size(), index_value_type(0));
      }

      index_type
      last(bool open_range=true) const
      {
        index_type result = origin();
        result += all_;
        if (!open_range) result -= index_value_type(1);
        return result;
      }

      bool
      has_focus() const { return focus_.size() != 0; }

      index_type
      focus(bool open_range=true) const
      {
        if (has_focus()) {
          index_type result = focus_;
          if (!open_range) result -= index_value_type(1);
          return result;
        }
        return last(open_range);
      }

      std::size_t
      focus_size_1d() const
      {
        if (focus_.size() == 0) return size_1d();
        index_type n = focus();
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
        if (!has_focus()) return false;
        SCITBX_ASSERT(all_.size() == focus_.size());
        SCITBX_ASSERT(last().all_ge(focus_));
        return !last().all_eq(focus_);
      }

      flex_grid
      shift_origin() const
      {
        if (is_0_based()) return *this;
        if (!has_focus()) return flex_grid(all_);
        SCITBX_ASSERT(focus_.size() == all_.size());
        index_type result_focus = focus_; // step-wise to avoid
        result_focus -= origin();         // gcc 2.96 internal error
        return flex_grid(all_).set_focus(result_focus);
      }

      bool
      is_valid_index(index_type const& i) const
      {
        std::size_t n = nd();
        if (i.size() != n) return false;
        if (has_origin()) {
          for(std::size_t j=0;j<n;j++) {
            if (i[j] < origin_[j] || i[j] >= (origin_[j] + all_[j])) {
              return false;
            }
          }
        }
        else {
          for(std::size_t j=0;j<n;j++) {
            if (i[j] < 0 || i[j] >= all_[j]) {
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
              result *= all_[j];
            }
          }
          else {
            for(std::size_t j=0;;) {
              result += i[j];
              j++;
              if (j == n) break;
              result *= all_[j];
            }
          }
        }
        return result;
      }

      bool
      operator==(flex_grid<index_type> const& other) const
      {
        if (!all_.all_eq(other.all_)) return false;
        if (!origin_.all_eq(other.origin_)) return false;
        return focus_.all_eq(other.focus_);
      }

      bool
      operator!=(flex_grid<index_type> const& other) const
      {
        return !(*this == other);
      }

    protected:
      index_type all_;
      index_type origin_;
      index_type focus_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_FLEX_GRID_H
