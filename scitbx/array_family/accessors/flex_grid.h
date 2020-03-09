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
        if (origin_.all_eq(0)) origin_.clear();
      }

      flex_grid
      set_focus(index_type const& focus, bool open_range=true)
      {
        SCITBX_ASSERT(focus.size() == all_.size());
        focus_ = focus;
        if (!open_range && focus_.size() > 0) focus_ += index_value_type(1);
        set_focus_finalize();
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0)
      {
        SCITBX_ASSERT(all_.size() == 1);
        focus_.clear();
        focus_.push_back(focus_0);
        set_focus_finalize();
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1)
      {
        SCITBX_ASSERT(all_.size() == 2);
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        set_focus_finalize();
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2)
      {
        SCITBX_ASSERT(all_.size() == 3);
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        set_focus_finalize();
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2,
                index_value_type const& focus_3)
      {
        SCITBX_ASSERT(all_.size() == 4);
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        set_focus_finalize();
        return *this;
      }

      flex_grid
      set_focus(index_value_type const& focus_0,
                index_value_type const& focus_1,
                index_value_type const& focus_2,
                index_value_type const& focus_3,
                index_value_type const& focus_4)
      {
        SCITBX_ASSERT(all_.size() == 5);
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        focus_.push_back(focus_4);
        set_focus_finalize();
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
        SCITBX_ASSERT(all_.size() == 6);
        focus_.clear();
        focus_.push_back(focus_0);
        focus_.push_back(focus_1);
        focus_.push_back(focus_2);
        focus_.push_back(focus_3);
        focus_.push_back(focus_4);
        focus_.push_back(focus_5);
        set_focus_finalize();
        return *this;
      }

      std::size_t
      nd() const { return all_.size(); }

      std::size_t
      size_1d() const
      {
        SCITBX_ASSERT(all_.all_ge(0));
        return af::product(all_);
      }

      index_type const&
      all() const { return all_; }

      bool
      is_0_based() const { return origin_.size() == 0; }

      index_type
      origin() const
      {
        if (!is_0_based()) return origin_;
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
      is_padded() const { return focus_.size() != 0; }

      index_type
      focus(bool open_range=true) const
      {
        if (is_padded()) {
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
        SCITBX_ASSERT(n.all_ge(0));
        return af::product(n);
      }

      bool
      is_trivial_1d() const
      {
        if (nd() != 1) return false;
        if (!is_0_based()) return false;
        if (is_padded()) return false;
        return true;
      }

      bool
      is_square_matrix() const
      {
        if (nd() != 2) return false;
        if (all_[0] != all_[1]) return false;
        if (!is_0_based()) return false;
        if (is_padded()) return false;
        return true;
      }

      flex_grid
      shift_origin() const
      {
        if (is_0_based()) return *this;
        if (!is_padded()) return flex_grid(all_);
        index_type result_focus = focus_; // step-wise to avoid
        result_focus -= origin();         // gcc 2.96 internal error
        return flex_grid(all_).set_focus(result_focus);
      }

      bool
      is_valid_index(index_type const& i) const
      {
        std::size_t n = nd();
        if (i.size() != n) return false;
        if (!is_0_based()) {
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
          if (!is_0_based()) {
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

      // No assertions to enable maximum performance.
      std::size_t
      operator()(
        index_value_type const& i,
        index_value_type const& j) const
      {
        if (origin_.size() == 0) {
          return i * all_[1] + j;
        }
        return (i-origin_[0]) * all_[1] + (j-origin_[1]);
      }

      // No assertions to enable maximum performance.
      std::size_t
      operator()(
        index_value_type const& i,
        index_value_type const& j,
        index_value_type const& k) const
      {
        if (origin_.size() == 0) {
          return (i * all_[1] + j) * all_[2] + k;
        }
        return ((i-origin_[0]) * all_[1] + (j-origin_[1]))
             * all_[2] + (k-origin_[2]);
      }

      // No assertions to enable maximum performance.
      std::size_t
      operator()(
        index_value_type const& i,
        index_value_type const& j,
        index_value_type const& k,
        index_value_type const& l) const
      {
        if (origin_.size() == 0) {
          return ((i * all_[1] + j) * all_[2] + k) * all_[3] + l;
        }
        return (((i-origin_[0]) * all_[1] + (j-origin_[1]))
             * all_[2] + (k-origin_[2]))
             * all_[3] + (l-origin_[3]);
      }

      // No assertions to enable maximum performance.
      std::size_t
      operator()(
        index_value_type const& i,
        index_value_type const& j,
        index_value_type const& k,
        index_value_type const& l,
        index_value_type const& m) const
      {
        if (origin_.size() == 0) {
          return (((i * all_[1] + j) * all_[2] + k) * all_[3] + l) * all_[4] + m;
        }
        return ((((i-origin_[0]) * all_[1] + (j-origin_[1]))
             * all_[2] + (k-origin_[2]))
             * all_[3] + (l-origin_[3]))
             * all_[4] + (m-origin_[4]);
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

      void
      set_focus_finalize()
      {
        index_type last_ = last();
        if (focus_.all_eq(last_)) {
          focus_.clear();
        }
        else {
          SCITBX_ASSERT(last_.all_ge(focus_));
        }
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_FLEX_GRID_H
