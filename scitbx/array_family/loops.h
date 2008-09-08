#ifndef SCITBX_ARRAY_FAMILY_LOOPS_H
#define SCITBX_ARRAY_FAMILY_LOOPS_H

#include <scitbx/error.h>
#include <algorithm>
#include <cstddef>

#include <boost/config.hpp> // FIXES for broken compilers

namespace scitbx { namespace af {

  template <typename ArrayType>
  class nested_loop
  {
    public:
      nested_loop() : over_(true) {}

      explicit
      nested_loop(typename ArrayType::value_type const& end)
      :
        over_(end == 0)
      {
        SCITBX_ASSERT(end >= 0);
        std::fill(begin_.begin(), begin_.end(), 0);
        std::fill(end_.begin(), end_.end(), end);
        current_ = begin_;
      }

      nested_loop(
        typename ArrayType::value_type const& begin,
        typename ArrayType::value_type const& end)
      :
        over_(end == begin)
      {
        SCITBX_ASSERT(end >= begin);
        std::fill(begin_.begin(), begin_.end(), begin);
        std::fill(end_.begin(), end_.end(), end);
        current_ = begin_;
      }

      explicit
      nested_loop(
        ArrayType const& end,
        bool open_range=true)
      :
        begin_(end), end_(end), current_(end), over_(true)
      {
        std::fill(begin_.begin(), begin_.end(), 0);
        current_ = begin_;
        adjust_end_and_over(open_range);
      }

      nested_loop(
        ArrayType const& begin,
        ArrayType const& end,
        bool open_range=true)
      :
        begin_(begin), end_(end), current_(begin), over_(true)
      {
        SCITBX_ASSERT(begin_.size() == end_.size());
        adjust_end_and_over(open_range);
      }

      bool
      incr()
      {
        for (std::size_t i = current_.size(); i != 0;) {
          i--;
          current_[i]++;
          if (current_[i] < end_[i]) return true;
          current_[i] = begin_[i];
        }
        over_ = true;
        return false;
      }

      ArrayType const&
      begin() const { return begin_; }

      ArrayType const&
      end() const { return end_; }

      ArrayType const&
      operator()() const { return current_; }

      bool over()
      const { return over_; }

    protected:
      ArrayType begin_;
      ArrayType end_;
      ArrayType current_;
      bool over_;

      void
      adjust_end_and_over(bool open_range)
      {
        if (!open_range) {
          for(std::size_t i=0;i<end_.size();i++) {
            end_[i] += 1;
          }
        }
        for(std::size_t i=0;i<end_.size();i++) {
          SCITBX_ASSERT(end_[i] >= begin_[i]);
          if (end_[i] > begin_[i]) over_ = false;
        }
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_LOOPS_H
