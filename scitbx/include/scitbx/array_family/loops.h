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
      nested_loop() : m_over(1) {}

      explicit
      nested_loop(typename ArrayType::value_type const& end)
      :
        m_over(0)
      {
        std::fill(begin_.begin(), begin_.end(), 0);
        std::fill(end_.begin(), end_.end(), end);
        current_ = begin_;
      }

      explicit
      nested_loop(
        ArrayType const& end,
        bool open_range=true)
      :
        begin_(end), end_(end), current_(end), m_over(0)
      {
        std::fill(begin_.begin(), begin_.end(), 0);
        current_ = begin_;
        adjust_end(open_range);
      }

      nested_loop(
        ArrayType const& begin,
        ArrayType const& end,
        bool open_range=true)
      :
        begin_(begin), end_(end), current_(begin), m_over(0)
      {
        SCITBX_ASSERT(begin_.size() == end_.size());
        adjust_end(open_range);
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
        m_over++;
        return false;
      }

      ArrayType const&
      begin() const { return begin_; }

      ArrayType const&
      end() const { return end_; }

      ArrayType const&
      operator()() const { return current_; }

      std::size_t over()
      const { return m_over; }

    protected:
      ArrayType begin_;
      ArrayType end_;
      ArrayType current_;
      std::size_t m_over;

      void
      adjust_end(bool open_range)
      {
        if (!open_range) {
          for(std::size_t i=0;i<end_.size();i++) {
            end_[i] += 1;
          }
        }
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_LOOPS_H
