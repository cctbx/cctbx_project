// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Copied from cctbx/loops.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_LOOPS_H
#define CCTBX_ARRAY_FAMILY_LOOPS_H

#include <algorithm>
#include <cstddef>
#include <cctbx/error.h>

#include <boost/config.hpp> // FIXES for broken compilers

namespace cctbx { namespace af {

  template <typename ArrayType>
  class nested_loop
  {
    public:
      nested_loop() : m_over(1) {}
      explicit
      nested_loop(const ArrayType& end,
                  bool open_range = true)
        : m_begin(end), m_end(end), m_current(end), m_over(0)
      {
        std::fill(m_begin.begin(), m_begin.end(), 0);
        m_current = m_begin;
        adjust_end(open_range);
      }
      nested_loop(const ArrayType& begin, const ArrayType& end,
                  bool open_range = true)
        : m_begin(begin), m_end(end), m_current(begin), m_over(0)
      {
        cctbx_assert(m_begin.size() == m_end.size());
        adjust_end(open_range);
      }
      bool incr()
      {
        for (std::size_t i = m_current.size(); i != 0;) {
          i--;
          m_current[i]++;
          if (m_current[i] < m_end[i]) return true;
          m_current[i] = m_begin[i];
        }
        m_over++;
        return false;
      }
      const ArrayType& begin() const { return m_begin; }
      const ArrayType& end() const { return m_end; }
      const ArrayType& operator()() const { return m_current; }
      std::size_t over() const { return m_over; }
    private:
      ArrayType m_begin;
      ArrayType m_end;
      ArrayType m_current;
      std::size_t m_over;

      void adjust_end(bool open_range)
      {
        if (!open_range) {
          for(std::size_t i=0;i<m_end.size();i++) {
            m_end[i] += 1;
          }
        }
      }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_LOOPS_H
