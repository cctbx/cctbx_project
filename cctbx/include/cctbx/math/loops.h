// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Created (parts of cctbx/loops.h) (rwgk)
     Jan 2002: Created (parts of cctbx/sgtbx/utils.h) (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MATH_LOOPS_H
#define CCTBX_MATH_LOOPS_H

#include <algorithm>
#include <cstddef>
#include <cctbx/error.h>

#include <boost/config.hpp> // FIXES for broken compilers

namespace cctbx {

  template <std::size_t MaxN>
  class loop_n_from_m
  {
    public:
      loop_n_from_m() : m_m(0), m_n(0), m_over(1) {}
      loop_n_from_m(std::size_t m, std::size_t n) : m_m(m), m_n(n), m_over(0) {
        cctbx_assert(m_m >= m_n);
        cctbx_assert(MaxN >= m_n);
        for(std::size_t i=0;i<m_n;i++) m_current[i] = i;
      }
      bool incr()
      {
        if (m_n > 0) {
          std::size_t p, l;
          p = l = m_n - 1;
          for (;;) {
                m_current[p]++;
            if (m_current[p] + l == m_m + p) {
              if (p == 0) break;
              p--;
            }
            else if (p < l) {
              m_current[p + 1] = m_current[p];
                        p++;
            }
            else {
              return true;
            }
          }
        }
        m_over++;
        return false;
      }
      std::size_t m() { return m_m; }
      std::size_t n() { return m_n; }
      std::size_t operator[](std::size_t i) { return m_current[i]; }
      std::size_t over() const { return m_over; }
    private:
      std::size_t m_m;
      std::size_t m_n;
      std::size_t m_current[MaxN];
      std::size_t m_over;
  };

} // namespace cctbx

#endif // CCTBX_MATH_LOOPS_H
