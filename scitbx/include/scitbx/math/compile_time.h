/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_MATH_COMPILE_TIME_H
#define SCITBX_MATH_COMPILE_TIME_H

#include <cstddef>

namespace scitbx { namespace math { namespace compile_time {

  template <std::size_t N>
  struct product
  {
    template <typename IndexType>
    static std::size_t
    get(IndexType const& i) { return i[N-1] * product<N-1>::get(i); }
  };

  template <>
  struct product<1>
  {
    template <typename IndexType>
    static std::size_t
    get(IndexType const& i) { return i[0]; }
  };

}}} // namespace scitbx::math::compile_time

#endif // SCITBX_MATH_COMPILE_TIME_H
