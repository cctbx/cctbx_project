/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored part of sgtbx/matrix.cpp (rwgk)
 */

#ifndef CCTBX_RATIONAL_H
#define CCTBX_RATIONAL_H

#include <boost/rational.hpp>
#include <string>
#include <cstdio>

namespace cctbx {

  //! Formatting of rational numbers.
  template <typename IntType>
  std::string
  format(boost::rational<IntType> const& v, bool decimal=false)
  {
    if (v.numerator() == 0) return std::string("0");
    char buf[128];
    if (decimal) {
      std::sprintf(buf, "%.6g", double(v.numerator()) / v.denominator());
      char* cp = buf;
      if (*cp == '-') cp++;
      if (*cp == '0') {
        char* cpp = cp + 1; while (*cp) *cp++ = *cpp++;
      }
    }
    else if (v.denominator() == 1) {
      std::sprintf(buf, "%ld", static_cast<long>(v.numerator()));
    }
    else {
      std::sprintf(buf, "%ld/%ld", static_cast<long>(v.numerator()),
                                   static_cast<long>(v.denominator()));
    }
    return std::string(buf);
  }

} // namespace cctbx

#endif // CCTBX_RATIONAL_H
