#ifndef SCITBX_RATIONAL_H
#define SCITBX_RATIONAL_H

#include <boost/rational.hpp>
#include <string>
#include <cstdio>

namespace scitbx {

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

  template <typename ArrayType>
  typename ArrayType::value_type
  array_lcm(ArrayType const& a)
  {
    typename ArrayType::value_type result;
    if (a.size() > 0) {
      result = a[0];
      for(std::size_t i=1;i<a.size();i++) {
        result = boost::lcm(result, a[i]);
      }
    }
    return result;
  }

} // namespace scitbx

#endif // SCITBX_RATIONAL_H
