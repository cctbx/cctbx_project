#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_UTIL_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_UTIL_H

namespace iotbx { namespace boost_spirit {

template <typename IntegralType, int Width>
struct fortran_signed_int_extractor
{
  IntegralType value;
  bool negative;
  bool failed;

  template <typename ScannerType>
  fortran_signed_int_extractor(ScannerType const & scan)
    : failed(true)
  {
    int n = 0;
    for(; n < Width; ++n, ++scan) {
      if (scan.at_end()) return;
      if (*scan != ' ') break;
    }
    if (n == Width) return;
    negative = false;
    char ch = *scan;
    if (ch == '-') negative = true;
    if (ch == '-' || ch == '+') {
        if (++n == Width) return;
        ++scan;
        if (scan.at_end()) return;
    }
    value = 0;
    for (; n < Width; ++n, ++scan) {
      char ch = *scan;
      if (ch < '0' || ch > '9') return;
      value = 10*value + int(ch - '0');
    }
    failed = false;
    if (negative) value = -value;
  }
};


}} // iotbx::boost_spirit

#endif // GUARD
