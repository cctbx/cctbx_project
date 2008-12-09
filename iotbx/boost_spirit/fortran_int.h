#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_INT_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_INT_H

#include <boost/spirit/core/parser.hpp>
#include <boost/spirit/core/scanner/scanner.hpp>
#include <boost/spirit/core/composite/actions.hpp>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit;

template <int Width>
struct fortran_int_parser : parser<fortran_int_parser<Width> >
{
  BOOST_STATIC_ASSERT(Width > 0);

  typedef fortran_int_parser<Width> self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, int>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t first = scan.first;
    int n = 0;
    for(; n < Width; ++n, ++scan) {
      if (scan.at_end()) return scan.no_match();
      if (*scan != ' ') break;
    }
    if (n == Width) return scan.no_match();
    int sign = +1;
    char ch = *scan;
    if (ch == '-') sign = -1;
    if (ch == '-' || ch == '+') {
        if (++n == Width) return scan.no_match();
        ++scan;
        if (scan.at_end()) return scan.no_match();
    }
    int result = 0;
    for (; n < Width; ++n, ++scan) {
      char ch = *scan;
      if (ch < '0' || ch > '9') return scan.no_match();
      result = 10*result + int(ch - '0');
    }
    if (sign == -1) result = -result;
    return scan.create_match(Width, result, first, first + Width);
  }
};


}} // iotbx::boost_spirit

#endif // GUARD
