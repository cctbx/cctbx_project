#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_INT_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_INT_H

#include <boost/spirit/include/classic_parser.hpp>
#include <boost/spirit/include/classic_scanner.hpp>
#include <boost/spirit/include/classic_actions.hpp>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit::classic;

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
    int res = 0;
    for (; n < Width; ++n, ++scan) {
      char ch = *scan;
      if (ch < '0' || ch > '9') return scan.no_match();
      res = 10*res + int(ch - '0');
    }
    if (sign == -1) res = -res;
    return scan.create_match(Width, res, first, first + Width);
  }
};


}} // iotbx::boost_spirit

#endif // GUARD
