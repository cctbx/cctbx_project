#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_INT_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_INT_H

#include <boost/spirit/include/classic_parser.hpp>
#include <iotbx/boost_spirit/fortran_util.h>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit::classic;

template <typename IntegralType, int Width>
struct fortran_int_parser : parser<fortran_int_parser<IntegralType, Width> >
{
  BOOST_STATIC_ASSERT(Width > 0);

  typedef fortran_int_parser<IntegralType, Width> self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, IntegralType>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t first = scan.first;
    fortran_signed_int_extractor<IntegralType, Width> extractor(scan);
    if (extractor.failed) return scan.no_match();
    return scan.create_match(Width, extractor.value, first, first + Width);
  }
};


}} // iotbx::boost_spirit

#endif // GUARD
