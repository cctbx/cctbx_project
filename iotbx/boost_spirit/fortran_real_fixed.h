#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_REAL_FIXED_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_REAL_FIXED_H

#include <boost/spirit/include/classic_numerics.hpp>
#include <iotbx/boost_spirit/fortran_util.h>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit;

template <typename FloatingPointType, int Width, int FracDigits>
struct fortran_real_fixed_parser
  : parser<fortran_real_fixed_parser<FloatingPointType, Width, FracDigits> >
{
  BOOST_STATIC_ASSERT(Width > 0);
  static const int WholeDigits = Width - FracDigits - 1; // dot always printed
  BOOST_STATIC_ASSERT(WholeDigits > 0); // Never .123 style

  typedef fortran_real_fixed_parser<FloatingPointType, Width, FracDigits>
          self_t;

  typedef uint_parser<unsigned, /*Radix=*/10,
                      /*MinDigits=*/FracDigits, /*MaxDigits=*/FracDigits>
          frac_part_parser;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, FloatingPointType>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t first = scan.first;
    fortran_signed_int_extractor<int, WholeDigits> whole_part_extractor(scan);
    if (whole_part_extractor.failed) return scan.no_match();
    FloatingPointType whole_part = whole_part_extractor.value;
    bool negative = whole_part_extractor.negative;

    if (*scan != '.') return scan.no_match();
    ++scan;

    frac_part_parser frac_part_p;
    typename parser_result<frac_part_parser, ScannerType>::type
      frac_part_result = frac_part_p.parse(scan);
    if (!frac_part_result) return scan.no_match();
    FloatingPointType frac_part = FloatingPointType(frac_part_result.value())
                                    * std::pow(10., -FracDigits);
    FloatingPointType res = negative ? whole_part - frac_part
                                     : whole_part + frac_part;
    return scan.create_match(Width, res, first, first + Width);
  }
};

}} // iotbx::boost_spirit

#endif // GUARD
