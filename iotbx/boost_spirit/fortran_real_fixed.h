#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_REAL_FIXED_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_REAL_FIXED_H

#include <iotbx/boost_spirit/fortran_int.h>
#include <boost/spirit/include/classic_numerics.hpp>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit;

template <int Width, int FracDigits>
struct fortran_real_fixed : parser<fortran_real_fixed<Width, FracDigits> >
{
  BOOST_STATIC_ASSERT(Width > 0);
  static const int WholeDigits = Width - FracDigits - 1; // dot always printed
  BOOST_STATIC_ASSERT(WholeDigits > 0); // Never .123 style

  typedef fortran_real_fixed<Width, FracDigits> self_t;

  typedef fortran_int_parser<WholeDigits> whole_part_parser;
  typedef uint_parser<unsigned, /*Radix=*/10,
                      /*MinDigits=*/FracDigits, /*MaxDigits=*/FracDigits>
          frac_part_parser;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, double>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t first = scan.first;
    whole_part_parser whole_part_p;
    typename parser_result<self_t, ScannerType>::type
      whole_part_result = whole_part_p.parse(scan);
    if (!whole_part_result) return whole_part_result;
    double whole_part = double(whole_part_result.value());

    if (*scan != '.') return scan.no_match();
    ++scan;

    frac_part_parser frac_part_p;
    typename parser_result<self_t, ScannerType>::type
      frac_part_result = frac_part_p.parse(scan);
    if (!frac_part_result) return frac_part_result;
    double frac_part = double(frac_part_result.value())
                        * std::pow(10., -FracDigits);
    double res = whole_part >= 0 ? whole_part + frac_part
                                 : whole_part - frac_part;
    return scan.create_match(Width, res, first, first + Width);
  }
};

}} // iotbx::boost_spirit

#endif // GUARD
