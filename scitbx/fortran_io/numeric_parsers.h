#ifndef SCITBX_FORTRAN_IO_NUMERIC_PARSERS_H
#define SCITBX_FORTRAN_IO_NUMERIC_PARSERS_H

#include <boost/spirit/include/classic_parser.hpp>
#include <cmath>
#include <scitbx/error.h>
#include <scitbx/fortran_io/details/numeric_extractors.h>

namespace scitbx { namespace fortran_io { namespace parsers {

using namespace boost::spirit::classic;

/// Parser following FORTRAN Iw editing of input fields
/** The width w is fixed at compile-time by a template parameter.
    Reference: FORTRAN standard, section 13.5.9.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <typename IntegralType, int FieldWidth, bool StrictWidth>
struct fortran_int_parser
  : parser<fortran_int_parser<IntegralType, FieldWidth, StrictWidth> >
{
  BOOST_STATIC_ASSERT(FieldWidth > 0);

  typedef fortran_int_parser<IntegralType, FieldWidth, StrictWidth>
          self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, IntegralType>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t start = scan.first;
    details::fortran_int_extractor<IntegralType> extract(
      StrictWidth, FieldWidth);
    if (!extract(scan)) return scan.no_match();
    return scan.create_match(scan.first - start, extract.value,
                             start, scan.first);
  }
};


/// Parser following FORTRAN Fw.d editing of input fields
/** The width w and number d of fractional digits are fixed
    at compile-time by template parameters.
    Reference: FORTRAN standard, section 13.5.9.2.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <
  typename FloatType, int FieldWidth, int FracDigits, bool StrictWidth>
struct fortran_real_parser
  : parser<fortran_real_parser<FloatType, FieldWidth, FracDigits, StrictWidth> >
{
  BOOST_STATIC_ASSERT(FieldWidth > 0);
  static const int MaxWholeDigits = FieldWidth - FracDigits - 1;
  BOOST_STATIC_ASSERT(MaxWholeDigits > 0);

  typedef fortran_real_parser<FloatType, FieldWidth, FracDigits, StrictWidth>
          self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, FloatType>::type type;
  };

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t start = scan.first;
    details::fortran_real_extractor<FloatType> extract(
      StrictWidth, FieldWidth, FracDigits);
    if (!extract(scan)) return scan.no_match();
    return scan.create_match(scan.first - start, extract.value,
                             start, scan.first);
  }
};

namespace dynamic {

/// Parser following FORTRAN Iw editing of input fields
/** The width w is fixed at run-time by a constructor parameter.
    Reference: FORTRAN standard, section 13.5.9.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <typename IntegralType>
struct fortran_int_parser
  : parser<fortran_int_parser<IntegralType> >
{
  typedef fortran_int_parser<IntegralType> self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, IntegralType>::type type;
  };

  bool const strict_width;
  int const field_width;

  fortran_int_parser(int field_width_, bool strict_width_)
    : strict_width(strict_width_), field_width(field_width_)
  {
    SCITBX_ASSERT(field_width > 0)(field_width);
  }

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t start = scan.first;
    details::fortran_int_extractor<IntegralType> extract(
      strict_width, field_width);
    if (!extract(scan)) return scan.no_match();
    return scan.create_match(scan.first - start, extract.value,
                             start, scan.first);
  }
};


/// Parser following FORTRAN Fw.d editing of input fields
/** The width w and number d of fractional digits are fixed
    at run-time by constructor parameters.
    Reference: FORTRAN standard, section 13.5.9.2.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <typename FloatType>
struct fortran_real_parser
  : parser<fortran_real_parser<FloatType> >
{
  typedef fortran_real_parser<FloatType> self_t;

  template <typename ScannerType>
  struct result
  {
    typedef typename match_result<ScannerType, FloatType>::type type;
  };

  bool const strict_width;
  int const field_width;
  int const frac_digits;

  fortran_real_parser(int field_width_, int frac_digits_, bool strict_width_)
    : strict_width(strict_width_), field_width(field_width_),
      frac_digits(frac_digits_)
  {
    SCITBX_ASSERT(field_width > 0)(field_width);
    int max_whole_digits = field_width - frac_digits - 1;
    SCITBX_ASSERT(max_whole_digits > 0)(max_whole_digits);
  }

  template <typename ScannerType>
  typename parser_result<self_t, ScannerType>::type
  parse(ScannerType const &scan) const {
    typename ScannerType::iterator_t start = scan.first;
    details::fortran_real_extractor<FloatType> extract(
      strict_width, field_width, frac_digits);
    if (!extract(scan)) return scan.no_match();
    return scan.create_match(scan.first - start, extract.value,
                             start, scan.first);
  }
};

}

}}} // scitbx::fortran_io::parsers

#endif // GUARD
