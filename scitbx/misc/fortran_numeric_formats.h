#ifndef SCITBX_MISC_FORTRAN_NUMERIC_H
#define SCITBX_MISC_FORTRAN_NUMERIC_H

#include <boost/spirit/include/classic_parser.hpp>
#include <cmath>
#include <scitbx/error.h>

namespace scitbx { namespace boost_spirit_classic {

using namespace boost::spirit::classic;

namespace details {

  template <typename NumericType, class Heir>
  struct fortran_number_extractor
  {
    int n;
    bool negative;
    bool const strict_width;
    int const field_width;

    fortran_number_extractor(bool strict_width_, int field_width_,
                             int position_in_field=0)
      : strict_width(strict_width_), 
        field_width(field_width_),
        n(position_in_field), negative(false) 
    {}

    bool shall_fail_at_end() const {
      return strict_width && n-1 < field_width - 1;
    }

    template <typename ScannerType>
    bool operator()(ScannerType const &scan) {
      Heir &heir = static_cast<Heir &>(*this);
      for (; n < field_width; ++n, ++scan) {
        if (scan.at_end()) {
          if (shall_fail_at_end()) return scan.no_match();
          else break;
        }
        char ch = *scan;
        if (ch == ' ') continue;
        if (ch != '-' && ch != '+') break;
        if (ch == '-') negative = true;
        ++n; ++scan;
        break;
      }
      bool hit = true;
      for (; n < field_width; ++n, ++scan) {
        if (scan.at_end()) {
          if (shall_fail_at_end()) return scan.no_match();
          else break;
        }
        char ch = *scan;
        if (ch == ' ') continue;
        if ('0' <= ch && ch <= '9') {
          heir.digit(NumericType(ch - '0'));
        }
        else if (!heir.non_digit(ch)) {
          hit = false;
          break;
        }
      }
      heir.finish();
      return hit;
    }
  };

  template <typename IntegralType>
  struct fortran_int_extractor
    : fortran_number_extractor<
        IntegralType,
        fortran_int_extractor<IntegralType> >
  {
    typedef fortran_number_extractor<
              IntegralType,
              fortran_int_extractor<IntegralType> >
            base_t;

    IntegralType value;

    fortran_int_extractor(bool strict_width_, int field_width_,
                          int position_in_field=0)
      : value(0), base_t(strict_width_, field_width_, position_in_field)
    {}

    void digit(IntegralType d) { value = 10*value + d; }

    bool non_digit(char ch) { return false; }

    void finish() { if (this->negative) value = -value; }
  };

  template <typename FloatType>
  struct fortran_decimal_extractor
    : fortran_number_extractor<FloatType,
        fortran_decimal_extractor<FloatType> >
  {
    long mantissa;
    int exponent;
    int n_digits, n_whole_digits;
    
    typedef fortran_number_extractor<
              FloatType,
              fortran_decimal_extractor<FloatType> >
            base_t;

    fortran_decimal_extractor(bool strict_width_, int field_width_,
                              int frac_digits_)
      : mantissa(0), exponent(-frac_digits_), n_digits(0), n_whole_digits(-1),
        base_t(strict_width_, field_width_)
    {}

    void digit(FloatType d) {
      mantissa = 10*mantissa + d;
      ++n_digits;
    }

    bool non_digit(char ch) {
      bool result = false;
      if (ch == '.' && n_whole_digits == -1) {
        n_whole_digits = n_digits;
        result = true;
      }
      return result;
    }

    void finish() {
      if (this->negative) mantissa = - mantissa;
      if (n_whole_digits != -1) exponent = -(n_digits - n_whole_digits);
    }
  };
  
  template <typename FloatType>
  struct fortran_real_extractor
  {
    FloatType value;
    bool const strict_width;
    int const field_width;
    int const frac_digits;
    
    fortran_real_extractor(bool strict_width_, int field_width_,
                           int frac_digits_)
      : strict_width(strict_width_), field_width(field_width_),
        frac_digits(frac_digits_)
    {}
    
    template <typename ScannerType>
    bool operator()(ScannerType const &scan) {
      fortran_decimal_extractor<FloatType> extract(
        strict_width, field_width, frac_digits);
      int exponent;
      if (!extract(scan)) {
        int n = extract.n;
        char ch = *scan;
        if (ch == 'e' || ch == 'd') {
          ++n, ++scan;
          if (n == field_width || scan.at_end()) return false;
        }
        fortran_int_extractor<int> extract_exponent(
          strict_width, field_width, extract.n);
        if (!extract_exponent(scan)) return false;
        exponent = extract_exponent.value + extract.exponent;
      }
      else exponent = extract.exponent;
      value = FloatType(extract.mantissa)*std::pow(10., exponent);
      return true;
    }
  };
}


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


}} // scitbx::boost_spirit_classic

#endif // GUARD
