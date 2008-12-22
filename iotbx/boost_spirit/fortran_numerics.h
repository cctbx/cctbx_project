#ifndef IOTBX_BOOST_SPIRIT_FORTRAN_NUMERIC_H
#define IOTBX_BOOST_SPIRIT_FORTRAN_NUMERIC_H

#include <boost/spirit/include/classic_parser.hpp>
#include <cmath>

namespace iotbx { namespace boost_spirit {

using namespace boost::spirit::classic;

namespace details {

  template <typename NumericType, int FieldWidth, class Heir>
  struct fortran_number_extractor
  {
    int n;
    bool negative;

    fortran_number_extractor(int position_in_field=0)
      : n(position_in_field), negative(false) {}

    template <typename ScannerType>
    bool operator()(ScannerType const &scan) {
      Heir &heir = static_cast<Heir &>(*this);
      for (; n < FieldWidth; ++n, ++scan) {
        if (scan.at_end()) break;
        char ch = *scan;
        if (ch == ' ') continue;
        if (ch != '-' && ch != '+') break;
        if (ch == '-') negative = true;
        ++n; ++scan;
        break;
      }
      bool hit = true;
      for (; n < FieldWidth; ++n, ++scan) {
        if (scan.at_end()) break;
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

  template <typename IntegralType, int FieldWidth>
  struct fortran_int_extractor
    : fortran_number_extractor<
        IntegralType, FieldWidth,
        fortran_int_extractor<IntegralType, FieldWidth> >
  {
    typedef fortran_number_extractor<
              IntegralType, FieldWidth,
              fortran_int_extractor<IntegralType, FieldWidth> >
            base_t;

    IntegralType value;

    fortran_int_extractor(int position_in_field=0)
      : value(0), base_t(position_in_field)
    {}

    void digit(IntegralType d) { value = 10*value + d; }

    bool non_digit(char ch) { return false; }

    void finish() { if (this->negative) value = -value; }
  };

  template <typename FloatType, int FieldWidth, int FracDigits>
  struct fortran_decimal_extractor
    : fortran_number_extractor<
        FloatType, FieldWidth,
        fortran_decimal_extractor<FloatType, FieldWidth, FracDigits> >
  {
    long mantissa;
    int exponent;
    int n_digits, n_whole_digits;

    fortran_decimal_extractor()
      : mantissa(0), exponent(-FracDigits), n_digits(0), n_whole_digits(-1)
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
}


/// Parser following FORTRAN Iw editing of input fields
/** Reference: FORTRAN standard, section 13.5.9.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <typename IntegralType, int FieldWidth>
struct fortran_int_parser
  : parser<fortran_int_parser<IntegralType, FieldWidth> >
{
  BOOST_STATIC_ASSERT(FieldWidth > 0);

  typedef fortran_int_parser<IntegralType, FieldWidth>
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
    details::fortran_int_extractor<IntegralType, FieldWidth> extract;
    if (!extract(scan)) return scan.no_match();
    return scan.create_match(scan.first - start, extract.value,
                             start, scan.first);
  }
};


/// Parser following FORTRAN Fw.d editing of input fields
/** Reference: FORTRAN standard, section 13.5.9.2.1
    http://www.fortran.com/fortran/F77_std/rjcnf0001.html
*/
template <typename FloatType, int FieldWidth, int FracDigits>
struct fortran_real_parser
  : parser<fortran_real_parser<FloatType, FieldWidth, FracDigits> >
{
  BOOST_STATIC_ASSERT(FieldWidth > 0);
  static const int MaxWholeDigits = FieldWidth - FracDigits - 1;
  BOOST_STATIC_ASSERT(MaxWholeDigits > 0);

  typedef fortran_real_parser<FloatType, FieldWidth, FracDigits>
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
    details::fortran_decimal_extractor<
      FloatType, FieldWidth, FracDigits> extract;
    int exponent;
    if (!extract(scan)) {
      int n = extract.n;
      char ch = *scan;
      if (ch == 'e' || ch == 'd') {
        ++n, ++scan;
        if (n == FieldWidth || scan.at_end()) return scan.no_match();
      }
      details::fortran_int_extractor<
        int, FieldWidth> extract_exponent(extract.n);
      if (!extract_exponent(scan)) return scan.no_match();
      exponent = extract_exponent.value + extract.exponent;
    }
    else exponent = extract.exponent;
    FloatType res = FloatType(extract.mantissa)*std::pow(10., exponent);
    return scan.create_match(scan.first - start, res, start, scan.first);
  }
};


}} // iotbx::boost_spirit

#endif // GUARD
