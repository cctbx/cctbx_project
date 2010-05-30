#ifndef SCITBX_FORTRAN_IO_DETAILS_NUMERIC_EXTRACTORS_H
#define SCITBX_FORTRAN_IO_DETAILS_NUMERIC_EXTRACTORS_H

#include <cmath>

namespace scitbx { namespace fortran_io { namespace details {

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
    bool operator()(ScannerType &scan) {
      Heir &heir = static_cast<Heir &>(*this);
      for (; n < field_width; ++n, ++scan) {
        if (scan.at_end()) {
          if (shall_fail_at_end()) return false;
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
          if (shall_fail_at_end()) return false;
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
    bool operator()(ScannerType &scan) {
      fortran_decimal_extractor<FloatType> extract(
        strict_width, field_width, frac_digits);
      int exponent;
      if (!extract(scan)) {
        int n = extract.n;
        char ch = *scan;
        if (ch == 'e' || ch == 'd' || ch == 'E' || ch == 'D') {
          ++n, ++scan;
          if (n == field_width || scan.at_end()) return false;
        }
        fortran_int_extractor<int> extract_exponent(
          strict_width, field_width, n);
        if (!extract_exponent(scan)) return false;
        exponent = extract_exponent.value + extract.exponent;
      }
      else exponent = extract.exponent;
      value = FloatType(extract.mantissa)*std::pow(10., exponent);
      return true;
    }
  };

}}} // scitbx::fortran_io::details

#endif // GUARD
