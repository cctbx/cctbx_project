#ifndef SCITBX_FORTRAN_IO_NUMERIC_MANIPULATORS_H
#define SCITBX_FORTRAN_IO_NUMERIC_MANIPULATORS_H

#include <scitbx/fortran_io/details/istream_scanner.h>
#include <scitbx/fortran_io/details/numeric_extractors.h>

namespace scitbx { namespace fortran_io { namespace manipulators {

/// Stream manipulator following FORTRAN Iw editing
/** This modifies only the next input or output operation */
struct fortran_int {
  int width;
  bool strict_width;
  fortran_int(int w, bool strict=true) : width(w), strict_width(strict) {}
};

/// Stream manipulator following FORTRAN Iw editing
/** This modifies all subsequent input or output operations */
struct sticky_fortran_int : fortran_int  {
  sticky_fortran_int(fortran_int const &fmt) : fortran_int(fmt) {}
};

/// One of several overloaded functions to easily create sticky version
/// of a given manipulator
sticky_fortran_int sticky(fortran_int const &format) {
  return sticky_fortran_int(format);
}

/// Stream manipulator following FORTRAN Fw.d editing
/** This modifies only the next input or output operation */
struct fortran_real {
  int width, fractional_digits;
  bool strict_width;

  fortran_real(int w, int d, bool strict=true)
    : width(w), fractional_digits(d), strict_width(strict) {}
};

/// Stream manipulator following FORTRAN Fw.d editing
/** This modifies all subsequent input or output operations */
struct sticky_fortran_real : fortran_real {
  sticky_fortran_real(fortran_real const &fmt) : fortran_real(fmt) {}
};

/// One of several overloaded functions to easily create sticky version
/// of a given manipulator
sticky_fortran_real sticky(fortran_real const &format) {
  return sticky_fortran_real(format);
}

namespace details {

  template <class FormatType, class BoundFormatType>
  struct sticky_bound : BoundFormatType
  {
    typedef typename BoundFormatType::char_type char_type;
    sticky_bound(std::basic_istream<char_type> &input_stream,
                 FormatType const &format)
      : BoundFormatType(input_stream, format)
    {}

    template <typename IntegralType>
    sticky_bound &operator>>(IntegralType &i) {
      BoundFormatType::operator>>(i);
      return *this;
    }
  };

  template <typename CharType>
  struct input_bound_fortran_int
  {
    typedef CharType char_type;

    std::basic_istream<CharType> &input;
    fortran_int const &fmt;

    input_bound_fortran_int(std::basic_istream<char_type> &input_stream,
                            fortran_int const &format)
      : input(input_stream), fmt(format)
    {}

    template <typename IntegralType>
    std::basic_istream<char_type> &operator>>(IntegralType &i) {
      fortran_io::details::istream_scanner<char_type> scan(input);
      fortran_io::details::fortran_int_extractor<IntegralType> extract(
        fmt.strict_width, fmt.width);
      if (!extract(scan)) input.setstate(std::ios_base::badbit);
      else {
        i = extract.value;
        --scan;
      }
      return input;
    }
  };

  template <typename CharType>
  struct input_bound_fortran_real
  {
    typedef CharType char_type;

    std::basic_istream<CharType> &input;
    fortran_real const &fmt;

    input_bound_fortran_real(std::basic_istream<char_type> &input_stream,
                            fortran_real const &format)
      : input(input_stream), fmt(format)
    {}

    template <typename IntegralType>
    std::basic_istream<char_type> &operator>>(IntegralType &i) {
      fortran_io::details::istream_scanner<char_type> scan(input);
      fortran_io::details::fortran_real_extractor<IntegralType> extract(
        fmt.strict_width, fmt.width, fmt.fractional_digits);
      if (!extract(scan)) input.setstate(std::ios_base::badbit);
      else {
        i = extract.value;
        --scan;
      }
      return input;
    }
  };

} // scitbx::fortran_io::manipulators::details

template <typename CharType>
details::input_bound_fortran_int<CharType>
operator>>(std::basic_istream<CharType> &input,
            fortran_int const &format)
{
  return details::input_bound_fortran_int<CharType>(input, format);
}

template <typename CharType>
details::sticky_bound<sticky_fortran_int,
                      details::input_bound_fortran_int<CharType> >
operator>>(std::basic_istream<CharType> &input,
           sticky_fortran_int const &format)
{
  return details::sticky_bound<sticky_fortran_int,
                               details::input_bound_fortran_int<CharType> >(
    input, format);
}

template <typename CharType>
details::input_bound_fortran_real<CharType>
operator>>(std::basic_istream<CharType> &input,
            fortran_real const &format)
{
  return details::input_bound_fortran_real<CharType>(input, format);
}

template <typename CharType>
details::sticky_bound<sticky_fortran_real,
                      details::input_bound_fortran_real<CharType> >
operator>>(std::basic_istream<CharType> &input,
           sticky_fortran_real const &format)
{
  return details::sticky_bound<sticky_fortran_real,
                               details::input_bound_fortran_real<CharType> >(
    input, format);
}

}}} // scitbx::fortran_io::manipulators

#endif // GUARD
