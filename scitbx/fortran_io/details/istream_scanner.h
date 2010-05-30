#ifndef SCITBX_FORTRAN_IO_DETAILS_ISTREAM_SCANNER_H
#define SCITBX_FORTRAN_IO_DETAILS_ISTREAM_SCANNER_H

#include <istream>

namespace scitbx { namespace fortran_io { namespace details {

template <typename CharType>
class istream_scanner
{
  public:
    typedef CharType value_t;
    typedef value_t &ref_t;

    bool at_end() const { return input.eof(); }

    value_t operator*() const { return current; }

    istream_scanner const &operator++() {
      current = input.get();
      return *this;
    }

    istream_scanner const &operator--() {
      input.putback(current);
      return *this;
    }

    istream_scanner(std::basic_istream<CharType> &input_stream)
      : input(input_stream) {
      ++(*this);
    }

  private:
    std::basic_istream<CharType> &input;
    value_t current;
};

}}} // scitbx::fortran_io::details

#endif // GUARD
