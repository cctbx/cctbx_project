/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/arraytbx/flex_picklers.cpp (rwgk)
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2002 Aug: Created, based on shared_picklers.cpp (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_BOOST_PYTHON_PICKLE_SINGLE_BUFFERED_H
#define SCITBX_BOOST_PYTHON_PICKLE_SINGLE_BUFFERED_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <scitbx/error.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  namespace detail {

    inline char* o_advance(char *ptr)
    {
      while (*ptr != ',') ptr++;
      return ptr + 1;
    }

  }

  template <typename ValueType>
  struct from_string {};

  inline
  char* to_string(char* start, bool const& value)
  {
    if (value == true) *start = '1';
    else               *start = '0';
    return start + 1;
  }

  template <>
  struct from_string<bool>
  {
    from_string(const char* start) : end(start)
    {
      if (*end++ == '1') value = true;
      else               value = false;
    }

    bool value;
    const char* end;
  };

  inline
  char* to_string(char* start, int const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%d,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<int>
  {
    from_string(const char* start)
    {
      value = static_cast<int>(strtol(start, &end, 10));
      SCITBX_ASSERT(*end++ == ',');
    }

    int value;
    char* end;
  };

  inline
  char* to_string(char* start, unsigned int const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%u,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<unsigned int>
  {
    from_string(const char* start)
    {
      value = static_cast<unsigned int>(strtoul(start, &end, 10));
      SCITBX_ASSERT(*end++ == ',');
    }

    unsigned int value;
    char* end;
  };

  inline
  char* to_string(char* start, long const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%ld,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<long>
  {
    from_string(const char* start)
    {
      value = strtol(start, &end, 10);
      SCITBX_ASSERT(*end++ == ',');
    }

    long value;
    char* end;
  };

  inline
  char* to_string(char* start, unsigned long const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%lu,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<unsigned long>
  {
    from_string(const char* start)
    {
      value = strtoul(start, &end, 10);
      SCITBX_ASSERT(*end++ == ',');
    }

    unsigned long value;
    char* end;
  };

  inline
  char* to_string(char* start, float const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%.6g,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<float>
  {
    from_string(const char* start)
    {
      value = strtod(start, &end);
      SCITBX_ASSERT(*end++ == ',');
    }

    float value;
    char* end;
  };

  inline
  char* to_string(char* start, double const& value)
  {
    SCITBX_ASSERT(sprintf(start, "%.12g,", value) > 0);
    return detail::o_advance(start);
  }

  template <>
  struct from_string<double>
  {
    from_string(const char* start)
    {
      value = strtod(start, &end);
      SCITBX_ASSERT(*end++ == ',');
    }

    double value;
    char* end;
  };

  template <typename FloatType>
  inline
  char* to_string(char* start, std::complex<FloatType> const& value)
  {
    return to_string(to_string(start, value.real()), value.imag());
  }

  template <>
  struct from_string<std::complex<double> >
  {
    from_string(const char* start)
    {
      from_string<double> proxy_r(start);
      from_string<double> proxy_i(proxy_r.end);
      value = std::complex<double>(proxy_r.value, proxy_i.value);
      end = proxy_i.end;
    }

    std::complex<double> value;
    char* end;
  };

  inline
  char* to_string(char* start, std::string const& value)
  {
    start = to_string(start, value.size());
    for(std::size_t i=0;i<value.size();i++) *start++ = value[i];
    *start++ = ',';
    return start;
  }

  template <>
  struct from_string<std::string>
  {
    from_string(const char* start)
    {
      from_string<std::size_t> proxy_len(start);
      value.append(proxy_len.end, proxy_len.value);
      end = proxy_len.end + proxy_len.value;
      SCITBX_ASSERT(*end++ == ',');
    }

    std::string value;
    char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#endif // SCITBX_BOOST_PYTHON_PICKLE_SINGLE_BUFFERED_H
