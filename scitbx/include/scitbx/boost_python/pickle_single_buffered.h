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

#include <scitbx/serialization/base_256.h>
#include <complex>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  namespace base_256 = serialization::base_256;

  template <typename ValueType>
  struct from_string {};

  inline
  char*
  to_string(char* start, bool const& value)
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

    const char* end;
    bool value;
  };

  inline
  char*
  to_string(char* start, int const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<int> : base_256::from_string<int>
  {
    from_string(const char* start)
    : base_256::from_string<int>(start)
    {}
  };

  inline
  char*
  to_string(char* start, unsigned int const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<unsigned int> : base_256::from_string<unsigned int>
  {
    from_string(const char* start)
    : base_256::from_string<unsigned int>(start)
    {}
  };

  inline
  char*
  to_string(char* start, long const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<long> : base_256::from_string<long>
  {
    from_string(const char* start)
    : base_256::from_string<long>(start)
    {}
  };

  inline
  char*
  to_string(char* start, unsigned long const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<unsigned long> : base_256::from_string<unsigned long>
  {
    from_string(const char* start)
    : base_256::from_string<unsigned long>(start)
    {}
  };

  inline
  char*
  to_string(char* start, float const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<float> : base_256::from_string<float>
  {
    from_string(const char* start)
    : base_256::from_string<float>(start)
    {}
  };

  inline
  char*
  to_string(char* start, double const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<double> : base_256::from_string<double>
  {
    from_string(const char* start)
    : base_256::from_string<double>(start)
    {}
  };

  template <typename FloatType>
  inline
  char*
  to_string(char* start, std::complex<FloatType> const& value)
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

    const char* end;
    std::complex<double> value;
  };

  inline
  char*
  to_string(char* start, std::string const& value)
  {
    start = to_string(start, value.size());
    for(std::size_t i=0;i<value.size();i++) *start++ = value[i];
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
    }

    const char* end;
    std::string value;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#endif // SCITBX_BOOST_PYTHON_PICKLE_SINGLE_BUFFERED_H
