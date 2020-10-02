#ifndef SCITBX_SERIALIZATION_SINGLE_BUFFERED_H
#define SCITBX_SERIALIZATION_SINGLE_BUFFERED_H

#include <stdint.h>
#include <scitbx/serialization/base_256.h>
#include <complex>

namespace scitbx { namespace serialization { namespace single_buffered {

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
  to_string(char* start, short const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<short> : base_256::from_string<short>
  {
    from_string(const char* start)
    : base_256::from_string<short>(start)
    {}
  };

  inline
  char*
  to_string(char* start, int8_t const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<int8_t> : base_256::from_string<int8_t>
  {
    from_string(const char* start)
    : base_256::from_string<int8_t>(start)
    {}
  };

  inline
  char*
  to_string(char* start, uint8_t const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<uint8_t> : base_256::from_string<uint8_t>
  {
    from_string(const char* start)
    : base_256::from_string<uint8_t>(start)
    {}
  };

  inline
  char*
  to_string(char* start, unsigned short const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<unsigned short> : base_256::from_string<unsigned short>
  {
    from_string(const char* start)
    : base_256::from_string<unsigned short>(start)
    {}
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
  to_string(char* start, long long const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<long long> : base_256::from_string<long long>
  {
    from_string(const char* start)
    : base_256::from_string<long long>(start)
    {}
  };

  inline
  char*
  to_string(char* start, unsigned long long const& value)
  {
    return base_256::to_string(start, value);
  }

  template <>
  struct from_string<unsigned long long>
  : base_256::from_string<unsigned long long>
  {
    from_string(const char* start)
    : base_256::from_string<unsigned long long>(start)
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

}}} // namespace scitbx::serialization::single_buffered

#endif // SCITBX_SERIALIZATION_SINGLE_BUFFERED_H
