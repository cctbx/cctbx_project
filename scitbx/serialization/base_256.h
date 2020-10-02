#include <stdint.h>
#include <scitbx/error.h>
#include <cmath>

#if defined(__i386__) \
 || defined(__x86_64__) \
 || defined(__alpha__) \
 || defined(__host_mips) \
 || defined(__APPLE_CC__)
# define SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS
#endif

//#define SCITBX_SERIALIZATION_BASE_256_USE_LDEXP_POSITIVE
//#define SCITBX_SERIALIZATION_BASE_256_USE_LDEXP_NEGATIVE

namespace scitbx { namespace serialization { namespace base_256 {

  namespace integer {

    namespace signed_ {

      template <typename T>
      inline
      char*
      to_string(char* buf, T value)
      {
        unsigned char* u_buf = reinterpret_cast<unsigned char*>(buf);
        unsigned char* u_buf0 = u_buf;
        *u_buf++ = 0;
        if (value != 0) {
          if (value < 0) {
            *u_buf0 = 128;
            value = -value;
          }
          for(unsigned i=0;i<sizeof(T);i++) {
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
            *u_buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
#else
            *u_buf++ = static_cast<unsigned char>(value & 0xFFU);
            value >>= 8;
#endif
            if (value == 0) break;
          }
          *u_buf0 += static_cast<unsigned char>(u_buf - u_buf0);
        }
        return reinterpret_cast<char*>(u_buf);
      }

      template<typename T>
      struct from_string
      {
        from_string(const char* buf)
        :
          end(buf),
          value(0)
        {
          const unsigned char*
            u_end = reinterpret_cast<const unsigned char*>(end);
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
          int len = *u_end % 128;
#else
          int len = *u_end & 0x7FU;
#endif
          if (len == 0) {
            end++;
            return;
          }
          const unsigned char*
            u_buf = reinterpret_cast<const unsigned char*>(buf);
          u_end += len-1;
          for(;;) {
            value += static_cast<T>(*u_end--);
            if (u_end == u_buf) break;
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
            value *= 256;
#else
            value <<= 8;
#endif
          }
          if (*u_buf > 128) value = -value;
          end = buf + len;
        }

        const char* end;
        T value;
      };

    } // namespace signed_

    namespace unsigned_ {

      template <typename T>
      inline
      char*
      to_string(char* buf, T value)
      {
        unsigned char* u_buf = reinterpret_cast<unsigned char*>(buf);
        unsigned char* u_buf0 = u_buf;
        *u_buf++ = 0;
        if (value != 0) {
          for(unsigned i=0;i<sizeof(T);i++) {
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
            *u_buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
#else
            *u_buf++ = static_cast<unsigned char>(value & 0xFFU);
            value >>= 8;
#endif
            if (value == 0) break;
          }
          *u_buf0 += static_cast<unsigned char>(u_buf - u_buf0);
        }
        return reinterpret_cast<char*>(u_buf);
      }

      template<typename T>
      struct from_string
      {
        from_string(const char* buf)
        :
          end(buf),
          value(0)
        {
          const unsigned char*
            u_end = reinterpret_cast<const unsigned char*>(end);
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
          int len = *u_end % 128;
#else
          int len = *u_end & 0x7FU;
#endif
          if (len == 0) {
            end++;
            return;
          }
          const unsigned char*
            u_buf = reinterpret_cast<const unsigned char*>(buf);
          u_end += len-1;
          for(;;) {
            value += static_cast<T>(*u_end--);
            if (u_end == u_buf) break;
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
            value *= 256;
#else
            value <<= 8;
#endif
          }
          end = buf + len;
        }

        const char* end;
        T value;
      };

    } // namespace unsigned_

  } // namespace integer

  template <typename T>
  struct from_string;

  inline
  char*
  to_string(char* buf, short const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<short> : integer::signed_::from_string<short>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<short>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, int8_t const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<int8_t> : integer::signed_::from_string<int8_t>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<int8_t>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, uint8_t const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<uint8_t> : integer::signed_::from_string<uint8_t>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<uint8_t>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, unsigned short const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned short>
  : integer::unsigned_::from_string<unsigned short>
  {
    from_string(const char* buf)
    : integer::unsigned_::from_string<unsigned short>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, int const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<int> : integer::signed_::from_string<int>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<int>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, unsigned int const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned int>
  : integer::unsigned_::from_string<unsigned int>
  {
    from_string(const char* buf)
    : integer::unsigned_::from_string<unsigned int>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, long const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<long> : integer::signed_::from_string<long>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<long>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, unsigned long const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned long>
  : integer::unsigned_::from_string<unsigned long>
  {
    from_string(const char* buf)
    : integer::unsigned_::from_string<unsigned long>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, long long const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<long long> : integer::signed_::from_string<long long>
  {
    from_string(const char* buf)
    : integer::signed_::from_string<long long>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, unsigned long long const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned long long>
  : integer::unsigned_::from_string<unsigned long long>
  {
    from_string(const char* buf)
    : integer::unsigned_::from_string<unsigned long long>(buf)
    {}
  };

  namespace floating_point {

    struct decomposition
    {
      decomposition(double x)
      {
        f = std::frexp(x, &e);
      }

      double f;
      int e;
    };

    template <typename T>
    inline
    char*
    to_string(char* buf, T value)
    {
      unsigned char* u_buf = reinterpret_cast<unsigned char*>(buf);
      unsigned char* u_buf0 = u_buf;
      *u_buf++ = 0;
      if (value == 0) return reinterpret_cast<char*>(u_buf);
      if (value < 0) {
        *u_buf0 = 128;
        value = -value;
      }
      decomposition v(value);
      for(unsigned i=0;i<sizeof(double);i++) {
#ifndef SCITBX_SERIALIZATION_BASE_256_USE_LDEXP_POSITIVE
        v.f *= 256;
#else
        v.f = std::ldexp(v.f, 8);
#endif
        int d = static_cast<int>(v.f);
        SCITBX_ASSERT(d < 256);
        *u_buf++ = static_cast<unsigned char>(d);
        v.f -= d;
        if (v.f == 0) break;
      }
      *u_buf0 += static_cast<unsigned char>(u_buf - u_buf0);
      return base_256::to_string(reinterpret_cast<char*>(u_buf), v.e);
    }

    template<typename T>
    struct from_string
    {
      from_string(const char* buf)
      :
        end(buf)
      {
        const unsigned char*
          u_end = reinterpret_cast<const unsigned char*>(end);
#if !defined(SCITBX_SERIALIZATION_BASE_256_USE_BIT_SHIFTS)
        int len = *u_end % 128;
#else
        int len = *u_end & 0x7FU;
#endif
        if (len == 0) {
          value = 0;
          end++;
          return;
        }
        const unsigned char*
          u_buf = reinterpret_cast<const unsigned char*>(buf);
        u_end += len-1;
        double vd = 0;
        for(;;) {
          vd += static_cast<double>(*u_end--);
#ifndef SCITBX_SERIALIZATION_BASE_256_USE_LDEXP_NEGATIVE
          vd /= 256;
#else
          vd = std::ldexp(vd, -8);
#endif
          if (u_end == u_buf) break;
        }
        base_256::from_string<int> e_proxy(buf + len);
        value = static_cast<T>(std::ldexp(vd, e_proxy.value));
        if (*u_buf > 128) value = -value;
        end = e_proxy.end;
      }

      const char* end;
      T value;
    };

  } // namespace floating_point

  inline
  char*
  to_string(char* buf, float const& value)
  {
    return floating_point::to_string(buf, value);
  }

  template<>
  struct from_string<float> : floating_point::from_string<float>
  {
    from_string(const char* buf)
    : floating_point::from_string<float>(buf)
    {}
  };

  inline
  char*
  to_string(char* buf, double const& value)
  {
    return floating_point::to_string(buf, value);
  }

  template<>
  struct from_string<double> : floating_point::from_string<double>
  {
    from_string(const char* buf)
    : floating_point::from_string<double>(buf)
    {}
  };

}}} // namespace scitbx::serialization::base_256
