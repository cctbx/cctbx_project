/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/error.h>
#include <cmath>

namespace scitbx { namespace serialization { namespace base_256 {

  namespace detail {

    double
    pow_10(double x, int e)
    {
      if (e < 0) {
        e = -e;
        while (e >= 4) { x /= 10000; e -= 4; }
        while (e >= 3) { x /= 1000; e -= 3; }
        while (e >= 2) { x /= 100; e -= 2; }
        if (e) x /= 10;
      }
      else {
        while (e >= 4) { x *= 10000; e -= 4; }
        while (e >= 3) { x *= 1000; e -= 3; }
        while (e >= 2) { x *= 100; e -= 2; }
        if (e) x *= 10;
      }
      return x;
    }

  } // namespace detail

  namespace integer {

    namespace signed_ {

      template <typename T>
      inline
      unsigned char*
      to_string(unsigned char* buf, T value)
      {
        unsigned char* buf0 = buf;
        *buf++ = 0;
        if (value != 0) {
          if (value < 0) {
            *buf0 = 128;
            value = -value;
          }
          for(int i=0;i<sizeof(T);i++) {
            *buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
            if (value == 0) break;
          }
          *buf0 += static_cast<unsigned char>(buf - buf0);
        }
        return buf;
      }

      template<typename T>
      struct from_string
      {
        from_string(const unsigned char* buf)
        :
          end(buf),
          value(0)
        {
          int len = *end % 128;
          if (len == 0) {
            end++;
            return;
          }
          end += len-1;
          for(;;) {
            value += static_cast<T>(*end--);
            if (end == buf) break;
            value *= 256;
          }
          if (*buf > 128) value *= -1;
          end = buf + len;
        }

        const unsigned char* end;
        T value;
      };

    } // namespace signed_

    namespace unsigned_ {

      template <typename T>
      inline
      unsigned char*
      to_string(unsigned char* buf, T value)
      {
        unsigned char* buf0 = buf;
        *buf++ = 0;
        if (value != 0) {
          for(int i=0;i<sizeof(T);i++) {
            *buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
            if (value == 0) break;
          }
          *buf0 += static_cast<unsigned char>(buf - buf0);
        }
        return buf;
      }

      template<typename T>
      struct from_string
      {
        from_string(const unsigned char* buf)
        :
          end(buf),
          value(0)
        {
          int len = *end;
          if (len == 0) {
            end++;
            return;
          }
          end += len-1;
          for(;;) {
            value += static_cast<T>(*end--);
            if (end == buf) break;
            value *= 256;
          }
          end = buf + len;
        }

        const unsigned char* end;
        T value;
      };

    } // namespace unsigned_

  } // namespace integer

  template <typename T>
  struct from_string;

  inline
  unsigned char*
  to_string(unsigned char* buf, int const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<int> : integer::signed_::from_string<int>
  {
    from_string(const unsigned char* buf)
    : integer::signed_::from_string<int>(buf)
    {}
  };

  inline
  unsigned char*
  to_string(unsigned char* buf, unsigned int const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned int>
  : integer::unsigned_::from_string<unsigned int>
  {
    from_string(const unsigned char* buf)
    : integer::unsigned_::from_string<unsigned int>(buf)
    {}
  };

  inline
  unsigned char*
  to_string(unsigned char* buf, long const& value)
  {
    return integer::signed_::to_string(buf, value);
  }

  template<>
  struct from_string<long> : integer::signed_::from_string<long>
  {
    from_string(const unsigned char* buf)
    : integer::signed_::from_string<long>(buf)
    {}
  };

  inline
  unsigned char*
  to_string(unsigned char* buf, unsigned long const& value)
  {
    return integer::unsigned_::to_string(buf, value);
  }

  template<>
  struct from_string<unsigned long>
  : integer::unsigned_::from_string<unsigned long>
  {
    from_string(const unsigned char* buf)
    : integer::unsigned_::from_string<unsigned long>(buf)
    {}
  };

  namespace floating_point {

    template <typename T>
    inline
    unsigned char*
    to_string(unsigned char* buf, T value)
    {
      unsigned char* buf0 = buf;
      *buf++ = 0;
      if (value == 0) return buf;
      if (value < 0) {
        *buf0 = 128;
        value = -value;
      }
      int e = static_cast<int>(std::log10(value)) + 1;
      T m = detail::pow_10(value, -e);
      for(int i=0;i<sizeof(T);i++) {
        m *= 256;
        int d = static_cast<int>(m);
        SCITBX_ASSERT(d < 256);
        *buf++ = static_cast<unsigned char>(d);
        m -= d;
        if (m == 0) break;
      }
      *buf0 += static_cast<unsigned char>(buf - buf0);
      return base_256::to_string(buf, e);
    }

    template<typename T>
    struct from_string
    {
      from_string(const unsigned char* buf)
      :
        end(buf),
        value(0)
      {
        int len = *end % 128;
        if (len == 0) {
          end++;
          return;
        }
        end += len-1;
        for(;;) {
          value += static_cast<T>(*end--);
          value /= 256;
          if (end == buf) break;
        }
        if (*buf > 128) value = -value;
        base_256::from_string<int> e_proxy(buf + len);
        value = detail::pow_10(value, e_proxy.value);
        end = e_proxy.end;
      }

      const unsigned char* end;
      T value;
    };

  } // namespace floating_point

  inline
  unsigned char*
  to_string(unsigned char* buf, float const& value)
  {
    return floating_point::to_string(buf, value);
  }

  template<>
  struct from_string<float> : floating_point::from_string<float>
  {
    from_string(const unsigned char* buf)
    : floating_point::from_string<float>(buf)
    {}
  };

  inline
  unsigned char*
  to_string(unsigned char* buf, double const& value)
  {
    return floating_point::to_string(buf, value);
  }

  template<>
  struct from_string<double> : floating_point::from_string<double>
  {
    from_string(const unsigned char* buf)
    : floating_point::from_string<double>(buf)
    {}
  };

}}} // namespace scitbx::serialization::base_256
