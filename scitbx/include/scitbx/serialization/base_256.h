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
          for(int i=0;i<sizeof(T);i++) {
            *u_buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
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
          int len = *u_end % 128;
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
            value *= 256;
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
          for(int i=0;i<sizeof(T);i++) {
            *u_buf++ = static_cast<unsigned char>(value % 256);
            value /= 256;
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
          int len = *u_end % 128;
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
            value *= 256;
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

  namespace floating_point {

    template <typename T>
    struct decomposition
    {
      decomposition(double x)
      {
        f = static_cast<T>(std::frexp(x, &e));
      }

      T f;
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
      decomposition<T> v(value);
      for(int i=0;i<sizeof(T);i++) {
        v.f = static_cast<T>(std::ldexp(v.f, 8)); // v.f *= 256;
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
        end(buf),
        value(0)
      {
        const unsigned char*
          u_end = reinterpret_cast<const unsigned char*>(end);
        int len = *u_end % 128;
        if (len == 0) {
          end++;
          return;
        }
        const unsigned char*
          u_buf = reinterpret_cast<const unsigned char*>(buf);
        u_end += len-1;
        for(;;) {
          value += static_cast<T>(*u_end--);
          value = static_cast<T>(std::ldexp(value, -8)); // value /= 256;
          if (u_end == u_buf) break;
        }
        if (*u_buf > 128) value = -value;
        base_256::from_string<int> e_proxy(buf + len);
        value = static_cast<T>(std::ldexp(value, e_proxy.value));
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
