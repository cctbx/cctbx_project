#include <scitbx/serialization/base_256.h>
#include <cstddef>
#include <iostream>

using namespace scitbx;

namespace {

  static std::size_t n_delta_4 = 0;
  static std::size_t n_delta_8 = 0;

  template <typename T>
  inline
  void
  integer_core_exercise(T const& value)
  {
    char buf[128];
    char *end = serialization::base_256::to_string(buf, value);
    serialization::base_256::from_string<T> proxy(buf);
    SCITBX_ASSERT(proxy.end == end);
    SCITBX_ASSERT(proxy.value == value);
  }

  template <typename T>
  struct integer_max
  {
    static const int loop1 = 1000000;
    static const int loop2 = 2147481557-2357*2;
  };

  template <>
  struct integer_max<short>
  {
    static const int loop1 = 32767;
    static const int loop2 = 32767;
  };

  template <>
  struct integer_max<unsigned short> : integer_max<short>
  {};

  template <typename T>
  struct integer_signed
  {
    static void
    exercise()
    {
      for(T i=0;i<integer_max<T>::loop1;i++) {
        for (T s = 1; s >= -1; s -= 2) {
          integer_core_exercise(s*i);
        }
      }
      for(T i=0;i<integer_max<T>::loop2;i += 2357) {
        for (T s = 1; s >= -1; s -= 2) {
          integer_core_exercise(s*i);
        }
      }
    }
  };

  template <typename T>
  struct integer_unsigned
  {

    static void
    exercise()
    {
      for(T i=0;i<integer_max<T>::loop1;i++) integer_core_exercise(i);
      for(T i=0;i<integer_max<T>::loop2;i += 2357) integer_core_exercise(i);
    }
  };

  template <typename T>
  void
  floating_point_core_exercise(T const& value)
  {
    char buf[128];
    char *end = serialization::base_256::to_string(buf, value);
    serialization::base_256::from_string<T> proxy(buf);
    SCITBX_ASSERT(proxy.end == end);
    if (proxy.value != value) {
      T delta = proxy.value / value - 1;
      if (delta < 0) delta = -delta;
      if (sizeof(T) == 4) {
        n_delta_4++;
        SCITBX_ASSERT(delta < 1.e-5);
      }
      else {
        n_delta_8++;
        SCITBX_ASSERT(delta < 1.e-10);
      }
    }
  }

  template <typename T>
  struct floating_point
  {
    static void
    exercise()
    {
      floating_point_core_exercise(T(0));
      for(int i=1;i<10000;i++) {
        for (int s = 1; s >= -1; s -= 2) {
          floating_point_core_exercise(T(s*T(i)));
          floating_point_core_exercise(T(s/T(i)));
        }
      }
      for(int i=23571;i<2147481557-23571;i += 23571) {
        for (int s = 1; s >= -1; s -= 2) {
          floating_point_core_exercise(T(s*T(i)));
          floating_point_core_exercise(T(s/T(i)));
        }
      }
      T pi = 4 * std::atan(1.);
      T value = 1;
      for(int i=0;i<60;i++) {
        value *= pi;
        for (int s = 1; s >= -1; s -= 2) {
          floating_point_core_exercise(T(s*value));
          floating_point_core_exercise(T(s/value));
        }
      }
    }
  };
}

int main(int argc, char* /*argv*/[])
{
  try {
    for(;;)
    {
      integer_signed<short>::exercise();
      integer_unsigned<unsigned short>::exercise();
      integer_signed<int>::exercise();
      integer_unsigned<unsigned>::exercise();
      integer_signed<long>::exercise();
      integer_unsigned<unsigned long>::exercise();
      integer_signed<long long>::exercise();
      integer_unsigned<unsigned long long>::exercise();
      floating_point<float>::exercise();
      floating_point<double>::exercise();
      if (argc == 1) break;
    }
    if (n_delta_4 != 0) std::cout << "n_delta_4: " << n_delta_4 << std::endl;
    if (n_delta_8 != 0) std::cout << "n_delta_8: " << n_delta_8 << std::endl;
    std::cout << "OK" << std::endl;
  }
  catch (error const& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
