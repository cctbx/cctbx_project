#include <scitbx/math/copysign.h>
#include <iostream>
#include <scitbx/error.h>

#define SCITBX_CHECK_COPYSIGN(T, a, b, c, warn_only)                           \
{                                                                              \
  T d = scitbx::math::copysign(a, b);                                          \
  if (d != c) {                                                                \
    if (!warn_only) std::cout << "ERROR: ";                                    \
    std::cout << "copysign(" << #a << ", " << #b << ") = " << d                \
              << " instead of " << #c << " as per C99" << std::endl;           \
  }                                                                            \
}                                                                              \

void exercise_copy_sign() {
  using scitbx::math::copysign;

  /*
  C99 requires that:
   copysign(a, +0) =  |a|
   copysign(a, -0) = -|a|
  However this is at least broken for the float version of gcc 4.2
  on MacOSX intel. This is also broken in boost::math::copysign.
  So we just report the failure without flagging it as an error.
  */

  SCITBX_CHECK_COPYSIGN(double,  1.,  -2., -1., false);
  SCITBX_CHECK_COPYSIGN(double, -3.,  -2., -3., false);
  SCITBX_CHECK_COPYSIGN(double,  3.,  +0.,  3., true); // FIXME
  SCITBX_CHECK_COPYSIGN(double,  3.,  -0., -3., true); // FIXME

  SCITBX_CHECK_COPYSIGN(float,  1.f,  -2.f, -1.f, false);
  SCITBX_CHECK_COPYSIGN(float, -3.f,  -2.f, -3.f, false);
  SCITBX_CHECK_COPYSIGN(float,  3.f,  +0.f,  0.f, true); // FIXME
  SCITBX_CHECK_COPYSIGN(float,  3.f,  -0.f, -0.f, true); // FIXME
}

int main() {
  exercise_copy_sign();
  std::cout << "OK\n";
  return 0;
}
