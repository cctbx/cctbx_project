#include <scitbx/math/copysign.h>
#include <scitbx/math/imaginary.h>
#include <scitbx/math/utils.h>
#include <tbxx/pretty_type_name.hpp>
#include <iostream>
#include <iomanip>
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

void exercise_imaginary() {
  scitbx::math::imaginary_unit_t i;
  {
    typedef std::complex<double> c_t;

    SCITBX_ASSERT(2.*i == i*2.);
    SCITBX_ASSERT(2.*(2.*i) == 4.*i);
    SCITBX_ASSERT((2.*i)*2. == 4.*i);
    SCITBX_ASSERT(2./(2.*i) == -1.*i);
    SCITBX_ASSERT((2.*i)/2. == 1.*i);

    SCITBX_ASSERT(2.*i + 3.*i == 5.*i);
    SCITBX_ASSERT(2.*i - 3.*i == -1.*i);

    SCITBX_ASSERT(1.+2.*i == c_t(1., 2.));

    SCITBX_ASSERT(2.*i + c_t(3., -2.) == c_t(3., 0.));
    SCITBX_ASSERT(-2.*i - c_t(3., 2.) == c_t(-3., -4.));
    SCITBX_ASSERT(c_t(3., -2.) + 2.*i == c_t(3., 0.));
    SCITBX_ASSERT(c_t(3., 2.) - 2.*i == c_t(3., 0.));

    SCITBX_ASSERT(2.*i * c_t(1., 2.) == c_t(-4., 2.));
    SCITBX_ASSERT(c_t(1., 2.) * 2.*i == c_t(-4., 2.));

    double delta = std::abs( 5.*i / c_t(2., 1.) - c_t(1., 2.) );
    SCITBX_ASSERT(delta < 1e-15);
    SCITBX_ASSERT(c_t(3., 9.) / (3.*i) == c_t(3., -1.));
  }
}

template <typename UnsignedType>
struct unsigned_product_leads_to_overflow_test
{
  static void run() {
    using scitbx::math::unsigned_product_leads_to_overflow;
    {
      UnsignedType a[3];
      UnsignedType n = std::floor(std::pow(std::numeric_limits<UnsignedType>::max(), 1./3.));
      a[0] = n; a[1] = n; a[2] = n;
      a[0] = n+1; a[1] = n+1; a[2] = n+1;
      SCITBX_ASSERT(unsigned_product_leads_to_overflow(a, 3))
                   (tbxx::pretty_type_name<UnsignedType>())
                   (a[0])(a[1])(a[2]);
    }
    {
      UnsignedType a[3];
      UnsignedType n = std::numeric_limits<UnsignedType>::max();
      a[0] = n; a[1] = 1; a[2] = 1;
      SCITBX_ASSERT(!unsigned_product_leads_to_overflow(a, 3))
                   (tbxx::pretty_type_name<UnsignedType>())
                   (a[0])(a[1])(a[2]);
      a[0] = n; a[1] = 2; a[2] = 1;
      SCITBX_ASSERT(unsigned_product_leads_to_overflow(a, 3))
                   (tbxx::pretty_type_name<UnsignedType>())
                   (a[0])(a[1])(a[2]);
    }
  }
};

void exercise_unsigned_product_leads_to_overflow() {
  unsigned_product_leads_to_overflow_test<unsigned>::run();
  unsigned_product_leads_to_overflow_test<unsigned long>::run();
}

int main() {
  try {
    exercise_copy_sign();
    exercise_imaginary();
    exercise_unsigned_product_leads_to_overflow();
    std::cout << "OK\n";
    return 0;
  }
  catch(scitbx::error e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
