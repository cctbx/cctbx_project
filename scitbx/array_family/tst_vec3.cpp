#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
    vec3<int> a;
    vec3<int> b(0,0,0);
    vec3<int> c(af::tiny_plain<int,3>(0,0,0));
    int id[] = {1,1,1};
    vec3<int> d(id);
  }
  {
    vec3<int> a(0,0,0);
    vec3<int> b(0,0,0);
    vec3<int> c(1,1,1);
    check_true(__LINE__, a == b);
    check_false(__LINE__, a == c);
    check_true(__LINE__, a == 0);
    check_false(__LINE__, a == 1);
    check_true(__LINE__, 0 == a);
    check_false(__LINE__, 1 == a);
    check_false(__LINE__, a != b);
    check_true(__LINE__, a != c);
    check_false(__LINE__, a != 0);
    check_true(__LINE__, a != 1);
    check_false(__LINE__, 0 != a);
    check_true(__LINE__, 1 != a);
  }
  {
    vec3<int> a(1,2,3);
    vec3<int> b(4,5,6);
    vec3<int> c(5,7,9);
    vec3<int> d(2,4,6);
    check_true(__LINE__, a + b == c);
    check_true(__LINE__, a + 3 == b);
    check_true(__LINE__, 3 + a == b);
    check_true(__LINE__, a - b == -3);
    check_true(__LINE__, b - 3 == a);
    check_true(__LINE__, c - b == a);
    check_true(__LINE__, a * b == 32);
    check_true(__LINE__, a * 2 == d);
    check_true(__LINE__, 2 * a == d);
    check_true(__LINE__, d / 2 == a);
    check_true(__LINE__, 12 / d == vec3<int>(6,3,2));
    check_true(__LINE__, a % 2 == vec3<int>(1,0,1));
    check_true(__LINE__, 4 % d == vec3<int>(0,0,4));
    {
      vec3<int> t(a);
      check_true(__LINE__, (t += b) == c);
    }
    {
      vec3<int> t(a);
      check_true(__LINE__, (t += 3) == b);
    }
    {
      vec3<int> t(a);
      check_true(__LINE__, (t -= b) == -3);
    }
    {
      vec3<int> t(b);
      check_true(__LINE__, (t -= 3) == a);
    }
    {
      vec3<int> t(a);
      check_true(__LINE__, (t *= 2) == d);
    }
    {
      vec3<int> t(d);
      check_true(__LINE__, (t /= 2) == a);
    }
    {
      vec3<int> t(a);
      check_true(__LINE__, (t %= 2) == vec3<int>(1,0,1));
    }
    check_true(__LINE__, vec3<int>(0,0,0).is_zero());
    check_false(__LINE__, vec3<int>(0,1,0).is_zero());
    check_true(__LINE__, -a == vec3<int>(-1,-2,-3));
    check_true(__LINE__, +a == a);
    check_true(__LINE__, a.cross(b) == vec3<int>(-3,6,-3));
  }
  {
    check_true(__LINE__, vec3<int>(-1,2,3).min() == -1);
    check_true(__LINE__, vec3<int>(1,-2,3).min() == -2);
    check_true(__LINE__, vec3<int>(1,2,-3).min() == -3);
    check_true(__LINE__, vec3<int>(1,-2,-3).max() == 1);
    check_true(__LINE__, vec3<int>(-1,2,-3).max() == 2);
    check_true(__LINE__, vec3<int>(-1,-2,3).max() == 3);
    check_true(__LINE__, vec3<int>(3,5,7).sum() == 3+5+7);
    check_true(__LINE__, vec3<int>(3,5,7).product() == 3*5*7);
  }
  {
    vec3<int> m(1,2,3);
    m.each_update_min(vec3<int>(0,3,3));
    check_true(__LINE__, m == vec3<int>(0,2,3));
    m.each_update_max(vec3<int>(0,3,4));
    check_true(__LINE__, m == vec3<int>(0,3,4));
  }
  {
    vec3<double> a(2., 3., 4.);
    check_true(__LINE__, std::fabs(a.length_sq() - 29.) < 1.e-6);
    check_true(__LINE__, std::fabs(a.length() - std::sqrt(29.)) < 1.e-6);
    check_true(__LINE__, std::fabs(abs(a) - std::sqrt(29.)) < 1.e-6);
    check_true(__LINE__, std::fabs(a.normalize().length() - 1.) < 1.e-6);
    check_true(__LINE__, std::fabs(
      vec3<double>(1.,.5,-1.8).angle(vec3<double>(-.3,.75,.5))
      - 1.99306755584) < 1.e-6);
    check_true(__LINE__, std::fabs(
      (  vec3<double>(1.,.5,-1.8).reflect(vec3<double>(1.,0.,1.))
       - vec3<double>(2.6,.5,-.2)).length()) < 1.e-6);
    check_true(__LINE__, std::fabs(
      (  vec3<double>(1.,-1.5,0.8).refract(vec3<double>(0.,1.,0.), 1.33)
       - vec3<double>(1.33,-1.79196,1.064)).length()) < 1.e-4);
  }
  {
    vec3<double> a(1., -2., 5.);
    check_true(__LINE__, std::fabs(a * a.ortho()) < 1.e-6);
  }
  {
    vec3<double> a(-2., 5., 1.);
    check_true(__LINE__, std::fabs(a * a.ortho()) < 1.e-6);
  }
  {
    vec3<double> a(5., 1., -2.);
    check_true(__LINE__, std::fabs(a * a.ortho()) < 1.e-6);
    verify(__LINE__, a, a.as_tiny());
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
