#include <cctbx/vec3.h>
#include <cctbx/array_family/simple_io.h>

using namespace cctbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
    vec3<int> va;
    vec3<int> vb(0,0,0);
    vec3<int> vc(af::tiny_plain<int,3>(0,0,0));
    int id[] = {1,1,1};
    vec3<int> vd(id);
  }
  {
    vec3<int> va(0,0,0);
    vec3<int> vb(0,0,0);
    vec3<int> vc(1,1,1);
    check_true(__LINE__, va == vb);
    check_false(__LINE__, va == vc);
    check_true(__LINE__, va == 0);
    check_false(__LINE__, va == 1);
    check_true(__LINE__, 0 == va);
    check_false(__LINE__, 1 == va);
    check_false(__LINE__, va != vb);
    check_true(__LINE__, va != vc);
    check_false(__LINE__, va != 0);
    check_true(__LINE__, va != 1);
    check_false(__LINE__, 0 != va);
    check_true(__LINE__, 1 != va);
  }
  {
    vec3<int> va(1,2,3);
    vec3<int> vb(4,5,6);
    vec3<int> vc(5,7,9);
    vec3<int> vd(2,4,6);
    check_true(__LINE__, va + vb == vc);
    check_true(__LINE__, va + 3 == vb);
    check_true(__LINE__, 3 + va == vb);
    check_true(__LINE__, va - vb == -3);
    check_true(__LINE__, vb - 3 == va);
    check_true(__LINE__, vc - vb == va);
    check_true(__LINE__, va * vb == 32);
    check_true(__LINE__, va * 2 == vd);
    check_true(__LINE__, 2 * va == vd);
    check_true(__LINE__, vd / 2 == va);
    check_true(__LINE__, 12 / vd == vec3<int>(6,3,2));
    check_true(__LINE__, va % 2 == vec3<int>(1,0,1));
    check_true(__LINE__, 4 % vd == vec3<int>(0,0,4));
    {
      vec3<int> t(va);
      check_true(__LINE__, (t += vb) == vc);
    }
    {
      vec3<int> t(va);
      check_true(__LINE__, (t += 3) == vb);
    }
    {
      vec3<int> t(va);
      check_true(__LINE__, (t -= vb) == -3);
    }
    {
      vec3<int> t(vb);
      check_true(__LINE__, (t -= 3) == va);
    }
    {
      vec3<int> t(va);
      check_true(__LINE__, (t *= 2) == vd);
    }
    {
      vec3<int> t(vd);
      check_true(__LINE__, (t /= 2) == va);
    }
    {
      vec3<int> t(va);
      check_true(__LINE__, (t %= 2) == vec3<int>(1,0,1));
    }
    check_true(__LINE__, -va == vec3<int>(-1,-2,-3));
    check_true(__LINE__, +va == va);
    check_true(__LINE__, va.cross(vb) == vec3<int>(-3,6,-3));
  }
  {
    vec3<double> va(2., 3., 4.);
    check_true(__LINE__, std::fabs(va.length() - std::sqrt(29.)) < 1.e-6);
    check_true(__LINE__, std::fabs(abs(va) - std::sqrt(29.)) < 1.e-6);
    check_true(__LINE__, std::fabs(va.normalize().length() - 1.) < 1.e-6);
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
    vec3<double> va(1., -2., 5.);
    check_true(__LINE__, std::fabs(va * va.ortho()) < 1.e-6);
  }
  {
    vec3<double> va(-2., 5., 1.);
    check_true(__LINE__, std::fabs(va * va.ortho()) < 1.e-6);
  }
  {
    vec3<double> va(5., 1., -2.);
    check_true(__LINE__, std::fabs(va * va.ortho()) < 1.e-6);
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
