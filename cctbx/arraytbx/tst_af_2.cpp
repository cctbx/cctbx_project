#include <cctbx/array_family/tiny.h>
#include <cctbx/array_family/tiny_algebra.h>
#include <cctbx/array_family/small.h>
#include <cctbx/array_family/small_algebra.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/shared_algebra.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/versa_algebra.h>

using namespace cctbx;

namespace {

}

int main(void)
{
  {
    af::tiny<int, 3> a1(1,2,3);
    af::tiny<int, 3> a2(4,5,6);
    a1 + a2;
    af::tiny<double, 3> a3(1.,2.,3.);
    a1 + a3;
  }
  {
    af::small<int, 3> a1;
    af::small<int, 3> a2;
    a1 + a2;
    af::small<double, 3> a3;
    a1 + a3;
  }
  {
    af::shared<int> a1;
    af::shared<int> a2;
    a1 + a2;
    af::shared<double> a3;
    a1 + a3;
  }
  {
    af::versa<int> a1;
    af::versa<int> a2;
    a1 + a2;
    af::versa<double> a3;
    a1 + a3;
  }
  return 0;
}
