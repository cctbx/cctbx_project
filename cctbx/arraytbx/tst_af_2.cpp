#include <cctbx/array_family/tiny.h>
#include <cctbx/array_family/tiny_algebra.h>

using namespace cctbx;

namespace {

}

int main(void)
{
  af::tiny<int, 3> a1;
  af::tiny<int, 3> a2;
  a1 + a2;
  af::tiny<double, 3> a3;
  a1 + a3;
  return 0;
}
