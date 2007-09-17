#include <iostream>
#include <scitbx/error.h>
#include <scitbx/array_family/tiny_algebra.h>

using namespace scitbx;

void int_float_tiny_mult() {
  af::tiny<double, 3> b(1./4, 1./8, 1./32);
  af::tiny<std::size_t, 3> a(2, 4, 16);
  af::tiny<double, 3> c = a*b;
  SCITBX_ASSERT(c[0] == 0.5 && c[1] == 0.5 && c[2] == 0.5)(c[0])(c[1])(c[2]);
}


int main() {
  int_float_tiny_mult();
  std::cout << "OK\n";
}
