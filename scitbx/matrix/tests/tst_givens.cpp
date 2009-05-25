#include <scitbx/array_family/initialiser.h>
#include <scitbx/matrix/givens.h>
#include <scitbx/matrix/tests/utils.h>
#include <scitbx/error.h>
#include <iostream>

namespace af = scitbx::af;
using namespace scitbx::matrix;

void exercise_product() {
  givens::rotation<double> g;
  matrix_t a0(dim(3, 3));
  af::init(a0) = 1, 1, 1,
                2, 2, 2,
                3, 3, 3;
  {
    givens::product<double> p(2, true);
    double x, y;
    x =  1; y =  1; g.zero_x1(x, y);
    matrix_t g1(dim(3,3));
    af::init(g1) = g.c, -g.s, 0,
                   g.s,  g.c, 0,
                     0,    0, 1;
    p.multiply_by(g);
    x = -1; y =  1; g.zero_x1(x, y);
    matrix_t g2(dim(3,3));
    af::init(g2) = 1,   0,    0,
                   0, g.c, -g.s,
                   0, g.s,  g.c;
    p.multiply_by(g);
    matrix_t a = a0.deep_copy();
    p.apply_downward_on_right(a.ref(), 0);
    matrix_t a_g1 = af::matrix_multiply(a0.ref(), g1.ref());
    matrix_t a_g1_g2 = af::matrix_multiply(a_g1.ref(), g2.ref());
    SCITBX_ASSERT( a_g1_g2.all_approx_equal(a, 1e-15) );
    SCITBX_ASSERT( p.end == 0 );
  }

  {
    givens::product<double> p(2, false);
    double x, y;
    x =  1; y =  1; g.zero_x1(x, y);
    p.multiply_by(g);
    x = -1; y =  1; g.zero_x1(x, y);
    p.multiply_by(g);
    matrix_t a = a0.deep_copy();
    p.apply_downward_on_right(a.ref(), 0);
    SCITBX_ASSERT( a0.all_eq(a) );
    SCITBX_ASSERT( p.end == 0 );
  }
}


int main() {
  exercise_product();
  std::cout << "OK\n";
}
