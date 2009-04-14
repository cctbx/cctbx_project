#include <iostream>
#include <scitbx/array_family/initialiser.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/mat_ref/make.h>

#include <scitbx/matrix/norms.h>


using scitbx::fn::approx_equal;
namespace af = scitbx::af;

typedef af::c_grid<2> dim;
typedef af::versa<double, dim> matrix_t;
typedef af::ref<double, dim> matrix_ref_t;

void exercise_norms() {
  using namespace scitbx::matrix;
  using scitbx::mat_ref_to;
  matrix_t a(dim(4,3));
  af::init(a) =  1,  2, -3,
                 4, -5,  6,
                 7,  8,  9,
                -1, -2, -3;
  SCITBX_ASSERT(norm_1(mat_ref_to(a)) == 21);
  SCITBX_ASSERT(norm_inf(mat_ref_to(a)) == 24);
  SCITBX_ASSERT(approx_equal(norm_frobenius(mat_ref_to(a)), std::sqrt(299.),
                             1e-12));
}

int main() {
  exercise_norms();
  std::cout << "OK\n";
  return 0;
}
