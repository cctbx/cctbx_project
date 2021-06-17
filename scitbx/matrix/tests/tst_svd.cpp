#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/initialiser.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/matrix/special_matrices.h>
#include <scitbx/error.h>
#include <scitbx/constants.h>
#include <iostream>

#include <scitbx/matrix/svd.h>
#include <scitbx/matrix/tests.h>
#include <scitbx/matrix/tests/utils.h>

namespace af = scitbx::af;
using scitbx::fn::approx_equal;
using namespace scitbx::matrix;

template<class SVD2x2Type>
void check_2x2_decomposition(double f, double g, double h,
                             SVD2x2Type const &svd,
                             double tol)
{
  double cu=svd.c_u, su=svd.s_u, cv=svd.c_v, sv=svd.s_v;
  SCITBX_ASSERT(approx_equal(cu*cu + su*su, 1., tol));
  SCITBX_ASSERT(approx_equal(cv*cv + sv*sv, 1., tol));

  matrix_t u(dim(2,2));
  af::init(u) =  cu, -su,
                 su,  cu;
  matrix_t v(dim(2,2));
  af::init(v) =  cv, -sv,
                 sv,  cv;
  matrix_t a(dim(2,2));
  af::init(a) = f, g,
                0, h;
  vec_t sigma(2);
  af::init(sigma) = svd.s_max, svd.s_min;
  matrix_t a1 = svd::reconstruct(u.ref(), v.ref(), sigma.ref());
  SCITBX_ASSERT( std::abs((a1(0,0) -f)/f) < tol );
  SCITBX_ASSERT( std::abs((a1(0,1) -g)/g) < tol );
  SCITBX_ASSERT( std::abs(a1(1,0)) < tol );
  SCITBX_ASSERT( std::abs((a1(1,1) -h)/h) < tol );
}


void exercise_2x2_decomposition() {
  double tol = 5*std::numeric_limits<double>::epsilon();

  double s_min = std::sqrt(3. - std::sqrt(5.));
  double s_max = std::sqrt(3. + std::sqrt(5.));
  double sign[2] = { -1., 1. };
  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      for (int k=0; k<2; ++k) {
        double f = scitbx::math::copysign(2., sign[i]);
        double g = scitbx::math::copysign(1., sign[j]);
        double h = scitbx::math::copysign(1., sign[k]);
        svd::bidiagonal_2x2_decomposition<double> svd(f, g, h, true);
        check_2x2_decomposition(f, g, h, svd, tol);
        SCITBX_ASSERT( approx_equal(std::abs(svd.s_min), s_min, tol) );
        SCITBX_ASSERT( approx_equal(std::abs(svd.s_max), s_max, tol) );
  }

  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      for (int k=0; k<2; ++k) {
        double f = scitbx::math::copysign(1.,   sign[i]);
        double g = scitbx::math::copysign(10.,  sign[j]);
        double h = scitbx::math::copysign(100., sign[k]);
        double elt[18] = { f, g, h,
                           f, h, g,
                           g, f, h,
                           g, h, f,
                           h, f, g,
                           h, g, f  };
        for (int p=0; p<18; p+=3) {
          svd::bidiagonal_2x2_decomposition<double> svd(
            elt[p], elt[p+1], elt[p+2], true);
          check_2x2_decomposition(elt[p], elt[p+1], elt[p+2], svd, tol);
        }
  }
}

void exercise_lower_bidiagonal_rectification() {
  int n = 5;
  af::shared<double> d(n), f(n-1);
  af::init(d) = -2, 5, 8, -3, 7;
  af::init(f) = -1, 4, 5, 2;
  matrix_t a = lower_bidiagonal(d.ref(), f.ref());
  matrix_const_ref_t a_ = a.const_ref();
  matrix_t identity_n = identity<double>(d.size());
  matrix_ref_t u_ = identity_n.ref(), v_;
  svd::bidiagonal_decomposition<double> svd(
        d.ref(), f.ref(), svd::lower_bidiagonal_kind, u_, true, v_, false);
  matrix_t b = upper_bidiagonal(d.ref(), f.ref());
  matrix_const_ref_t b_ = b.const_ref();
  matrix_t u_b = af::matrix_multiply(u_, b_);
  matrix_const_ref_t u_b_ = u_b.const_ref();
  SCITBX_ASSERT( normality_ratio(u_) < 10 );
  SCITBX_ASSERT( equality_ratio(a_, u_b_) < 10 );
}

struct golub_kahan_iteration_test_case
{
  af::shared<double> diagonal, superdiagonal;
  af::shared<double> target_diagonal, target_superdiagonal;
  int r, s;
  typedef svd::bidiagonal_decomposition<double> svd_t;
  typedef void (svd_t::*tested_member_func_t)(bool);
  tested_member_func_t tested;

  golub_kahan_iteration_test_case(int n, int r_, int s_,
                                  tested_member_func_t tested_)
    : diagonal(n), superdiagonal(n-1),
      target_diagonal(n), target_superdiagonal(n-1),
      r(r_), s(s_), tested(tested_)
  {}

  void exercise() {
    matrix_t a = upper_bidiagonal(diagonal.ref(), superdiagonal.ref());
    matrix_ref_t a_ = a.ref();

    matrix_t identity_n = identity<double>(diagonal.size());

    matrix_t u = identity_n.deep_copy();
    matrix_ref_t u_ = u.ref();

    matrix_t v = identity_n.deep_copy();
    matrix_ref_t v_ = v.ref();

    svd::bidiagonal_decomposition<double> svd(
      diagonal.ref(), superdiagonal.ref(), svd::upper_bidiagonal_kind,
      u_, true, v_, true);

    svd.r = r; svd.s = s;
    (svd.*tested)(/*compute_shift=*/true);

    SCITBX_ASSERT( normality_ratio(u_) < 10 );
    SCITBX_ASSERT( normality_ratio(v_) < 10 );

    matrix_t ut_a_v = product_UT_M_V(u_, a_, v_);
    matrix_t a1 = upper_bidiagonal(diagonal.ref(), superdiagonal.ref());
    SCITBX_ASSERT( equality_ratio(ut_a_v.ref(), a1.ref()) < 10 );
  }
};

struct golub_kahan_iteration_test_case_1 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_1()
    : golub_kahan_iteration_test_case(
                4+3, 2, 6, &svd_t::do_implicit_shift_qr_iteration_downward)
  {
    af::init(diagonal) = 222, 222, 1, 2, 3, 4, 222;
    af::init(superdiagonal) = 333, 0, 1, 1, 1, 0;
    // Carefully obtain by hand using Mathematica (c.f. tst_svd.nb)
    af::init(target_diagonal) = 222, 222,
                                0.894427190999916, -2.1525179054617247,
                                3.8902675436754994, -3.204350435852646,
                                222;
    af::init(target_superdiagonal) = 333, 0,
                                     0.5477225575051663, 1.0963204208751276,
                                     0.8140669040777953,
                                     0;
  }
};

struct golub_kahan_iteration_test_case_2 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_2()
    : golub_kahan_iteration_test_case(
                3, 0, 3, &svd_t::do_implicit_shift_qr_iteration_downward)
  {
    af::init(diagonal) = 1, 2, 3;
    af::init(superdiagonal) = 2, 4;
    // Carefully obtain by hand using Mathematica (c.f. tst_svd.nb)
    af::init(target_diagonal) = 2.694486452681983, 3.625541410822385,
                                0.6141894835903716;
    af::init(target_superdiagonal) = 3.1321583436450946, -1.8459543936150495;
  }
};

struct golub_kahan_iteration_test_case_3 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_3()
    : golub_kahan_iteration_test_case(
                3, 0, 3, &svd_t::do_implicit_shift_qr_iteration_downward)
  {
    af::init(diagonal) = 1, 0, 3;
    af::init(superdiagonal) = 1, 1;
    // Carefully obtain by hand using Mathematica (c.f. tst_svd.nb)
    af::init(target_diagonal) = std::sqrt(2.), 0., 3.;
    af::init(target_superdiagonal) = 0., 1.;
  }
};

struct golub_kahan_iteration_test_case_4 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_4()
    : golub_kahan_iteration_test_case(
                3+3, 2, 5, &svd_t::do_implicit_shift_qr_iteration_upward)
  {
    af::init(diagonal) = 222, 222, 1, 2, 3, 222;
    af::init(superdiagonal) = 0, 0, 2, 1, 0;
    // Carefully obtain by hand using Mathematica (c.f. tst_svd.nb)
    af::init(target_diagonal) = 222, 222,
                                  0.6632629851554485, 2.800057205029313,
                                  3.2307145315672567,
                                222;
    af::init(target_superdiagonal) = 0, 0,
                                      -0.012544160425610518, 0.5311196858011309,
                                     0;
  }
};

struct golub_kahan_iteration_test_case_5 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_5()
    : golub_kahan_iteration_test_case(
                4+2, 1, 5, &svd_t::do_implicit_zero_shift_qr_iteration_downward)
  {
    af::init(diagonal) = 222, 1, 2, 3, 4, 222;
    af::init(superdiagonal) = 0, -1, -1, -1, 0;
    // Carefully obtain by hand using Mathematica (c.f. tst_svd.nb)
    af::init(target_diagonal) = 222,
                                2., -2.1213203435596433, -2.14919697074224,
                                     2.6320780861415174,
                                222;
    af::init(target_superdiagonal) = 0,
                                     -1.2247448713915892, 2.160246899469287,
                                      2.605081699820434,
                                     0;
  }
};

struct golub_kahan_iteration_test_case_6 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_6()
    : golub_kahan_iteration_test_case(
                4+2, 1, 5, &svd_t::do_implicit_zero_shift_qr_iteration_upward)
  {
    af::init(diagonal) = 222, 4, 3, 2, 1, 222;
    af::init(superdiagonal) = 0, -1, -1, -1, 0;
    // Revert the target elements from test_case_5:
    af::init(target_diagonal) = 222,
                                 2.6320780861415174, -2.1491969707422398,
                                -2.1213203435596433,  2.0,
                                222;
    af::init(target_superdiagonal) = 0,
                                      2.605081699820434 ,2.160246899469287,
                                     -1.2247448713915892,
                                     0;
  }
};

struct golub_kahan_iteration_test_case_7 : golub_kahan_iteration_test_case
{
  golub_kahan_iteration_test_case_7()
    : golub_kahan_iteration_test_case(
                9, 0, 9, &svd_t::do_implicit_zero_shift_qr_iteration_upward)
  {
    for (int i=r; i<s; ++i) {
      diagonal[i] = i+1;
      if (i < s-1) superdiagonal[i] = 2*diagonal[i];
    }
  }
};

template <class grading_func_t>
void exercise_golub_kahan_iterations(grading_func_t grading_func,
                                     double superdiagonal_multiplier,
                                     double ratio_threshold=10) {
  double thresh = ratio_threshold;

  for (int sorting=0; sorting < 2; ++sorting)
  for (int grading=0; grading < 3; ++grading)
  for (int zero_on_diag=0; zero_on_diag<2; ++zero_on_diag)
  for (int n = 3; n < 10; ++n)
  {
    af::shared<double> diagonal(n);
    af::shared<double> superdiagonal(n-1);
    for (int i=0; i<n; ++i) {
      diagonal[i] =   grading == 0 ? grading_func(i+1)
                    : grading == 1 ? grading_func(n-i)
                    : grading == 2 ? grading_func(n/2 - i)
                    :                0;
      if (i < n-1) superdiagonal[i] = superdiagonal_multiplier*diagonal[i];
    }
    if (zero_on_diag) diagonal[n/2] = 0;
    matrix_t a = upper_bidiagonal(diagonal.ref(), superdiagonal.ref());

    matrix_t identity_n = identity<double>(n);

    matrix_t u = identity_n.deep_copy();
    matrix_ref_t u_ = u.ref();

    matrix_t v = identity_n.deep_copy();
    matrix_ref_t v_ = v.ref();

    svd::bidiagonal_decomposition<double> svd(
      diagonal.ref(), superdiagonal.ref(), svd::upper_bidiagonal_kind,
      u.ref(), true, v.ref(), true);
    svd.compute();
    if (sorting) svd.sort();

    SCITBX_ASSERT(superdiagonal.all_eq(0));
    SCITBX_ASSERT(normality_ratio(u.const_ref(), svd.tol) < thresh);
    SCITBX_ASSERT(normality_ratio(v.const_ref(), svd.tol) < thresh);

    SCITBX_ASSERT(diagonal.all_ge(0));

    matrix_t a1 = svd::reconstruct(u_, v_, diagonal.const_ref());
    SCITBX_ASSERT(equality_ratio(a.const_ref(), a1.const_ref(), svd.tol)
                    < thresh);
  }
}

void exercise_singular_values_accuracy(double x) {
  int n = 20;

  // Some of the tests suggested in the ref [4]
  // (c.f. comments for svd::bidiagonal_decomposition)
  vec_t d0(n), f0(n-1);
  for (int i=0; i<n; ++i) {
    d0[i] = std::pow(x, i);
    if (i < n-1) f0[i] = d0[i];
  }
  matrix_ref_t u_, v_;

  {
    // graded from small at the upper left corner to large at the lower right
    vec_t d = d0.deep_copy(), f = f0.deep_copy(), sigma=d0;
    svd::bidiagonal_decomposition<double> svd(d.ref(), f.ref(),
                                              svd::upper_bidiagonal_kind,
                                              u_, false, v_, false);
    svd.compute();
    SCITBX_ASSERT(f.all_eq(0));
    svd.sort();
    SCITBX_ASSERT(svd.numerical_rank((1.+ x)/2) == n-1);
    SCITBX_ASSERT(svd.numerical_rank((x + x*x)/2) == n-2);
    std::reverse(d.begin(), d.end());
    vec_t delta = af::abs(d - sigma)/sigma;
    SCITBX_ASSERT( delta.all_lt(svd.tol) );
  }
  {
    // graded from large at the upper left corner to small at the lower right
    vec_t d = d0.deep_copy(), f = f0.deep_copy(), sigma=d0;
    std::reverse(d.begin(), d.end());
    std::reverse(f.begin(), f.end());
    svd::bidiagonal_decomposition<double> svd(d.ref(), f.ref(),
                                              svd::upper_bidiagonal_kind,
                                              u_, false, v_, false);
    svd.compute();
    SCITBX_ASSERT(f.all_eq(0));
    svd.sort();
    std::reverse(d.begin(), d.end());
    vec_t delta = af::abs(d - sigma)/sigma;
    SCITBX_ASSERT( delta.all_lt(svd.tol) );
  }
  {
    // graded from large at upper left corner to small at center
    // to large at lower right corner
    int n = 10;
    vec_t d0(2*n), f0(2*n-1);
    for (int i=0; i<2*n; ++i) {
      d0[i] = std::pow(x, std::abs(n-1-i));
      if (i < 2*n-1) f0[i] = d0[i];
    }
    vec_t sigma(2*n);
    af::init(sigma) = 1.e100,
                      1.41421356237309505e90, 1.e90,
                      1.22474487139158905e80, 1.e80,
                      1.15470053837925153e70, 1.e70,
                      1.11803398874989485e60, 1.e60,
                      1.09544511501033223e50, 1.e50,
                      1.08012344973464337e40, 1.e40,
                      1.06904496764969754e30, 1.e30,
                      1.06066017177982129e20, 1.e20,
                      1.05409255338945978e10, 1.e10,
                      0.316227766016837933; // Thanks Mathematica!
    vec_t d=d0.deep_copy(), f=f0.deep_copy();
    svd::bidiagonal_decomposition<double> svd(d.ref(), f.ref(),
                                              svd::upper_bidiagonal_kind,
                                              u_, false, v_, false);
    svd.compute();
    svd.sort();
    SCITBX_ASSERT(svd.numerical_rank(1.) == 2*n-1);
    SCITBX_ASSERT(svd.numerical_rank(1e95) == 1);
    SCITBX_ASSERT(svd.numerical_rank(1e85) == 3);
    SCITBX_ASSERT(f.all_eq(0));
    vec_t delta = af::abs(d - sigma)/sigma;
    SCITBX_ASSERT( delta.all_lt(svd.tol) );
  }
}

double linear_law(int i) { return i; }

struct power_law
{
  double m;
  power_law(double mult) : m(mult) {}
  // The absolute value allows for the grading to change direction
  // in the middle of the diagonal in case 2 above
  double operator()(int i) { return std::pow(m, std::abs(i)); }
};




int main() {
  exercise_2x2_decomposition();
  exercise_lower_bidiagonal_rectification();
  golub_kahan_iteration_test_case_1().exercise();
  golub_kahan_iteration_test_case_2().exercise();
  golub_kahan_iteration_test_case_3().exercise();
  golub_kahan_iteration_test_case_4().exercise();
  golub_kahan_iteration_test_case_5().exercise();
  golub_kahan_iteration_test_case_7().exercise();
  exercise_golub_kahan_iterations(linear_law, 2);
  exercise_golub_kahan_iterations(power_law(10), 2);
  exercise_golub_kahan_iterations(power_law(0.1), 2);
  exercise_golub_kahan_iterations(power_law(1.e10), 2);
  exercise_golub_kahan_iterations(power_law(1.e-10), 2);
  exercise_singular_values_accuracy(1e10);
  std::cout << "OK\n";
  return 0;
}
