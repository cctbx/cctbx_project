#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/initialiser.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/matrix/special_matrices.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/error.h>
#include <scitbx/array_family/simple_io.h>
#include <iostream>

#include <scitbx/matrix/householder.h>
#include <scitbx/matrix/tests.h>

namespace af = scitbx::af;
using scitbx::fn::approx_equal;
using namespace scitbx::matrix;

typedef af::c_grid<2> dim;
typedef af::versa<double, dim> matrix_t;
typedef af::ref<double, dim> matrix_ref_t;
typedef af::const_ref<double, dim> matrix_const_ref_t;

double tol = 1e-12;

struct test_case
{
  matrix_t a__;
  double thresh;

  test_case(int m, int n, double ratio_threshold=10)
    : a__(dim(m,n)), thresh(ratio_threshold)
  {}

  void check_qr(bool thin_q) {
    matrix_t a = a__.deep_copy();
    matrix_ref_t a_ = a.ref();

    int m = a.accessor()[0], n=a.accessor()[1];

    householder::qr_decomposition<double> qr(a_,
                                             true, // accumulate_q
                                             thin_q);
    matrix_ref_t q_ = qr.q.ref();

    // Check Q is orthogonal
    SCITBX_ASSERT(normality_ratio(qr.q.const_ref()) < thresh);

    // Get R
    matrix_t r(dim(thin_q ? std::min(m,n) : m, n));
    matrix_ref_t r_ = r.ref();
    for (int i=0; i<std::min(m,n); ++i)
    for (int j=i; j<n; ++j) {
      r_(i,j) = a_(i,j);
    }

    // Check A = QR
    matrix_t q_r = af::matrix_multiply(q_, r_);
    double delta_r = equality_ratio(a__.const_ref(), q_r.const_ref());
    SCITBX_ASSERT(delta_r < thresh);
  }

  void check_bidiagonalisation(bool debug=false) {
    matrix_t a = a__.deep_copy();
    matrix_ref_t a_ = a.ref();

    int m = a.accessor()[0], n=a.accessor()[1];

    householder::bidiagonalisation<double> bidiag(a_);
    matrix_ref_t u_ = bidiag.u.ref();
    matrix_ref_t v_ = bidiag.v.ref();

    // Check that U and V are orthogonal
    SCITBX_ASSERT(normality_ratio(bidiag.u.const_ref()) < thresh);
    SCITBX_ASSERT(normality_ratio(bidiag.v.const_ref()) < thresh);

    // Get B
    matrix_t b(dim(m,n));
    matrix_ref_t b_ = b.ref();
    for (int i=0; i<n; ++i) for (int j=i; j<i+2 && j<n; ++j) b_(i,j) = a_(i,j);

    // Check that A = U B V^T
    matrix_t u_b = af::matrix_multiply(u_, b_);
    matrix_t vt = af::matrix_transpose(v_);
    matrix_t u_b_vt = af::matrix_multiply(u_b.ref(), vt.ref());
    SCITBX_ASSERT(equality_ratio(a__.const_ref(), u_b_vt.const_ref()) < thresh);
  }
};

struct lotkin_test_case : test_case
{
  lotkin_test_case(int m, int n) : test_case(m,n)
  {
    matrix_ref_t a_ = this->a__.ref();
    for (int i=0; i<m; ++i) for (int j=0; j<n; ++j) {
      a_(i,j) = i > 0 ? 1./(i+j+1) : 1;
    }
  }
};

int main() {
  // Householder zeroing vector (1)
  {
    af::tiny<double, 4> x;
    af::init(x) = 3, 1, 5, 1;
    householder::reflection<double> p(x.ref());
    af::tiny<double, 3> expected;
    af::init(expected) = -1, -5, -1;
    expected /= 3.;
    af::ref<double> obtained(&p.v[0], 3);
    SCITBX_ASSERT(approx_equal(p.beta, 0.5, tol))(p.beta);
    SCITBX_ASSERT(expected.ref().all_approx_equal(obtained, tol))
                 (obtained);
    af::ref<double> overwritten(&x[1], 3);
    SCITBX_ASSERT(expected.ref().all_approx_equal(overwritten, tol))
                 (overwritten);
  }

  // Householder zeroing vector (2)
  {
    af::tiny<double, 4> x;
    af::init(x) = -3, 1, 5, 1;
    householder::reflection<double> p(x.ref());
    af::tiny<double, 3> expected;
    af::init(expected) = -1, -5, -1;
    expected /= 9.;
    af::ref<double> obtained(&p.v[0], 3);
    SCITBX_ASSERT(approx_equal(p.beta, 3./2, tol))(p.beta);
    SCITBX_ASSERT(expected.ref().all_approx_equal(obtained, tol))
                 (obtained);
    af::ref<double> overwritten(&x[1], 3);
    SCITBX_ASSERT(expected.ref().all_approx_equal(overwritten, tol))
                 (overwritten);
  }

  // Householder zeroing matrix columns or rows
  {
    int const m=6, n=7;
    matrix_t a__(dim(m,n));
    af::init(a__) = 11, 12, 13, 14, 15, 16, 17,
                    21, 22,  1, 24, 25, 26, 27,
                    31, 32, -1, 34, 35, 36, 37,
                    41, 42,  2, 44, 45, 46, 47,
                    51, 52,  1, 54, 55, 56, 57,
                    61, 62,  3, 64, 65, 66, 67;

    matrix_t a = a__.deep_copy();
    householder::reflection<double> p(m, n,
                                      householder::applied_on_left_tag(),
                                      false);
    p.zero_vector(af::column_below(a.ref(), 1, 2));
    p.apply_on_left_to_lower_right_block(a.ref(), 1, 3);
    /* Mathematica:
      a = Table[10 i + j, {i, 6}, {j, 7}];
      a[[2 ;;, 3]] = {1, -1, 2, 1, 3};
      a // MatrixForm
      v = {0, 1, 1/3, -2/3, -1/3, -1};
      beta = 3/4;
      p = IdentityMatrix[6] - beta KroneckerProduct[v, v];
      b = p.a;
      b // MatrixForm
    */
    matrix_t expected(dim(m,n));
    af::init(expected) =  11, 12, 13   ,  14,  15   ,  16   ,  17   ,
                          21, 22,  4   ,  81, 165./2,  84   , 171./2,
                          31, 32,  1./3,  53, 325./6, 166./3, 113./2,
                          41, 42, -2./3,   6,  20./3,  22./3,   8   ,
                          51, 52, -1./3,  35, 215./6, 110./3,  75./2,
                          61, 62, -1.  ,   7,  15./2,   8   ,  17./2;
    SCITBX_ASSERT(expected.all_approx_equal(a, tol));

    matrix_t a_t = af::matrix_transpose(a__.ref());

    householder::reflection<double> q(n, m,
                                      householder::applied_on_right_tag(),
                                      false);
    q.zero_vector(af::row_right_of(a_t.ref(), 2, 1));
    q.apply_on_right_to_lower_right_block(a_t.ref(), 3, 1);
    SCITBX_ASSERT(matrix_transpose(expected.ref()).all_approx_equal(a_t, tol));
  }

  // Householder QR (Lotkin matrix: ill-conditioned)
  {
    lotkin_test_case t(3,2);
    t.check_qr(false); // full QR
    t.check_qr(true);  // thin QR
  }
  {
    lotkin_test_case t(3,3);
    t.check_qr(false); // full QR
    t.check_qr(true);  // thin QR
  }
  {
    lotkin_test_case t(4,3);
    t.check_qr(false); // full QR
    t.check_qr(true);  // thin QR
  }
  {
    lotkin_test_case t(5,3);
    t.check_qr(false); // full QR
    t.check_qr(true);  // thin QR
  }

  /* Householder bidiagonalisation (Golub and Van Loan example 5.4.2)
     The last diagonal entry in B is zero which makes it interesting.
     In the book, the authors used higher a precision than available on standard
     hardware: hence the difference between their result and those obtained
     here. But the SVD is still correct to machine precision as we assert it is.
  */
  {
    matrix_t a__(dim(4,3));
    af::init(a__) =  1,  2,  3,
                     4,  5,  6,
                     7,  8,  9,
                    10, 11, 12;
    matrix_t a = a__.deep_copy();
    householder::bidiagonalisation<double> bidiag(a.ref(), true, true);
    matrix_t b(dim(4,3));
    af::init(b) = 12.8840987267251, 21.876432827428,  0             ,
                   0              ,  2.246235240294, -0.613281332054,
                   0              ,  0             ,  0             ,
                   0              ,  0             ,  0             ;
    matrix_ref_t u_ = bidiag.u.ref();
    matrix_ref_t v_ = bidiag.v.ref();
    SCITBX_ASSERT(af::matrix_multiply(u_, af::matrix_transpose(u_).ref())
                  .all_approx_equal(identity<double>(4), tol));
    SCITBX_ASSERT(af::matrix_multiply(v_, af::matrix_transpose(v_).ref())
                  .all_approx_equal(identity<double>(3), tol));
    matrix_t ut = af::matrix_transpose(u_);
    matrix_t ut_a = af::matrix_multiply(ut.ref(), a__.ref());
    matrix_t ut_a_v = af::matrix_multiply(ut_a.ref(), v_);
    SCITBX_ASSERT(b.all_approx_equal(ut_a_v, tol));
  }

  // Householder bidiagonalisation (Lotkin matrix)
  {
    lotkin_test_case t(7,5);
    t.check_bidiagonalisation();
  }
  {
    lotkin_test_case t(10,5);
    t.check_bidiagonalisation();
  }
  {
    lotkin_test_case t(5,5);
    t.check_bidiagonalisation();
  }
  std::cout << "OK\n";
  return 0;
}
