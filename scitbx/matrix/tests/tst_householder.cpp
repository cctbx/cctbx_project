#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/initialiser.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/matrix/special_matrices.h>
#include <scitbx/matrix/move.h>
#include <scitbx/error.h>
#include <scitbx/random.h>
#include <scitbx/array_family/simple_io.h>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#include <scitbx/matrix/householder.h>
#include <scitbx/matrix/tests.h>
#include <scitbx/matrix/tests/utils.h>

namespace af = scitbx::af;
using scitbx::fn::approx_equal;
using namespace scitbx::matrix;

#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ == 1 \
 && defined(__i386__) && defined(__linux)
// avoid internal compiler error
#undef SCITBX_ASSERT
#define SCITBX_ASSERT(cond)
#endif

double tol = 1e-12;

struct test_case
{
  matrix_t a0;
  double thresh;

  test_case(int m, int n, double ratio_threshold=10)
    : a0(dim(m,n)), thresh(ratio_threshold)
  {}

  void check_qr(bool thin_q) {
    matrix_const_ref_t a0_ = a0.const_ref();

    matrix_t a = a0.deep_copy();
    matrix_ref_t a_ = a.ref();

    int m = a_.n_rows(), n=a_.n_columns();

    householder::qr_decomposition<double> qr(a_);

    // Get R
    matrix_t r(dim(thin_q ? std::min(m,n) : m, n));
    matrix_ref_t r_ = r.ref();
    for (int i=0; i<std::min(m,n); ++i)
    for (int j=i; j<n; ++j) {
      r_(i,j) = a_(i,j);
    }

    // Accumulate Q out-of-place
    matrix_t q = qr.q(thin_q);
    matrix_const_ref_t q_ = q.const_ref();
    SCITBX_ASSERT( normality_ratio(q_) < thresh );

    // Accumulate Q in-place if it makes sense
    if (thin_q && m >= n) {
      qr.accumulate_q_in_place();
      SCITBX_ASSERT( equality_ratio(q_, a_) < thresh );
    }

    // Check A = QR
    matrix_t q_r = af::matrix_multiply(q_, r_);
    matrix_const_ref_t q_r_ = q_r.const_ref();
    SCITBX_ASSERT( equality_ratio(a0_, q_r_) < thresh );
  }

  void check_lq(bool thin_q) {
    matrix_const_ref_t a0_ = a0.const_ref();

    matrix_t a = a0.deep_copy();
    matrix_ref_t a_ = a.ref();

    int m = a_.n_rows(), n=a_.n_columns();

    householder::lq_decomposition<double> lq(a_);

    // Get L
    matrix_t l(dim(m, thin_q ? std::min(m,n) : n));
    matrix_ref_t l_ = l.ref();
    for (int j=0; j<std::min(m,n); ++j) {
    for (int i=j; i<m; ++i)
      l_(i,j) = a_(i,j);
    }

    // Accumulate Q out-of-place
    matrix_t q = lq.q(thin_q);
    matrix_const_ref_t q_ = q.const_ref();
    SCITBX_ASSERT( normality_ratio(q_) < thresh );

    // Accumulate Q in-place if it makes sense
    if (thin_q && m <= n) {
      lq.accumulate_q_in_place();
      SCITBX_ASSERT( equality_ratio(q_, a_) < thresh );
    }

    // Check A = LQ
    matrix_t l_q = af::matrix_multiply(l_, q_);
    matrix_const_ref_t l_q_ = l_q.const_ref();
    SCITBX_ASSERT( equality_ratio(a0_, l_q_) < thresh );
  }

  void check_bidiagonalisation(bool thin) {
    matrix_const_ref_t a0_ = a0.const_ref();

    matrix_t a = a0.deep_copy();
    matrix_ref_t a_ = a.ref();

    int m = a_.n_rows(), n=a_.n_columns();

    householder::bidiagonalisation<double> bidiag(a_);

    // Get B
    matrix_t b(dim(thin ? std::min(m,n) : m, thin ? std::min(m,n) : n));
    matrix_ref_t b_ = b.ref();
    if (m >=n) copy_upper_bidiagonal(b_, a_);
    else copy_lower_bidiagonal(b_, a_);

    matrix_t u = bidiag.u(thin);
    matrix_ref_t u_ = u.ref();

    matrix_t v = bidiag.v(thin);
    matrix_ref_t v_ = v.ref();

    // Check that U and V are orthogonal
    SCITBX_ASSERT( normality_ratio(u_) < thresh );
    SCITBX_ASSERT( normality_ratio(v_) < thresh );

    // Check that A = U B V^T
    matrix_t u_b_vt = product_U_M_VT(u_, b_, v_);
    matrix_const_ref_t u_b_vt_ = u_b_vt.const_ref();
    SCITBX_ASSERT( equality_ratio(a0_, u_b_vt_) < thresh );
  }
};

struct lotkin_test_case : test_case
{
  lotkin_test_case(int m, int n) : test_case(m,n)
  {
    matrix_ref_t a_ = this->a0.ref();
    for (int i=0; i<m; ++i) for (int j=0; j<n; ++j) {
      a_(i,j) = i > 0 ? 1./(i+j+1) : 1;
    }
  }
};

struct graded_test_case : test_case
{
  householder::random_normal_matrix_generator<
    double, scitbx::boost_random::mt19937> gen;

  graded_test_case(int n, double x)
    : test_case(n,n), gen(n,n)
  {
    vec_t d(n), f(n-1);
    for (int i=0; i<n; ++i) {
      d[i] = std::pow(x, i);
      if (i < n-1) f[i] = d[i];
    }
    matrix_t u = gen.normal_matrix();
    matrix_t v = gen.normal_matrix();
    matrix_ref_t vt_ = v.ref();
    vt_.transpose_in_place();
    matrix_t b = upper_bidiagonal(d.ref(), f.ref());
    a0 = af::matrix_multiply(u.ref(), b.ref());
    a0 = af::matrix_multiply(a0.ref(), vt_);
  }
};

void exercise_householder_zeroing_vector() {
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
    SCITBX_ASSERT(expected.ref().all_approx_equal(obtained, tol));
    af::ref<double> overwritten(&x[1], 3);
    SCITBX_ASSERT(expected.ref().all_approx_equal(overwritten, tol));
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
    SCITBX_ASSERT(expected.ref().all_approx_equal(obtained, tol));
    af::ref<double> overwritten(&x[1], 3);
    SCITBX_ASSERT(expected.ref().all_approx_equal(overwritten, tol));
  }

  // Householder zeroing vector in-place
  {
    householder::reflection<double> p(4);
    af::ref<double> v(&p.v[0], 4);
    af::init(v) = 3, 1, 5, 1;
    af::tiny<double, 3> expected;
    af::init(expected) = -1, -5, -1;
    expected /= 3.;
    p.zero_vector(4);
    af::ref<double> obtained(&p.v[0], 3);
    SCITBX_ASSERT(approx_equal(p.beta, 0.5, tol))(p.beta);
    SCITBX_ASSERT(expected.ref().all_approx_equal(obtained, tol));
  }

  // Householder zeroing matrix columns or rows
  {
    int const m=6, n=7;
    matrix_t a0(dim(m,n));
    af::init(a0) = 11, 12, 13, 14, 15, 16, 17,
                    21, 22,  1, 24, 25, 26, 27,
                    31, 32, -1, 34, 35, 36, 37,
                    41, 42,  2, 44, 45, 46, 47,
                    51, 52,  1, 54, 55, 56, 57,
                    61, 62,  3, 64, 65, 66, 67;

    matrix_t a = a0.deep_copy();
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

    matrix_t a_t = af::matrix_transpose(a0.ref());

    householder::reflection<double> q(n, m,
                                      householder::applied_on_right_tag(),
                                      false);
    q.zero_vector(af::row_right_of(a_t.ref(), 2, 1));
    q.apply_on_right_to_lower_right_block(a_t.ref(), 3, 1);
    SCITBX_ASSERT(matrix_transpose(expected.ref()).all_approx_equal(a_t, tol));
  }
}

void exercise_householder_applied_to_symmetric_matrix() {
  double tol = 1e-12;
  scitbx::random::mersenne_twister rnd;
  for (int n=2; n<10; ++n) {
    vec_t x = rnd.random_double(n);
    householder::reflection<double> p(
      n, n, householder::applied_on_left_and_right_tag());
    p.zero_vector(x.ref(), false);
    matrix_t id_n = identity<double>(n);
    matrix_ref_t pp = id_n.ref();
    p.apply_on_left_to_lower_right_block(pp, 0, 0);

    symmetric_matrix_packed_u_t a(n);
    double *a_ = a.begin();
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) *a_++ = rnd.random_double();
    af::const_ref<double> r(a.begin(), a.size());
    matrix_t a1 = packed_u_as_symmetric(r);

    matrix_t b1 = product_U_M_VT(pp, a1.const_ref(), pp);
    symmetric_matrix_packed_u_t b(symmetric_as_packed_u(b1.ref(), tol), n);

    p.apply_to_lower_right_block(a.ref(), 0);
    SCITBX_ASSERT(a.all_approx_equal(b, tol));
  }
}

void exercise_householder_accumulation() {
  boost::mt19937 urng;
  boost::uniform_int<> gen(0, 64);
  for (int m=3; m<=6; ++m) for (int n=3; n<=6; ++n) {
    // set up test case
    matrix_t a(dim(m,n), -1.);
    matrix_ref_t a_ = a.ref();
    vec_t beta(std::min(m-1, n));
    for (int j=0; j<std::min(m-1, n); ++j) {
      double v2 = 1;
      for (int i=j+1; i<m; ++i) {
        a(i,j) = gen(urng);
        v2 += a(i,j)*a(i,j);
      }
      beta[j] = 2./v2;
    }
    matrix_t at = af::matrix_transpose(a.ref());
    matrix_ref_t at_ = at.ref();

    // test accumulation of the thin Q if m >= n or the full Q otherwise

    // on the right
    {
      matrix_t q(dim(m, std::min(m,n)));
      matrix_ref_t q_ = q.ref();
      householder::reflection<double> h(m, n,
                                        householder::applied_on_left_tag(),
                                        true);
      h.accumulate_factored_form_in_columns(q_, a_, beta.ref());
      h.accumulate_in_place_factored_form_in_columns(a_, beta.ref());
      for (int i=0; i<q_.n_rows(); ++i)
      for (int j=0; j<q_.n_columns(); ++j) {
        approx_equal(q_(i,j), a_(i,j), 1e-15);
      }
    }

    // on the left
    {
      std::swap(m,n);
      int p = std::min(m,n); // <|workaround for gcc 3.3.4 compilation error
      matrix_t q(dim(p, n)); // <|
      matrix_ref_t q_ = q.ref();
      householder::reflection<double> h(m, n,
                                        householder::applied_on_right_tag(),
                                        true);
      h.accumulate_factored_form_in_rows(
        q_, at_, beta.ref(), householder::product_in_reverse_row_order);
      h.accumulate_in_place_factored_form_in_rows(at_, beta.ref());
      for (int i=0; i<q_.n_rows(); ++i)
      for (int j=0; j<q_.n_columns(); ++j) {
        approx_equal(q_(i,j), at_(i,j), 1e-15);
      }
    }
  }
}

void exercise_householder() {
  // Householder QR (Lotkin matrix: ill-conditioned)
  std::vector<lotkin_test_case> cases;

  cases.push_back( lotkin_test_case(3,2) );
  cases.push_back( lotkin_test_case(2,3) );
  cases.push_back( lotkin_test_case(3,3) );

  cases.push_back( lotkin_test_case(4,3) );
  cases.push_back( lotkin_test_case(3,4) );
  cases.push_back( lotkin_test_case(4,4) );

  cases.push_back( lotkin_test_case(5,3) );
  cases.push_back( lotkin_test_case(3,5) );

  for (int i=0; i<cases.size(); ++i) {
    lotkin_test_case &t = cases[i];
    t.check_qr(false); // full QR
    t.check_qr(true);  // thin QR
    t.check_lq(false); // full LQ
    t.check_lq(true);  // thin LQ
  }

  graded_test_case t(10, 1.e-10);
  t.check_qr(true);
}

void exercise_bidiagonalisation() {
  /* Householder bidiagonalisation (Golub and Van Loan example 5.4.2)
     The last diagonal entry in B is zero which makes it interesting.
     In the book, the authors used higher a precision than available on standard
     hardware: hence the difference between their result and those obtained
     here. But the SVD is still correct to machine precision as we assert it is.
  */
  {
    matrix_t a0(dim(4,3));
    af::init(a0) =  1,  2,  3,
                    4,  5,  6,
                    7,  8,  9,
                   10, 11, 12;
    matrix_const_ref_t a0_ = a0.const_ref();
    matrix_t a = a0.deep_copy();
    householder::bidiagonalisation<double> bidiag(a.ref());
    matrix_t b(dim(4,3));
    af::init(b) = 12.8840987267251, 21.876432827428,  0             ,
                   0              ,  2.246235240294, -0.613281332054,
                   0              ,  0             ,  0             ,
                   0              ,  0             ,  0             ;
    matrix_const_ref_t b_ = b.const_ref();
    matrix_t u = bidiag.u(false);
    matrix_ref_t u_ = u.ref();
    matrix_t v = bidiag.v(false);
    matrix_ref_t v_ = v.ref();
    SCITBX_ASSERT( normality_ratio(u_) < 10 );
    SCITBX_ASSERT( normality_ratio(v_) < 10 );
    matrix_t u_b_vt = product_U_M_VT(u_, b_, v_);
    matrix_const_ref_t u_b_vt_ = u_b_vt.const_ref();
    SCITBX_ASSERT( equality_ratio(u_b_vt_, a0_, 1e-13) < 10 );
  }
  {
    matrix_t a0(dim(3,3));
    af::init(a0) = 24.3104915623, 0.0          , 0.0                ,
                   26.4083512403, 1.26450969612, 0.0                ,
                   28.5062109182, 2.52901939223, 3.74525472711e-15 ;
    matrix_const_ref_t a0_ = a0.const_ref();
    matrix_t a = a0.deep_copy();
    matrix_ref_t a_ = a.ref();
    householder::bidiagonalisation<double> bidiag(a.ref());
    matrix_t u = bidiag.u();
    matrix_ref_t u_ = u.ref();
    matrix_t v = bidiag.v();
    matrix_ref_t v_ = v.ref();
    SCITBX_ASSERT( normality_ratio(u_) < 10 );
    SCITBX_ASSERT( normality_ratio(v_) < 10 );
    matrix_t b(dim(3, 3));
    matrix_ref_t b_ = b.ref();
    copy_upper_bidiagonal(b_, a_);
    matrix_t u_b_vt = product_U_M_VT(u_, b_, v_);
    matrix_const_ref_t u_b_vt_ = u_b_vt.const_ref();
    SCITBX_ASSERT( equality_ratio(u_b_vt_, a0_, 1e-13) < 10 );
  }

  // Householder bidiagonalisation (Lotkin matrix)
  {
    lotkin_test_case t(7,5);
    t.check_bidiagonalisation(true);
    t.check_bidiagonalisation(false);
  }
  {
    lotkin_test_case t(10,5);
    t.check_bidiagonalisation(true);
    t.check_bidiagonalisation(false);
  }
  {
    lotkin_test_case t(5,5);
    t.check_bidiagonalisation(true);
    t.check_bidiagonalisation(false);
  }
  {
    lotkin_test_case t(5,7);
    t.check_bidiagonalisation(true);
    t.check_bidiagonalisation(false);
  }
  {
    lotkin_test_case t(5,10);
    t.check_bidiagonalisation(true);
    t.check_bidiagonalisation(false);
  }
  {
    graded_test_case t(10, 1.e10);
    t.check_bidiagonalisation(false);
  }
}

int main() {
  exercise_householder_zeroing_vector();
  exercise_householder_applied_to_symmetric_matrix();
  exercise_householder_accumulation();
  exercise_householder();
  exercise_bidiagonalisation();
  std::cout << "OK\n";
  return 0;
}
