#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/shared_reductions.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/matrix/tests.h>
#include <scitbx/matrix/tests/utils.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/matrix/householder.h>
#include <scitbx/random.h>
#include <iostream>

namespace af = scitbx::af;
using namespace scitbx::matrix;

template <class TestType>
struct test_case
{
  symmetric_matrix_packed_l_t a0_l;
  symmetric_matrix_packed_u_t a0_u;
  double thresh;

  TestType test;

  test_case(int n, double ratio_threshold=10.)
    : a0_l(n), thresh(ratio_threshold)
  {
    a0_u = test.create_packed_u(n);
    symmetric_matrix_packed_u_ref_t a0_u_ = a0_u.ref();
    symmetric_matrix_packed_l_ref_t a0_l_ = a0_l.ref();
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
      a0_l_(j,i) = a0_u_(i,j);
    }
  }

  void check_L_L_transpose() {
    symmetric_matrix_packed_l_t a = a0_l.deep_copy();
    symmetric_matrix_packed_l_ref_t a_ = a.ref(), &l_ = a_;
    int n = a_.accessor().n;
    cholesky::l_l_transpose_decomposition_in_place<double> chol(a_);
    cholesky::failure_info<double> fail = chol.failure;
    if (fail) {
      test.check_failure(n, fail);
      return;
    }
    double *a__ = a.begin(), *a0__ = a0_l.begin();
    matrix_t l(dim(n,n)), l_trans(dim(n,n)), a1(dim(n,n));
    for (int i=0; i<n; ++i) for (int j=0; j<=i; ++j) {
      l(i,j) = *a__;
      l_trans(j,i) = *a__++;
      a1(i,j) = a1(j,i) = *a0__++;
    }

    matrix_t l_l_trans = af::matrix_multiply(l.ref(), l_trans.ref());
    double e = equality_ratio(l_l_trans.ref(), a1.ref());
    SCITBX_ASSERT(e < thresh);

    scitbx::random::mersenne_twister gen;
    for (int k=0; k<5; ++k) {
      vec_t x = gen.random_double(n);
      vec_t b = x.deep_copy();
      cholesky::solve_in_place::using_l_l_transpose(l_, x.ref());
      double e = cholesky_test_ratio(a1.ref(), x.ref(), b.ref());
      SCITBX_ASSERT( e < thresh )(e);
    }
  }

  void check_U_transpose_U() {
    symmetric_matrix_packed_u_t a = a0_u.deep_copy();
    symmetric_matrix_packed_u_ref_t a_ = a.ref(), &u_ = a_;
    int n = a_.accessor().n;
    cholesky::u_transpose_u_decomposition_in_place<double> chol(a_);
    cholesky::failure_info<double> fail = chol.failure;
    if (fail) {
      test.check_failure(n, fail);
      return;
    }
    double *a__ = a.begin(), *a0__ = a0_u.begin();
    matrix_t u(dim(n,n)), u_trans(dim(n,n)), a1(dim(n,n));
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
      u(i,j) = *a__;
      u_trans(j,i) = *a__++;
      a1(i,j) = a1(j,i) = *a0__++;
    }

    matrix_t u_trans_u = af::matrix_multiply(u_trans.ref(), u.ref());
    double e = equality_ratio(u_trans_u.ref(), a1.ref());
    SCITBX_ASSERT(e < thresh);

    symmetric_matrix_packed_u_t inv_a = cholesky::inverse_of_u_transpose_u(u_);
    double f = residual_of_symmetric(a0_u.ref(), inv_a.const_ref());
    SCITBX_ASSERT( f < thresh )(f);

    scitbx::random::mersenne_twister gen;
    for (int k=0; k<5; ++k) {
      vec_t x = gen.random_double(n);
      vec_t b = x.deep_copy();
      cholesky::solve_in_place::using_u_transpose_u(u_, x.ref());
      double e = cholesky_test_ratio(a1.ref(), x.ref(), b.ref());
      SCITBX_ASSERT( e < thresh )(e);
      SCITBX_ASSERT( e < thresh );
    }
  }

  void exercise() {
    check_L_L_transpose();
    check_U_transpose_U();
  }
};


struct hilbert
{
  symmetric_matrix_packed_u_t create_packed_u(int n)
  {
    symmetric_matrix_packed_u_t result(n);
    double *a = result.begin();
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
      *a++ = 1./(i+j+1);
    }
    return result;
  }

  void check_failure(int n, cholesky::failure_info<double> const &fail) {
    SCITBX_ASSERT(n >= 14);
    SCITBX_ASSERT(fail.index == 13 || fail.index == 14);
  }
};


struct random_test
{
  typedef householder::random_normal_matrix_generator<
            double,
            scitbx::boost_random::mt19937>
          random_householder_gen_t;

  scitbx::boost_random::mt19937 urng;
  vec_t lambda;

  symmetric_matrix_packed_u_t create_packed_u(int n)
  {
    random_householder_gen_t gen(urng, n, n);
    fill_eigenvalues();
    return gen.symmetric_matrix_with_eigenvalues(lambda.ref());
  }

  virtual void fill_eigenvalues()=0;

  virtual void check_failure(int n, cholesky::failure_info<double> const &fail)
  {}
};

struct condition_1 : random_test
{
  virtual void fill_eigenvalues() {
    vec_ref_t l = lambda.ref();
    for (int i=0; i<l.size(); ++i) {
      l[i] = std::pow(10., i);
    }
  }
};

struct condition_2 : random_test
{
  virtual void fill_eigenvalues() {
    vec_ref_t l = lambda.ref();
    for (int i=0; i<l.size(); ++i) {
      l[i] = std::pow(10., -i);
    }
  }
};

struct condition_3 : random_test
{
  virtual void fill_eigenvalues() {
    vec_ref_t l = lambda.ref();
    for (int i=0; i<l.size()/2; ++i) {
      l[i] = std::pow(10., -i);
    }
    for (int i=l.size()/2; i<l.size(); ++i) {
      l[i] = std::pow(10., i);
    }
  }
};

struct condition_4 : random_test
{
  virtual void fill_eigenvalues() {
    vec_ref_t l = lambda.ref();
    for (int i=0; i<l.size(); ++i) {
      l[i] = std::pow(10., i);
    }
    scitbx::random::mersenne_twister gen;
    af::shared<std::size_t> perm = gen.random_permutation(l.size());
    lambda = af::select(l, perm.const_ref());
  }
};

int main() {
  {
    for (int n=1; n<20; ++n) {
      test_case<hilbert> t(n);
      t.exercise();
    }
  }
  {
    for (int n=2; n<20; ++n) {
      test_case<condition_1> t(n);
      t.exercise();
    }
  }
  {
    for (int n=2; n<20; ++n) {
      test_case<condition_2> t(n);
      t.exercise();
    }
  }
  {
    for (int n=2; n<20; ++n) {
      test_case<condition_3> t(n);
      t.exercise();
    }
  }
  {
    for (int n=2; n<20; ++n) {
      test_case<condition_4> t(n);
      t.exercise();
    }
  }
  std::cout << "OK\n";
}
