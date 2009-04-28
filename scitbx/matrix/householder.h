#ifndef SCITBX_MATRIX_HOUSEHOLDER_H
#define SCITBX_MATRIX_HOUSEHOLDER_H

#include <scitbx/error.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/accessors/row_and_column.h>
#include <scitbx/math/accumulators.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace scitbx { namespace matrix { namespace householder {

struct applied_on_left_tag {};
struct applied_on_right_tag {};
struct applied_on_left_and_right_tag {};

/// Householder reflection P = I - beta v v^T
/** Reference: Golub and Van Loan 5.1.2 to 5.1.4
    Implementation note: only row-major storage is supported and
    this code is optimised accordingly (only relevant for big matrices).
*/
template <typename FloatType>
struct reflection
{
  typedef FloatType scalar_t;

  /// Normalisation of the Householder vector
  scalar_t beta;

  /// Normalisation of the vector x passed to the constructor
  scalar_t norm_x;

  /// Essential part of the Householder vector
  std::vector<scalar_t> v;

  /// Working vector for applying the Householder reflection to a matrix
  std::vector<scalar_t> w;

  template <class AccessorType>
  reflection(af::ref<scalar_t, AccessorType> const &x)
    : v(x.size()-1)
  {
    zero_vector(x);
  }

  reflection(int m, int n, applied_on_left_tag, bool accumulate)
    : v(m), w(accumulate ? std::max(m,n) : n)
  {}

  reflection(int m, int n, applied_on_right_tag, bool accumulate)
    : v(n), w(accumulate ? std::max(m,n) : m)
  {}

  reflection(int m, int n, applied_on_left_and_right_tag)
    : v(std::max(m,n)), w(std::max(m,n))
  {}

  /// Construct the Householder reflection P s.t. Px = ||x||_2 e_1
  /** x(0) is overwritten with the first element of Px
      whereas x(1:) is overwritten with the essential part of the Householder
      vector. Since operator() is used to access the elements of x, that code
      works even if x is stored with a stride different of 1
      Reference: Algorithm 5.1.1
  */
  template <class AccessorType>
  void zero_vector(af::ref<scalar_t, AccessorType> const &x) {
    using namespace math::accumulator;
    // compute beta and v
    int n = x.size();
    norm_accumulator<scalar_t> norm_accu;
    for (int i=1; i<n; ++i) {
      norm_accu(x(i));
      v[i-1] = x(i);
    }
    /* If I was a good boy, I would use norm_accu.norm here to get sqrt(sigma)
       and then compute sqrt(x(0)^2 + sqrt(sigma)^2) in a safe manner */
    scalar_t sigma = norm_accu.sum_sq();
    if (sigma == 0) {
      beta = 0;
      return;
    }
    scalar_t mu = norm_x = std::sqrt(x(0)*x(0) + sigma);
    scalar_t v0 = x(0) <= 0 ? x(0) - mu : -sigma/(x(0) + mu);
    beta = 2*v0*v0/(sigma + v0*v0);

    // overwrite
    x(0) = norm_x;
    for (int i=1; i<n; ++i) {
      v[i-1] /= v0;
      x(i) = v[i-1];
    }
  }

  /// Replace A(i:, j:) by PA(i:, j:)
  void apply_on_left_to_lower_right_block(
    af::ref<scalar_t, af::mat_grid> const &a, int i, int j)
  {
    int m = a.n_rows(), n=a.n_columns();
    // w = beta A(i:, j:)^T v
    for (int jj=j; jj < n; ++jj) w[jj-j] = a(i,jj);
    for (int ii=i+1; ii < m; ++ii)
      for (int jj=j; jj < n; ++jj) {
        w[jj-j] += a(ii,jj)*v[ii-i-1];
    }
    for (int k=0; k < n-j; ++k) w[k] *= beta;
    // A(i:, j:) = PA(i:, j:)
    for (int jj=j; jj < n; ++jj) a(i,jj) -= w[jj-j];
    for (int ii=i+1; ii < m; ++ii) {
      for (int jj=j; jj < n; ++jj) a(ii,jj) -= v[ii-i-1]*w[jj-j];
    }
  }

  /// Replace A(i:, j:) by A(i:, j:)P
  void apply_on_right_to_lower_right_block(
    af::ref<scalar_t, af::mat_grid> const &a, int i, int j)
  {
    int m = a.n_rows(), n=a.n_columns();
    // w = beta A(i:, j) v
    for (int ii=i; ii < m; ++ii) {
      w[ii-i] = a(ii, j);
      for (int jj=j+1; jj < n; ++jj) w[ii-i] += a(ii,jj)*v[jj-j-1];
      w[ii-i] *= beta;
    }
    // A(i:, j:) = A(i:, j:)P
    for (int ii=i; ii < m; ++ii) {
      a(ii,j) -= w[ii-i];
      for (int jj=j+1; jj < n; ++jj) a(ii,jj) -= w[ii-i]*v[jj-j-1];
    }
  }

  /// Accumulate the Householder reflection stored in factored form in
  /// the columns of the matrix \c a given the corresponding beta's
  /** C.f. Golub and Van Loan section 5.1.6 */
  void
  accumulate_factored_form_in_columns(af::ref<scalar_t, af::mat_grid> const &q,
                                      af::const_ref<scalar_t, af::mat_grid> const &a,
                                      af::const_ref<scalar_t> const &beta)
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(q.is_square());
    SCITBX_ASSERT(q.n_rows() == m)(q.n_rows())(m); // A=QR
    SCITBX_ASSERT(beta.size() == n);
    q.set_identity();
    for (int j=n-1; j >= 0; --j) {
      for (int i=j+1; i < m; ++i) v[i-j-1] = a(i,j);
      this->beta = beta[j];
      apply_on_left_to_lower_right_block(q, j, j);
    }
  }

  /// Accumulate the Householder reflection stored in factored form in
  /// the rows of the matrix \c a given the corresponding beta's
  /** off_diag is how far to the right of the diagonal was the first element
      of the row zeroed by the reflection now stored in that row.
      (off_diag = 0 corresponds to the same behaviour as
      \c accumulate_factored_form_in_columns but for rows)
      C.f. Golub and Van Loan section 5.1.6
  */
  void
  accumulate_factored_form_in_rows(af::ref<scalar_t, af::mat_grid> const &q,
                                   af::const_ref<scalar_t, af::mat_grid> const &a,
                                   af::const_ref<scalar_t> const &beta,
                                   int const off_diag=0)
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(q.is_square());
    SCITBX_ASSERT(q.n_columns() == n)(q.n_columns())(n); // A=RQ
    SCITBX_ASSERT(beta.size() == n - off_diag - 1);
    q.set_identity();
    for (int i=n - off_diag - 2; i >= 0; --i) {
      for (int j=i + off_diag + 1; j < n; ++j) v[j - i - off_diag - 1] = a(i,j);
      this->beta = beta[i];
      apply_on_left_to_lower_right_block(q, i+off_diag, i+off_diag);
    }
  }
};


/// Householder QR decomposition
/** Reference: Golub and Van Loan, Algorithm 5.2.1 */
template <typename FloatType>
struct qr_decomposition
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;

  reflection<scalar_t> p;
  std::vector<scalar_t> beta;
  matrix_t q;

  qr_decomposition(af::ref<scalar_t, af::mat_grid> const &a,
                   bool accumulate_q=false)
    : p(a.n_rows(), a.n_columns(), applied_on_left_tag(), accumulate_q)
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(m >= n)(m)(n);
    beta.reserve(n);
    for (int j=0; j<n; ++j) {
      p.zero_vector(af::column_below(a, j, j));
      beta.push_back(p.beta);
      p.apply_on_left_to_lower_right_block(a, j, j+1);
    }

    if (accumulate_q) {
      af::const_ref<scalar_t> beta_(&beta[0], n);
      q = matrix_t(dim(m,m), af::init_functor_null<scalar_t>());
      p.accumulate_factored_form_in_columns(q.ref(), a, beta_);
    }
  }

};


/// Decomposition U^T A V = B where B is upper bidiagonal and U,V orthogonal
/** Reference: Golub and Van Loan, Algorithm 5.4.2 */
template <typename FloatType>
struct bidiagonalisation
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;

  reflection<FloatType> p;
  std::vector<scalar_t> beta_left, beta_right;
  matrix_t u, v;

  bidiagonalisation(af::ref<scalar_t, af::mat_grid> const &a,
                    bool accumulate_u=false,
                    bool accumulate_v=false)
    : p(a.n_rows(), a.n_columns(), applied_on_left_and_right_tag())
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(m >= n)(m)(n);
    int r_left = n, r_right = n-2; // number of Householder reflections
    beta_left.reserve(r_left);
    beta_right.reserve(r_right);
    for (int j=0; j<n; ++j) {
      p.zero_vector(af::column_below(a, j, j));
      beta_left.push_back(p.beta);
      p.apply_on_left_to_lower_right_block(a, j, j+1);
      if (j < r_right) {
        p.zero_vector(af::row_right_of(a, j, j+1));
        beta_right.push_back(p.beta);
        p.apply_on_right_to_lower_right_block(a, j+1, j+1);
      }
    }
    if (accumulate_u) {
      af::const_ref<scalar_t> beta_left_(&beta_left[0], beta_left.size());
      u = matrix_t(dim(m,m), af::init_functor_null<scalar_t>());
      p.accumulate_factored_form_in_columns(u.ref(), a, beta_left_);
    }
    if (accumulate_v) {
      af::const_ref<scalar_t> beta_right_(&beta_right[0], beta_right.size());
      v = matrix_t(dim(n,n), af::init_functor_null<scalar_t>());
      p.accumulate_factored_form_in_rows(v.ref(), a, beta_right_, 1);
    }
  }
};


}}} // scitbx::matrix::householder

#endif // GUARD
