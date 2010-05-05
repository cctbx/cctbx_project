#ifndef SCITBX_MATRIX_HOUSEHOLDER_H
#define SCITBX_MATRIX_HOUSEHOLDER_H

#include <scitbx/error.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/accessors/row_and_column.h>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/accumulators.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace scitbx { namespace matrix { namespace householder {

struct applied_on_left_tag {};
struct applied_on_right_tag {};
struct applied_on_left_and_right_tag {};

enum { product_in_row_order, product_in_reverse_row_order };

/// Householder reflection P = I - beta v v^T
/** Reference: Golub and Van Loan 5.1.2 to 5.1.4
    Implementation note: only row-major storage is supported and
    this code is optimised accordingly (only relevant for big matrices).
*/
template <typename FloatType>
struct reflection
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;
  typedef af::ref<scalar_t, af::mat_grid> matrix_ref_t;
  typedef af::const_ref<scalar_t, af::mat_grid> matrix_const_ref_t;
  typedef af::versa<double, af::packed_u_accessor> symmetric_matrix_packed_u_t;
  typedef af::ref<double, af::packed_u_accessor> symmetric_matrix_packed_u_ref_t;

  /// Normalisation of the Householder vector
  scalar_t beta;

  /// Normalisation of the vector x passed to the constructor
  scalar_t norm_x;

  /// Essential part of the Householder vector
  std::vector<scalar_t> v;

  /// Working vector for applying the Householder reflection to a matrix
  std::vector<scalar_t> w;

  /** @name Constructors
      One may think that the working array v should be of size n-1 for a problem
      of size n since the real vector v is actually [ 1 v[0] v[1] ... ]
      but we actually use transformations in place in several places where
      v is loaded with [ x[0] x[1] ... ] and then transformed to
      [ 1 v[0] v[1] ... ]. Hence the need for one more element.
   */
  //@{
  template <class AccessorType>
  reflection(af::ref<scalar_t, AccessorType> const &x)
    : v(x.size())
  {
    zero_vector(x);
  }

  reflection(int n)
    : v(n)
  {}

  reflection(int m, int n, applied_on_left_tag, bool accumulate)
    : v(m), w(accumulate ? std::max(m,n) : n)
  {}

  reflection(int m, int n, applied_on_right_tag, bool accumulate)
    : v(n), w(accumulate ? std::max(m,n) : m)
  {}

  reflection(int m, int n, applied_on_left_and_right_tag)
    : v(std::max(m,n)), w(std::max(m,n))
  {}

  //@}

  /// Construct the Householder reflection P s.t. Px = ||x||_2 e_1
  /** If requested, x(0) is overwritten with the first element of Px
      whereas x(1:) is overwritten with the essential part of the Householder
      vector. Since operator() is used to access the elements of x, that code
      works even if x is stored with a stride different of 1

      v is filled with the essential part of the Householder vector.

      Reference: Algorithm 5.1.1 (with the substitution sigma -> sqrt(sigma))
  */
  template <class AccessorType>
  void zero_vector(af::ref<scalar_t, AccessorType> const &x, bool overwrite=true)
  {
    using namespace math::accumulator;
    // compute beta and v(0)
    int n = x.size();
    norm_accumulator<scalar_t> norm_accu;
    for (int i=1; i<n; ++i) norm_accu(x(i));
    scalar_t sigma = norm_accu.norm();
    if (sigma == 0) {
      beta = 0;
      return;
    }
    scalar_t mu = norm_x = math::norm(x(0), sigma);
    scalar_t v0 = x(0) <= 0 ? x(0) - mu : -sigma/(x(0) + mu) * sigma;
    scalar_t sigma_1 = sigma/v0;
    beta = 2/(sigma_1*sigma_1 + 1);

    // compute v(1:) and overwrite if requested
    if (overwrite) {
      x(0) = norm_x;
      for (int i=1; i<n; ++i) {
        x(i) /= v0;
        v[i-1] = x(i);
      }
    }
    else {
      for (int i=1; i<n; ++i) v[i-1] = x(i)/v0;
    }
  }

  /// Overload taking the vector x from v(0:n), i.e. operating in-place on v
  void zero_vector(int n) {
    zero_vector(af::ref<scalar_t>(&v[0], n), false);
  }

  /// Replace A(i:, j:) by PA(i:, j:)
  void apply_on_left_to_lower_right_block(matrix_ref_t const &a, int i, int j)
  {
    int m = a.n_rows(), n=a.n_columns();
    // w = beta A(i:, j:)^T v
    for (int jj=j; jj < n; ++jj) w[jj-j] = a(i,jj);
    for (int ii=i+1; ii < m; ++ii)
      for (int jj=j; jj < n; ++jj) {
        w[jj-j] += a(ii,jj)*v[ii-i-1];
    }
    for (int k=0; k < n-j; ++k) w[k] *= beta;
    // A(i:, j:) -= v(i:) w(j:)^T
    for (int jj=j; jj < n; ++jj) a(i,jj) -= w[jj-j];
    for (int ii=i+1; ii < m; ++ii) {
      for (int jj=j; jj < n; ++jj) a(ii,jj) -= v[ii-i-1]*w[jj-j];
    }
  }

  /// Replace A(i:, j:) by A(i:, j:)P
  void apply_on_right_to_lower_right_block(matrix_ref_t const &a, int i, int j)
  {
    int m = a.n_rows(), n=a.n_columns();
    // w = beta A(i:, j) v
    for (int ii=i; ii < m; ++ii) {
      w[ii-i] = a(ii, j);
      for (int jj=j+1; jj < n; ++jj) w[ii-i] += a(ii,jj)*v[jj-j-1];
      w[ii-i] *= beta;
    }
    // A(i:, j:) -= w(i:) v(j:)^T
    for (int ii=i; ii < m; ++ii) {
      a(ii,j) -= w[ii-i];
      for (int jj=j+1; jj < n; ++jj) a(ii,jj) -= w[ii-i]*v[jj-j-1];
    }
  }

  /// Replace A(i:, i:) by PA(i:,i:)P for a symmetric matrix A
  void
  apply_to_lower_right_block(symmetric_matrix_packed_u_ref_t const &a_, int i0)
  {
    int n = a_.accessor().n;
    scalar_t *a0 = &a_(i0,i0);

    // w = beta A(i:,i:) v
    scalar_t *a = a0;
    w[0] = *a++;
    for (int k=i0+1; k<n; ++k) {
      scalar_t a_i0_k = *a++, &a_k_i0 = a_i0_k;
      w[0] += a_i0_k * v[k-i0-1];
      w[k-i0] = a_k_i0;
    }
    symmetric_packed_u_vector(n-i0-1, a, &v[0], &w[1], 1., 1.);
    for (int k=0; k < n-i0; ++k) w[k] *= beta;

    // nu = v^T w
    scalar_t nu = w[0];
    for (int i=i0+1; i<n; ++i) nu += v[i-i0-1]*w[i-i0];
    scalar_t beta_nu = beta * nu;

    // PA(i:, i:)P = A(i:,i:) - (v w^T + w v^T) + beta nu v v^T
    a = a0;
    scalar_t a_i0_i0 = *a;
    a_i0_i0 -= 2*w[0];
    a_i0_i0 += beta_nu;
    *a++ = a_i0_i0;
    for (int j=i0+1; j<n; ++j) {
      scalar_t a_i0_j = *a;
      a_i0_j -= w[j-i0];
      a_i0_j -= w[0]*v[j-i0-1];
      a_i0_j += beta_nu * v[j-i0-1];
      *a++ = a_i0_j;
    }
    for (int i=i0+1; i<n; ++i) for (int j=i; j<n; ++j) {
      scalar_t a_ij = *a;
      a_ij -= v[i-i0-1]*w[j-i0];
      a_ij -= w[i-i0]*v[j-i0-1];
      a_ij += beta_nu * v[i-i0-1]*v[j-i0-1];
      *a++ = a_ij;
    }
  }

  /// Accumulate the k Householder reflection stored in factored form in
  /// the columns of the matrix \c a given the corresponding beta's
  /** \c q is filled with the product of Householder
      reflections H(0) H(1) .. H(k-1) where H(j) is stored
      in the j-th column of \c a.

      q may be rectangular, m x n, in which case q is filled with the
      first n columns of that product (useful for the thin QR).

      The argument off_diag is used to handle the both of the QR and
      the bidiagonalisation algorithms with the same code:
      it must be respectively 0 and 1 for those two cases.
  C.f. Golub and Van Loan section 5.1.6
  */
  void
  accumulate_factored_form_in_columns(matrix_ref_t const &q,
                                      matrix_const_ref_t const &a,
                                      af::const_ref<scalar_t> const &beta,
                                      int const off_diag=0)
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(q.n_rows() == m)(q.n_rows())(m); // A = Q x ...
    q.set_identity(false); // Q may be rectangular
    for (int j=beta.size()-1; j >= 0; --j) {
      for (int i=j + off_diag + 1; i < m; ++i) v[i - j - off_diag - 1] = a(i,j);
      this->beta = beta[j];
      apply_on_left_to_lower_right_block(q, j+off_diag, j+off_diag);
    }
  }

  /// Accumulate the k Householder reflection stored in factored form in the
  /// columns of the matrix a, in place into a
  /** The first n columns of the product H(0) H(1) .. H(k-1) are stored
      in the first n columns of \c a
      Reference: LAPACK DORG2R
  */
  void accumulate_in_place_factored_form_in_columns(
    matrix_ref_t const &a,
    af::const_ref<scalar_t> const &beta)
  {
    int m = a.n_rows(), n = a.n_columns();
    // Don't forget Q and A occupy the same memory locations!
    if (m <= n) {
      // the last reduced column of A does not contain any Householder vector
      // and it will therefore not be touched by the coming loop
      for (int i=0; i<m-1; ++i) a(i, m-1) = 0;
      a(m-1, m-1) = 1;
    }
    // Apply Householder reflections in order
    for (int j=beta.size()-1; j >= 0; --j) {
      // A(j+1:, j) starts as the essential part of the vector of H(j)
      // Q(j:, j) starts as [ 1 0 ... 0 ]^T
      // Copy the former to v and perform  Q(j:, j) = H(j) Q(j:, j)
      // in one loop.
      for (int i=j+1; i < m; ++i) {
        v[i-j-1] = a(i,j);
        a(i,j) *= -beta[j];
      }
      this->beta = beta[j];
      a(j,j) = 1 - beta[j];

      // Q(j:, j+1:) = H(j) Q(j:, j+1:)
      if (j < n-1) apply_on_left_to_lower_right_block(a, j, j+1);

      // Q(:j, j) = 0
      for (int i=0; i<j; ++i) a(i,j) = 0;
    }
  }

  /// Accumulate the Householder reflection stored in factored form in
  /// the rows of the matrix \c a given the corresponding beta's
  /** The product formed and stored in \c q is then
      either H(0) H(1) .. H(k-1) or H(k-1) H(k-2) .. H(0) depending
      on whether the argument \c reflection_order  is respectively
      \c product_in_row_order or \c product_in_reverse_row_order,
      where H(i) is the Householder reflection stored in the i-th row of \c a.

      q may be rectangular,

        - either m x n, then q is filled with the first m rows
          of H(k-1) H(k-2) .. H(0)  (used by thin LQ)

        - n x m, then q is filled with the first m columns
          of H(0) H(1) .. H(k-1)    (used by thin bidiagonalisation)

      The argument off_diag is used to handle the both of the LQ and
      the bidiagonalisation algorithms with the same code:
      it must be respectively 0 and 1 for those two cases.
      C.f. Golub and Van Loan section 5.1.6
  */
  void
  accumulate_factored_form_in_rows(matrix_ref_t const &q,
                                   matrix_const_ref_t const &a,
                                   af::const_ref<scalar_t> const &beta,
                                   int reflection_order,
                                   int const off_diag=0)
  {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(   reflection_order == product_in_row_order
                  || reflection_order == product_in_reverse_row_order);
    switch (reflection_order) {
      case product_in_row_order:
        // the reduction of A is done by A -> A Q
        SCITBX_ASSERT(q.n_rows() == n)(q.n_rows())(n);
        break;
      case product_in_reverse_row_order:
        // the reduction of A is done by A -> A Q^T
        SCITBX_ASSERT(q.n_columns() == n)(q.n_columns())(n);
        break;
    }
    q.set_identity(false); // Q may be rectangular
    for (int i=beta.size()-1; i >= 0; --i) {
      for (int j=i + off_diag + 1; j < n; ++j) v[j - i - off_diag - 1] = a(i,j);
      this->beta = beta[i];
      switch (reflection_order) {
        case product_in_row_order:
          apply_on_left_to_lower_right_block(q, i+off_diag, i+off_diag);
          break;
        case product_in_reverse_row_order:
          apply_on_right_to_lower_right_block(q, i+off_diag, i+off_diag);
          break;
      }
    }
  }

  /// Accumulate the k Householder reflection stored in factored form in the
  /// rows of the matrix a, in place into a
  /** The first m rows of the product H(k-1) H(k-2) .. H(0) are stored
      in the first m rows of \c a
  */
  void accumulate_in_place_factored_form_in_rows(
    matrix_ref_t const &a,
    af::const_ref<scalar_t> const &beta)
  {
    int m = a.n_rows(), n = a.n_columns();
    // Don't forget Q and A occupy the same memory locations!
    if (m >= n) {
      // the last reduced row of A does not contain any Householder vector
      // and it will therefore not be touched by the coming loop
      for (int j=0; j<n-1; ++j) a(n-1, j) = 0;
      a(n-1, n-1) = 1;
    }
    for (int i=beta.size()-1; i >= 0; --i) {
      // A(i, i+1:) starts as the essential part of the vector of H(i)
      // Q(i, i:) starts as [ 1 0 ... 0 ]^T
      // Copy the former to v and perform  Q(i, i:) = H(i) Q(i, i:)
      // in one loop.
      for (int j=i+1; j < n; ++j) {
        v[j-i-1] = a(i,j);
        a(i,j) *= -beta[i];
      }
      this->beta = beta[i];
      a(i,i) = 1 - beta[i];

      // Q(i+1:, i:) = Q(i+1:, i:) H(i)
      if (i < m-1) apply_on_right_to_lower_right_block(a, i+1, i);

      // Q(i, :i) = 0
      for (int j=0; j<i; ++j) a(i,j) = 0;
    }
  }

  /// Accumulate a normal random matrix.
  /** Q being m x n,
      - if m == n, then Q is a random normal matrix;
      - if m > n, then Q only consists of the first n columns
        of a random normal matrix;
      - if m < n, then Q only consists of the first m rows
        of a random normal matrix.

      The distribution is uniform on the set of normal matrices.
      (For mathematically minded people, it is the Haar measure).

   Reference:
    G.W. Stewart,
    The efficient generation of random orthogonal matrices
    with an application to condition estimators,
    SIAM Journal on Numerical Analysis 17 (1980), no. 3, 403-409.

    This is the method used in LAPACK test tool DLAGGE for example.
  */
  template <class UniformRandomNumberGenerator>
  void
  accumulate_random_normal_matrix(
    boost::variate_generator<UniformRandomNumberGenerator,
                             boost::normal_distribution<scalar_t> > &normal,
    matrix_ref_t const &q)
  {
    int m = q.n_rows(), n = q.n_columns();
    q.set_identity(false);
    for (int i=std::min(m, n) - 1; i >= 0; --i) {
      if (i < n-1) {
        for (int k=0; k<n-i; ++k) v[k] = normal();
        zero_vector(n-i);
        apply_on_right_to_lower_right_block(q, i, i);
      }
    }
  }

  /// Accumulate a random matrix with given singular values
  template <class UniformRandomNumberGenerator>
  void
  accumulate_random_matrix_with_singular_values(
    boost::variate_generator<UniformRandomNumberGenerator,
                             boost::normal_distribution<scalar_t> > &normal,
    af::const_ref<scalar_t> const &sigma,
    matrix_ref_t const &a)
  {
    int m = a.n_rows(), n = a.n_columns();
    a.set_diagonal(sigma, false);
    for (int i=std::min(m, n) - 1; i >= 0; --i) {
      if (i < m-1) {
        for (int k=0; k<m-i; ++k) v[k] = normal();
        zero_vector(m-i);
        apply_on_left_to_lower_right_block(a, i, i);
      }
      if (i < n-1) {
        for (int k=0; k<n-i; ++k) v[k] = normal();
        zero_vector(n-i);
        apply_on_right_to_lower_right_block(a, i, i);
      }
    }
  }

  /// Accumulate a random symmetric matrix with given eigenvalues
  template <class UniformRandomNumberGenerator>
  void
  accumulate_random_symmetric_matrix_with_eigenvalues(
    boost::variate_generator<UniformRandomNumberGenerator,
                             boost::normal_distribution<scalar_t> > &normal,
    af::const_ref<scalar_t> const &lambda,
    symmetric_matrix_packed_u_ref_t const &a)
  {
    int n = a.n_columns();
    a.set_diagonal(lambda);
    for (int i=n-2; i>=0; --i) {
      for (int k=0; k<n-i; ++k) v[k] = normal();
      zero_vector(n-i);
      apply_to_lower_right_block(a, i);
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
  typedef af::ref<scalar_t, af::mat_grid> matrix_ref_t;

  bool may_accumulate_q;
  matrix_ref_t a;
  reflection<scalar_t> p;
  std::vector<scalar_t> beta;

  /// Construct the QR decomposition of the m x n matrix A in place.
  /** The lower diagonal part of R is zero and its non-zero part is
      stored above and on the diagonal of A.

      Q is a product of Householder reflections whose vectors are stored
      below the diagonal of A.
      The argument may_accumulate_q_ must be set to true if Q is to be
      accumulated later.
  */
  qr_decomposition(matrix_ref_t const &a_, bool may_accumulate_q_=true)
    : a(a_),
      may_accumulate_q(may_accumulate_q_),
      p(a.n_rows(), a.n_columns(), applied_on_left_tag(), may_accumulate_q)
  {
    int m = a.n_rows(), n = a.n_columns(), k=(m > n ? n : m-1);
    beta.reserve(k);
    for (int j=0; j<k; ++j) {
      p.zero_vector(af::column_below(a, j, j));
      beta.push_back(p.beta);
      p.apply_on_left_to_lower_right_block(a, j, j+1);
    }
  }

  /// The Q in QR
  /** Either the full Q or only the thin Q, in a storage separate from A

        - full Q: Q is m x m and R is m x n
        - thin Q (only if m >= n): Q is m x n and R is n x n

      If \c thin is true, then the thin q is returned if m <= n
      and the full Q otherwise. If \c thin is false, the full Q is returned.
  */
  matrix_t q(bool thin=true) {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(may_accumulate_q);
    af::const_ref<scalar_t> beta_(&beta[0], beta.size());
    matrix_t q(dim(m, thin ? std::min(m,n) : m),
               af::init_functor_null<scalar_t>());
    p.accumulate_factored_form_in_columns(q.ref(), a, beta_);
    return q;
  }

  /// Accumulate the thin Q inside A in-place
  /** L is therefore lost */
  void accumulate_q_in_place() {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(may_accumulate_q);
    SCITBX_ASSERT(m >= n);
    af::const_ref<scalar_t> beta_(&beta[0], beta.size());
    p.accumulate_in_place_factored_form_in_columns(a, beta_);
  }
};

/// Householder LQ decomposition
template <typename FloatType>
struct lq_decomposition
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;
  typedef af::ref<scalar_t, af::mat_grid> matrix_ref_t;

  bool may_accumulate_q;
  matrix_ref_t a;
  reflection<scalar_t> p;
  std::vector<scalar_t> beta;

  /// Construct the LQ decomposition of A
  /** The upper diagonal part of L is zero and its non-zero part is
      stored below and on the diagonal of A.

      Q is a product of Householder reflections whose vectors are stored
      above the diagonal of A.
      The argument may_accumulate_q_ must be set to true if Q is to be
      accumulated later.
  */
  lq_decomposition(matrix_ref_t const &a_, bool may_accumulate_q_=true)
    : a(a_),
      may_accumulate_q(may_accumulate_q_),
      p(a.n_rows(), a.n_columns(), applied_on_right_tag(), may_accumulate_q)
  {
    int m = a.n_rows(), n = a.n_columns(), k=(n > m ? m : n-1);
    beta.reserve(k);
    for (int i=0; i<k; ++i) {
      p.zero_vector(af::row_right_of(a, i, i));
      beta.push_back(p.beta);
      p.apply_on_right_to_lower_right_block(a, i+1, i);
    }
  }

  /// The Q in LQ
  /** Either the full Q or only the thin Q in a storage separate from A

        - full Q: Q is n x n and L is m x n
        - thin Q (only if m <= n): Q is m x n and L is m x m

      If \c thin is true, then the thin q is returned if m <= n
      and the full Q otherwise. If \c thin is false, the full Q is returned.
  */
  matrix_t q(bool thin=true) {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(may_accumulate_q);
    af::const_ref<scalar_t> beta_(&beta[0], beta.size());
    matrix_t q(dim(thin ? std::min(m, n) : n, n),
               af::init_functor_null<scalar_t>());
    p.accumulate_factored_form_in_rows(q.ref(), a, beta_,
                                       product_in_reverse_row_order);
    return q;
  }

  /// Accumulate Q inside A in-place
  /** L is therefore lost */
  void accumulate_q_in_place() {
    int m = a.n_rows(), n = a.n_columns();
    SCITBX_ASSERT(may_accumulate_q);
    SCITBX_ASSERT(m <= n);
    af::const_ref<scalar_t> beta_(&beta[0], beta.size());
    p.accumulate_in_place_factored_form_in_rows(a, beta_);
  }
};


/// Decomposition U^T A V = B where B is bibidiagonal and U,V orthogonal
/** B is upper diagonal if m >= n and lower diagonal if m < n.
    Reference: Golub and Van Loan, Algorithm 5.4.2 */
template <typename FloatType>
struct bidiagonalisation
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;
  typedef af::ref<scalar_t, af::mat_grid> matrix_ref_t;

  matrix_ref_t a;
  reflection<FloatType> p;
  std::vector<scalar_t> beta_left, beta_right;

  /// Overwrite a in-place with its bidiagonalisation
  /** The bidiagonal is that of B whereas the rest of a stores
      the essential part of the Householder vectors making U and V
  */
  bidiagonalisation(matrix_ref_t const &a_)
    : a(a_),
      p(a.n_rows(), a.n_columns(), applied_on_left_and_right_tag())
  {
    int m = a.n_rows(), n = a.n_columns();
    if (m >= n) {
      int k_left = m > n ? n : n-1;
      int k_right = n-2;
      beta_left.reserve(k_left);
      beta_right.reserve(k_right);
      for (int j=0; j<k_left; ++j) {
        p.zero_vector(af::column_below(a, j, j));
        beta_left.push_back(p.beta);
        p.apply_on_left_to_lower_right_block(a, j, j+1);
        if (j < k_right) {
          p.zero_vector(af::row_right_of(a, j, j+1));
          beta_right.push_back(p.beta);
          p.apply_on_right_to_lower_right_block(a, j+1, j+1);
        }
      }
    }
    else {
      int k_right = m;
      int k_left = m-2;
      for (int i=0; i<k_right; ++i) {
        p.zero_vector(af::row_right_of(a, i, i));
        beta_right.push_back(p.beta);
        p.apply_on_right_to_lower_right_block(a, i+1, i);
        if (i < k_left) {
          p.zero_vector(af::column_below(a, i+1, i));
          beta_left.push_back(p.beta);
          p.apply_on_left_to_lower_right_block(a, i+1, i+1);
        }
      }
    }
  }

  /// The matrix U, either full or thin (c.f. class \c qr_decomposition)
  matrix_t u(bool thin=true) {
    int m = a.n_rows(), n = a.n_columns();
    af::const_ref<scalar_t> beta_left_(&beta_left[0], beta_left.size());
    matrix_t u(dim(m, thin ? std::min(m, n) : m),
               af::init_functor_null<scalar_t>());
    p.accumulate_factored_form_in_columns(u.ref(), a, beta_left_,
                                          m >= n ? 0 : 1);
    return u;
  }

  /// The matrix V, either full or thin (c.f. class \c lq_decomposition)
  matrix_t v(bool thin=true) {
    int m = a.n_rows(), n = a.n_columns();
    af::const_ref<scalar_t> beta_right_(&beta_right[0], beta_right.size());
    matrix_t v(dim(n, thin ? std::min(m,n) : n),
               af::init_functor_null<scalar_t>());
    p.accumulate_factored_form_in_rows(v.ref(), a, beta_right_,
                                       product_in_row_order,
                                       m >= n ? 1 : 0);
    return v;
  }
};


template <typename FloatType, class UniformRandomNumberGenerator>
struct random_normal_matrix_generator
{
  typedef FloatType scalar_t;
  typedef af::c_grid<2> dim;
  typedef af::versa<scalar_t, dim> matrix_t;
  typedef af::ref<scalar_t, af::mat_grid> matrix_ref_t;
  typedef af::versa<double, af::packed_u_accessor> symmetric_matrix_packed_u_t;
  typedef af::ref<double, af::packed_u_accessor> symmetric_matrix_packed_u_ref_t;

  UniformRandomNumberGenerator uniform_gen;
  boost::normal_distribution<scalar_t> normal_dist;

  boost::variate_generator<UniformRandomNumberGenerator,
                           boost::normal_distribution<scalar_t> > normal_gen;

  int m, n;
  reflection<scalar_t> p;

  random_normal_matrix_generator(UniformRandomNumberGenerator &uniform,
                                 int rows, int columns)
    : uniform_gen(uniform), normal_dist(0, 1),
      normal_gen(uniform_gen, normal_dist),
      m(rows), n(columns),
      p(m, n, applied_on_left_and_right_tag())
  {}

  random_normal_matrix_generator(int rows, int columns)
    : normal_dist(0, 1),
      normal_gen(uniform_gen, normal_dist),
      m(rows), n(columns),
      p(m, n, applied_on_left_and_right_tag())
  {}

  matrix_t normal_matrix() {
    matrix_t result(dim(m, n), af::init_functor_null<scalar_t>());
    matrix_ref_t q = result.ref();
    p.accumulate_random_normal_matrix(normal_gen, q);
    return result;
  }

  matrix_t matrix_with_singular_values(af::const_ref<scalar_t> const &sigma) {
    matrix_t result(dim(m, n), af::init_functor_null<scalar_t>());
    matrix_ref_t a = result.ref();
    p.accumulate_random_matrix_with_singular_values(normal_gen, sigma, a);
    return result;
  }

  symmetric_matrix_packed_u_t
  symmetric_matrix_with_eigenvalues(af::const_ref<scalar_t> const& lambda) {
    SCITBX_ASSERT(m == n)(m)(n);
    symmetric_matrix_packed_u_t result(n, af::init_functor_null<scalar_t>());
    symmetric_matrix_packed_u_ref_t a = result.ref();
    p.accumulate_random_symmetric_matrix_with_eigenvalues(normal_gen, lambda, a);
    return result;
  }

};


}}} // scitbx::matrix::householder

#endif // GUARD
