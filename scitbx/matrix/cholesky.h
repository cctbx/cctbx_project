#ifndef SCITBX_MATRIX_CHOLESKY_H
#define SCITBX_MATRIX_CHOLESKY_H

#include <scitbx/matrix/move.h>
#include <scitbx/matrix/triangular_systems.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <vector>
#include <algorithm>
#include <boost_adaptbx/floating_point_exceptions.h>
#include <fast_linalg/lapacke.h>

namespace scitbx { namespace matrix { namespace cholesky {

  /** @name Cholesky decomposition in place
      Those classes use the same techniques as LAPACK routine DPPTRF.
     The L L^T and U^T U cases are reversed since we use
     a packing of symmetric matrices by rows (whereas LAPACK assumes a packing
     by columns of course).
   */
  //@{

  /// Info about a failure during a Cholesky process
  template <typename FloatType>
  struct failure_info
  {
    failure_info() : failed(false) {}
    failure_info(int i, FloatType v) : failed(true), index(i), value(v) {}
    operator bool() {
      return failed;
    }
    int index;
    FloatType value;
    bool failed;
  };


  /// Umbrella for A x = b solvers in place
  struct solve_in_place
  {
    /// Solve using A = L L^T
    template <typename FloatType>
    static
    void using_l_l_transpose(af::const_ref<FloatType,
                                           af::packed_l_accessor> const &l,
                             af::ref<FloatType> const &b)
    {
      SCITBX_ASSERT(l.n_columns() == b.size());
      // x = L^{-T} L^{-1} b
      forward_substitution(b.size(), l.begin(), b.begin());
      back_substitution_given_transpose(b.size(), l.begin(), b.begin());
    }

    /// Solve using A = U^T U
    template <typename FloatType>
    static
    void using_u_transpose_u(af::const_ref<FloatType,
                                           af::packed_u_accessor> const &u,
                             af::ref<FloatType> const &b)
    {
      SCITBX_ASSERT(u.n_columns() == b.size());
      // x = U^{-1} U^{-T} b
      forward_substitution_given_transpose(b.size(), u.begin(), b.begin());
      back_substitution(b.size(), u.begin(), b.begin());
    }
  };


  /// Inverse of U^T U
  /** If LAPACKE is available (through module fast_linalg), then use XPPTRI.
      Otherwise, use an implementation of the alternative method presented
      at the end of section 2.8.3 in Ake Bjorck's classic book. Never use the
      latter in production code where speed matters because it is horrendously
      inefficient.
   */
  template <typename FloatType>
  af::versa<FloatType, af::packed_u_accessor>
  inverse_of_u_transpose_u(af::ref<FloatType, af::packed_u_accessor> const &u) {
    typedef FloatType f_t;
    if (fast_linalg::is_initialised()) {
      // Usual trick: U^T U = L L^T where L = U^T is stored columnwise
      // since U is stored rowwise
      using namespace fast_linalg;
      unsigned n = u.accessor().n;
      af::versa<f_t, af::packed_u_accessor> result(
        n, af::init_functor_null<f_t>());
      af::shared<f_t> l_rfp(n*(n+1)/2, af::init_functor_null<f_t>());
      tpttf(LAPACK_COL_MAJOR, 'N', 'L', n, u.begin(), l_rfp.begin());
      {
        using namespace boost_adaptbx::floating_point;
        exception_trapping guard(exception_trapping::dont_trap);
        lapack_int info =
          pftri(LAPACK_COL_MAJOR, 'N', 'L', n, l_rfp.begin());
        SCITBX_ASSERT(!info)(info);
      }
      tfttp(LAPACK_COL_MAJOR, 'N', 'L', n, l_rfp.begin(), result.begin());
      return result;
    }
    else {
      af::versa<f_t, af::packed_u_accessor> result(u.accessor(),
                                                 af::init_functor_null<f_t>());
      af::ref<f_t, af::packed_u_accessor> c = result.ref();

      int n = u.accessor().n;
      for (int k=n-1; k>=0; --k) {
        // formula (2.8.13)
        c(k,k) = 1/u(k,k);
        for (int j=k+1; j<n; ++j) c(k,k) -= u(k,j)*c(k,j);
        c(k,k) *= 1/u(k,k);

        // formula (2.8.14)
        for (int i=k-1; i>=0; --i) {
          c(i,k) = 0;
          for (int j=i+1; j<=k; ++j) c(i,k) += u(i,j)*c(j,k);
          for (int j=k+1; j<n; ++j) c(i,k) += u(i,j)*c(k,j);
          c(i,k) *= -1/u(i,i);
        }
      }
      return result;
    }
  }

  /// Cholesky decomposition A = L L^T in place
  template <typename FloatType>
  struct l_l_transpose_decomposition_in_place
  {
    typedef FloatType scalar_t;
    typedef af::packed_l_accessor accessor_type;
    typedef af::ref<scalar_t, accessor_type> matrix_l_ref;

    /// Information on failure
    failure_info<scalar_t> failure;

    /// L
    matrix_l_ref l;

    /// Compute the decomposition of the matrix A stored at a
    /**
     On entry, a refers to the lower diagonal of the symmetric matrix A
     packed by rows.
     On exit, it refers to the lower diagonal matrix L packed by rows.
    */
    l_l_transpose_decomposition_in_place(matrix_l_ref const &a_)
      : l(a_)
    {
      scalar_t *a = a_.begin();
      int n = a_.n_rows();
      scalar_t *ll = a; // L(:i, :i) throughout
      for (int i=0; i<n; ++i, ++a) {
        // Solve L(:i,:i)  L(i, :i)^T = A(i, :i)^T in place
        if (i > 0) forward_substitution(i, ll, a);

        // Compute L(i,i), checking for positive-definitiveness
        scalar_t *l = a; // L(i, :i)
        a += i;
        scalar_t a_ii = *a;
        scalar_t l_sq = af::sum_sq(af::ref<scalar_t>(l, i)); // ||L(i, :i)||^2
        scalar_t l_ii_sq = a_ii - l_sq;
        if (l_ii_sq <= 0) {
          failure = failure_info<scalar_t>(i, l_ii_sq);
          return;
        }
        scalar_t l_ii = std::sqrt(l_ii_sq);
        *a = l_ii;
      }
    }

    /// Solve A x = b in place
    void solve_in_place(af::ref<scalar_t> const &b) {
      solve_in_place::using_l_l_transpose(l, b);
    }

    /// Solve A x = b without overwriting b
    af::shared<scalar_t> solve(af::ref<scalar_t> const &b) {
      af::shared<scalar_t> result(b.begin(), b.end());
      solve_in_place(result.ref());
      return result;
    }
  };


  /// Cholesky decomposition A = U^T U in place
  /** If module fast_linalg is available, use  LAPACK dpftrf and consort.
      Otherwise, use a handcrafted implementation of algorithm 4.2.2
      in Golub and Van Loan. The former is orders of magnitude faster than
      the latter.
  */
  template <typename FloatType>
  struct u_transpose_u_decomposition_in_place
  {
    typedef FloatType scalar_t;
    typedef af::packed_u_accessor accessor_type;
    typedef af::ref<scalar_t, accessor_type> matrix_u_ref;

    /// Information on failure
    failure_info<scalar_t> failure;

    /// U
    matrix_u_ref u;

    /// Compute the decomposition of the matrix A stored at a
    /**
     On entry, a refers to the upper diagonal of the symmetric matrix A
     packed by rows.
     On exit, it refers to the upper diagonal matrix U packed by rows.
     */
    u_transpose_u_decomposition_in_place(matrix_u_ref const &a_)
      : u(a_)
    {
      if (fast_linalg::is_initialised()) {
        using namespace fast_linalg;
        unsigned n = u.accessor().n;
        af::shared<scalar_t> l(u.size(), af::init_functor_null<scalar_t>());
        tpttf(LAPACK_COL_MAJOR, 'N', 'L', n, a_.begin(), l.begin());
        lapack_int info = pftrf(LAPACK_COL_MAJOR, 'N', 'L', n, l.begin());
        SCITBX_ASSERT(info >= 0);
        if (info > 0) {
          failure = failure_info<scalar_t>(info, 0);
        }
        tfttp(LAPACK_COL_MAJOR, 'N', 'L', n, l.begin(), u.begin());
      }
      else {
        scalar_t *a = a_.begin();
        int n = a_.n_columns();
        for (int i = 0; i < n; ++i) {
          /// Compute U(i,i) while checking positive-definiteness
          scalar_t a_ii = *a;
          if (a_ii <= 0) {
            failure = failure_info<scalar_t>(i, a_ii);
            return;
          }
          a_ii = std::sqrt(a_ii);
          *a++ = a_ii;

          /// Compute U(i, i+1:)
          scalar_t *u = a;
          for (int j = i + 1; j < n; ++j) {
            scalar_t &a_ij = *a;
            *a++ = a_ij / a_ii;
          }

          symmetric_packed_u_rank_1_update(n - i - 1, a, u, scalar_t(-1));
        }
      }
    }

    /// Solve A x = b in place
    void solve_in_place(af::ref<scalar_t> const &b) {
      solve_in_place::using_u_transpose_u(u, b);
    }

    /// Solve A x = b without overwriting b
    af::shared<scalar_t> solve(af::ref<scalar_t> const &b) {
      af::shared<scalar_t> result(b.begin(), b.end());
      solve_in_place(result.ref());
      return result;
    }
  };

  //@}


  //! Computes PAP^{T} + E = LDL^{T}
  /** P is classically represented as a product E_{n-1} ... E_2 E_1
      where E_k swaps rows k and pivots[k]. Thus pivots[k] is the index
      of the chosen pivot in column k at the k-th step of the algorithm.
   */
  /*! P. E. Gill, W. Murray, and M. H. Wright:
      Practical Optimization.
      New York: Academic Press, 1981.

      See also:

        J. Nocedal and S. Wright:
        Numerical Optimization.
        Springer, New York, 1999, pp. 145-150.

        http://lib.stat.cmu.edu/S/glm
        file: dmcdc.r
   */
  template <typename FloatType=double, typename PivotType=std::size_t>
  struct gill_murray_wright_decomposition_in_place
  {
    FloatType epsilon;
    af::shared<FloatType> packed_u;
    af::shared<FloatType> e;
    af::shared<PivotType> pivots;

    gill_murray_wright_decomposition_in_place() {}

    gill_murray_wright_decomposition_in_place(
      af::shared<FloatType> const& packed_u_,
      FloatType epsilon_=0)
    :
      epsilon(epsilon_),
      packed_u(packed_u_)
    {
      typedef FloatType f_t;
      static FloatType
        floating_point_epsilon = math::floating_point_epsilon<f_t>::get();
      if (epsilon <= 0) epsilon = floating_point_epsilon;
      af::ref<f_t> u = packed_u.ref();
      unsigned n = af::dimension_from_packed_size(u.size());
      e.resize(n);
      pivots.resize(n);
      f_t gamma = 0;
      f_t chi = 0;
      unsigned ij = 0;
      for(unsigned i=0;i<n;i++) {
        gamma = std::max(gamma, fn::absolute(u[ij++]));
        for(unsigned j=i+1;j<n;j++) {
          chi = std::max(chi, fn::absolute(u[ij++]));
        }
      }
      f_t delta = epsilon * std::max(gamma+chi, static_cast<f_t>(1));
      f_t beta_sq = std::max(gamma, epsilon);
      if (n > 1) {
        beta_sq = std::max(beta_sq, chi/std::sqrt(static_cast<f_t>(n*n-1)));
      }
      unsigned jj = 0;
      for(unsigned j=0;j<n;j++) {
        // pivoting
        {
          unsigned jmax = j;
          unsigned ii = jj;
          f_t uii_max = fn::absolute(u[ii]);
          ii += (n-j);
          for(unsigned i=j+1;i<n;ii+=(n-i),i++) {
            f_t uii = fn::absolute(u[ii]);
            if (uii_max < uii) {
                uii_max = uii;
                jmax = i;
            }
          }
          if (jmax != j) {
            packed_u_swap_rows_and_columns_in_place(u, j, jmax);
          }
          pivots[j] = jmax;
        }
        // compute j-th column of L^{T}
        {
          unsigned ii = 0;
          unsigned ij = j;
          for(unsigned i=0;i<j;ii+=(n-i),i++,ij+=(n-i)) {
            u[ij] /= u[ii];
          }
        }
        // update j-th row
        {
          unsigned kj = j;
          for(unsigned k=0;k<j;k++,kj+=(n-k)) {
            f_t f = u[kj];
            unsigned ki = kj;
            unsigned ji = jj;
            for(unsigned i=j;++i<n;) {
              u[++ji] -= u[++ki] * f;
            }
          }
        }
        // specify theta
        f_t theta_sq;
        if (j+1 == n) {
          theta_sq = 0;
        }
        else {
          unsigned jk = jj;
          theta_sq = fn::absolute(u[++jk]);
          for(unsigned k=j+2;k<n;k++) {
            f_t aajk = fn::absolute(u[++jk]);
            if (theta_sq < aajk) {
                theta_sq = aajk;
            }
          }
          theta_sq *= theta_sq;
        }
        // compute diagonal (j,j)
        {
          f_t f =
            std::max(delta, std::max(fn::absolute(u[jj]), theta_sq/beta_sq));
          e[j] = f - u[jj];
          u[jj] = f;
          // update remaining diagonals
          unsigned ji = jj;
          jj += (n-j);
          unsigned ii = jj;
          for(unsigned i=j+1;i<n;ii+=(n-i),i++) {
            u[ii] -= fn::pow2(u[++ji]) / f;
          }
        }
      }
      // scale
      ij = 0;
      for(unsigned i=0;i<n;i++) {
        f_t ujj = u[ij] = std::sqrt(u[ij]);
        ij++;
        for(unsigned j=i+1;j<n;j++) {
          u[ij++] *= ujj;
        }
      }
    }

    af::shared<FloatType>
    solve(af::const_ref<FloatType> const& b) const
    {
      int n = pivots.size();
      af::const_ref<FloatType, af::packed_u_accessor> u(packed_u.begin(), n);
      // x = P^T (U^T U)^{-1} P b
      af::shared<FloatType> result(b.begin(), b.end());
      permutation_vector(n, result.begin(), pivots.begin());
      solve_in_place::using_u_transpose_u(u, result.ref());
      permutation_transposed_vector(n, result.begin(), pivots.begin());
      return result;
    }
  };

}}} // namespace scitbx::matrix::cholesky

#endif // SCITBX_MATRIX_CHOLESKY_H
