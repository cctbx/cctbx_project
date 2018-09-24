#ifndef SCITBX_MATRIX_SVD_H
#define SCITBX_MATRIX_SVD_H

#include <scitbx/error.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/math/copysign.h>
#include <scitbx/math/numeric_limits.h>
#include <scitbx/matrix/householder.h>
#include <scitbx/matrix/move.h>
#include <algorithm>
#include <functional>
#include <limits>

#include <scitbx/matrix/givens.h>

namespace scitbx { namespace matrix { namespace svd {

/// Special SVD for the 2x2 case
/** Reference: LAPACK DLASV2 */
template <typename FloatType>
struct bidiagonal_2x2_decomposition
{
  typedef FloatType scalar_t;

  /// singular values:  |s_min| < |s_max|
  /// U^T A V = [ s_max    0   ]
  ///           [   0    s_min ]
  scalar_t s_min, s_max;

  /// U = [  c_u  -s_u ]
  ///     [  s_u   c_u ]
  scalar_t c_u, s_u;

  /// V = [  c_v  -s_v ]
  ///     [  s_v   c_v ]
  scalar_t c_v, s_v;

  /// Singular values of
  ///     [  f  g  ]
  ///     [  0  h  ]
  /** If the argument \c compute_singular_vectors is false, only the
   singular values are computed and they are both non-negative (making that code
   more or less equivalent to DLAS2 in LAPACK). Otherwise, s_max and s_min have
   respectively the sign of f and g.
   */
   bidiagonal_2x2_decomposition(scalar_t f, scalar_t g, scalar_t h,
                                bool compute_singular_vectors)
   {
    using math::copysign;

    scalar_t fa = std::abs(f);
    scalar_t ga = std::abs(g);
    scalar_t ha = std::abs(h);

    // Make |f| >= |h| and remember whether we swapped
    bool swap_f_h = (ha > fa);
    if (swap_f_h) {
      std::swap(fa, ha);
      std::swap(f, h);
    }

    // diagonal matrix
    if (ga == 0) {
      s_min = ha;
      s_max = fa;
      return;
    }

    // very large |g|
    if (ga > fa && fa/ga < std::numeric_limits<scalar_t>::epsilon()) {
      s_max = ga;
      s_min = ha > 1 ? fa/(ga/ha) : (fa/ga)*ha;
      if (compute_singular_vectors) {
        c_u = 1;
        s_u = h/g;
        c_v = f/g;
        s_v = 1;
      }
      return;
    }

    // normal case
    scalar_t d = fa - ha;
    scalar_t l = d/fa; // 0 <= l <= 1
    scalar_t m = g/f; // |m| < 1/u
    scalar_t t = 2 - l; // t >= 1
    scalar_t mm = m*m, tt = t*t;
    scalar_t s = std::sqrt(tt + mm); // 1 <= s <= 1 + 1/u
    scalar_t r = l != 0 ? std::sqrt(l*l + mm) : std::abs(m); // 0 <= r <= 1 + 1/u
    scalar_t a = 0.5*(s+r); // 1 <= a <= 1 + |m|
    s_min = ha/a;
    s_max = fa*a;

    if (!compute_singular_vectors) return;

    // Compute singular vectors
    scalar_t tau; // second assignment to T in DLASV2
    if (mm != 0) {
      tau = ( m/(s+t) + m/(r+l) )*(1+a); // That's just -2/m (1-a^2)^2
    }
    else {
      // m very tiny
      tau = l == 0 ? copysign(scalar_t(2), f) * copysign(scalar_t(1), g)
                   : g/copysign(d, f) + m/t;
    }
    scalar_t lv = std::sqrt(tau*tau + 4); // second assignment to L in DLASV2
    c_v = 2/lv;
    s_v = tau/lv;
    c_u = (c_v + s_v*m)/a;
    s_u = (h/f)*s_v/a;
    /* Fix signs of the singular values to abide with the sign choices
       made to obtain the singular vectors
    */
    s_max = copysign(s_max, f);
    s_min = copysign(s_min, h);
    // Fix signs

    if (swap_f_h) {
      std::swap(c_u, s_v);
      std::swap(s_u, c_v);
    }
  }
};


enum bidiagonal_kind { upper_bidiagonal_kind, lower_bidiagonal_kind };

/// Golub-Kahan decomposition of a bidiagonal matrix B
/** See [1] for an introduction. Our code is based on [2,3,4].

    References:
    [1] Golub and Van Loan, Matrix Computations, 3rd edition, section 8.6.2
    [2] LAPACK source code (version 3.0)
    [3] Ake Bjorck, Numerical Methods for Least-Squares Problems
    [4] Accurate Singular Values of Bidiagonal Matrices, Demmel, J. and Kahan, W.
        SIAM J. Sci. Stat. Comput. 11, 873-912, 1990
*/
template <typename FloatType>
struct bidiagonal_decomposition
{
  typedef FloatType scalar_t;

  af::ref<scalar_t> d,f;
  scalar_t eps, tol, thresh;
  af::ref<scalar_t, af::mat_grid> u, v;
  givens::product<scalar_t> q_u, q_v;
  int n_iterations, n_max_iterations;
  bool has_iteration_converged, has_converged;

  // Whether the singular values and eventually vectors are sorted
  bool sorted;

  // the block B(s:, s:) is diagonal
  // the block B(r:s, r:s) has non-zero superdiagonal entries
  int r,s;
  int r0, s0; // previous values
  scalar_t s_lower, s_upper; // lower and upper bound for B(r:s, r:s)
                             // singular values

  // In which direction to chase the nonzeroes
  bool shall_chase_nonzero_down;

  // Wilkinson shift or 0
  scalar_t shift;

  /// Prepare for the decomposition of the given bidiagonal matrix
  /** The argument kind can be either upper_diagonal or lower_diagonal.

      \c u0 and \c v0 shall be the U and V matrices from the bidiagonalisation.
      Each time a Givens rotation is performed, it may be accumulated into
      either \c u0 and \c v0 to eventually produce the final U and V of the SVD.
      Whether accumulation takes place depends on the arguments
      \c accumulate_u and \c accumulate_v.
  */
  bidiagonal_decomposition(af::ref<scalar_t> const &diagonal,
                           af::ref<scalar_t> const &off_diagonal,
                           int kind,
                           af::ref<scalar_t, af::mat_grid> const &u0,
                           bool accumulate_u,
                           af::ref<scalar_t, af::mat_grid> const &v0,
                           bool accumulate_v,
                           scalar_t epsilon
                             =std::numeric_limits<scalar_t>::epsilon(),
                           int max_iteration_multiplier=6)
    : d(diagonal), f(off_diagonal),
      u(u0), v(v0),
      q_u(diagonal.size(), accumulate_u),
      q_v(diagonal.size(), accumulate_v),
      r(0), s(diagonal.size()),
      eps(epsilon),
      n_iterations(0),
      n_max_iterations(max_iteration_multiplier*diagonal.size()*diagonal.size()),
      sorted(false)
  {
    SCITBX_ASSERT(diagonal.size() >= 2);
    SCITBX_ASSERT(off_diagonal.size() == diagonal.size()-1);
    SCITBX_ASSERT(!accumulate_u || u0.n_columns() == diagonal.size());
    SCITBX_ASSERT(!accumulate_v || v0.n_columns() == diagonal.size());
    SCITBX_ASSERT(kind == upper_bidiagonal_kind || kind == lower_bidiagonal_kind);

    int n = diagonal.size();

    // Tolerance
    tol = epsilon;
    tol *= std::max(10., std::min(100., std::pow(epsilon, -1./8)));

    // if the matrix is lower diagonal, make it upper diagonal
    // with a series of Givens rotation on the left
    if (kind == lower_bidiagonal_kind) {
      givens::rotation<scalar_t> g;
      for (int i=0; i<n-1; ++i) {
        g.chase_nonzero_from_x1_to_y0(d[i]  , f[i],
                                      f[i], d[i+1]);
        q_u.multiply_by(g);
      }
      q_u.apply_downward_on_right(u, 0);
    }


    // B's singular values are bounded below by s_lower
    // c.f. [3] and [4, eq.(2.4)] and related
    scalar_t s_lower;
    scalar_t mu = std::abs(d[0]);
    if (mu != 0) {
      s_lower = mu;
      for (int i=1; i<n; ++i) {
        mu = std::abs(d[i])*( mu/(mu + std::abs(f[i-1])) );
        s_lower = std::min(s_lower, mu);
        if (mu == 0) break;
      }
    }
    else {
      s_lower = 0;
    }
    s_lower /= std::sqrt(scalar_t(n));

    scalar_t sfmin = math::numeric_limits<scalar_t>::safe_min();
    thresh = std::max(tol*s_lower, n_max_iterations*sfmin);
  }

  /// Compute the decomposition
  /**
      Reference: [2, DBDSQR] and [3, 2.6.4]
      A few important points:
      - we handle zeroes on the diagonal with Givens rotation theta = pi/2
        as in [2, DBDSQR] (that alleviates the need of the decoupling stages
        described in [1]);
      - the Wilkinson shift is computed as a singular value of the
        trailing 2x2 block of B(r:s, r:s), as advised in [3, 2.6.4]
        and implemented in [2, DLAS2] (instead of an eigenvalue of
        the trailing 2x2 block of B(r:s, r:s)^T B(r:s, r:s)
        in the original Golub-Kahan algorithm described in [1]).
      - the shift of the leading 2x2 block of B(r:s, r:s)^T B(r:s, r:s)
        is carefully done to avoid overflow and underflow [2, DBDSQR, line 577]
      - when B(r:s, r:s) is a 2x2 matrix, we compute its SVD decomposition
        analytically with bidiagonal_2x2_decomposition
        (as in [2, DBDSQR] with DLASV2)
      - we check the direction of grading on B(r:s, r:s) diagonal
        and apply the algorithm from top to bottom or from bottom to top,
        depending on which end is the biggest ([2] and [3, 2.6.4]).
        This is essential for efficient convergence (or convergence at all
        when there are zero entries on the diagonal).
  */
  void compute() {
    int n = d.size();

    // Main iteration
    while (s > 1 && n_iterations < n_max_iterations) {
      // Find diagonal block B2 to work on
      s_upper = std::abs(d[s-1]);
      r = 0;
      for (int i=s-2; i>=0; --i) {
        scalar_t abs_f = std::abs(f[i]);
        if (abs_f <= thresh) {
          f[i] = 0;
          r = i+1;
          break;
        }
        s_upper = std::max(s_upper, std::abs(d[i]));
        s_upper = std::max(s_upper, abs_f);
      }

      // Current B2 is done?
      if (s-r == 1) {
        s--;
        continue;
      }

      // 2x2 case?
      if (s-r == 2) {
        solve_2x2_case();
        s -= 2;
        continue;
      }

      // if B2 disjoint with previous one, ...
      if (n_iterations == 0 || r >= s0 || s <= r0) {
        // ... evaluate grading
        shall_chase_nonzero_down = std::abs(d[r]) >= std::abs(d[s-1]);
      }

      // Test convergence
      if (shall_chase_nonzero_down) {
        test_downward_iteration_convergence();
      }
      else {
        test_upward_iteration_convergence();
      }
      if (has_iteration_converged) continue;

      // Compute shift
      if (n*tol*(s_lower/s_upper) <= std::max(eps, 0.01*tol)) {
        // shifting would ruin relative accuracy
        shift = 0;
      }
      else {
        scalar_t shifted;

        // Wilkinson shift obtained from singular value
        // of 2x2 block at either end
        if (shall_chase_nonzero_down) {
          compute_trailing_wilkinson_shift();
          shifted = std::abs(d[r]);
        }
        else {
          compute_leading_wilkinson_shift();
          shifted = std::abs(d[s-1]);
        }

        // Set shift to zero if negligible
        if (shifted > 0 && std::pow(shift/shifted, 2) < eps) shift = 0;
      }

      // QR iteration
      if (shall_chase_nonzero_down) {
        if (shift == 0) {
          do_implicit_zero_shift_qr_iteration_downward();
        }
        else {
          do_implicit_shift_qr_iteration_downward();
        }
        if (std::abs(f[s-2]) <= thresh) f[s-2] = 0;
      }
      else {
        if (shift == 0) {
          do_implicit_zero_shift_qr_iteration_upward();
        }
        else {
          do_implicit_shift_qr_iteration_upward();
        }
        if (std::abs(f[r]) <= thresh) f[r] = 0;
      }

      // Done with this iteration
      r0 = r;
      s0 = s;
      n_iterations += s - r;
    }

    has_converged = s <= 1;

    // Make all singular values non-negative ...
    for (int j=0; j<n; ++j){
      if (d[j] >= 0) continue;
      d[j] = -d[j];
      if (q_v.effective) {
        // B = U Sigma V^T: so change sign of i-th row of V^T
        for (int i=0; i<v.n_rows(); ++i) v(i,j) = -v(i,j);
      }
    }
  }

  void test_downward_iteration_convergence() {
    has_iteration_converged = false;

    // quick test at bottom
    if (std::abs(f[s-2]) <= tol*std::abs(d[s-1])) {
      f[s-2] = 0;
      has_iteration_converged = true;
      return;
    }
    // thorough one applied forward
    scalar_t mu = std::abs(d[r]);
    s_lower = mu;
    for (int i=r; i<s-1; ++i) {
      if (std::abs(f[i]) <= tol*mu) {
        f[i] = 0;
        has_iteration_converged = true;
        return;
      }
      mu = std::abs(d[i+1])*( mu/(mu + std::abs(f[i])));
      s_lower = std::min(s_lower, mu);
    }
  }

  void test_upward_iteration_convergence() {
    has_iteration_converged = false;

    // quick test at top
    if (std::abs(f[r]) <= tol*std::abs(d[r])) {
      f[r] = 0;
      has_iteration_converged = true;
      return;
    }
    // thorough one applied backward
    scalar_t mu = std::abs(d[s-1]);
    s_lower = mu;
    for (int i=s-2; i>=r; --i) {
      if (std::abs(f[i]) <= tol*mu) {
        f[i] = 0;
        has_iteration_converged = true;
        return;
      }
      mu = std::abs(d[i])*( mu/(mu + std::abs(f[i])));
      s_lower = std::min(s_lower, mu);
    }
  }

  void solve_2x2_case() {
    bidiagonal_2x2_decomposition<scalar_t> svd(d[r], f[r], d[r+1], true);
    d[r]   = svd.s_max;
    d[r+1] = svd.s_min;
    f[r]   = 0;
    // same signs convention for sines in 2x2 decomposition and givens
    givens::rotation<scalar_t> g_u(svd.c_u, svd.s_u),
                               g_v(svd.c_v, svd.s_v);
    if (q_u.effective) g_u.apply_on_right(u, r, r+1);
    if (q_v.effective) g_v.apply_on_right(v, r, r+1);
  }

  void compute_trailing_wilkinson_shift() {
    bidiagonal_2x2_decomposition<scalar_t> trailing_svd(
                d[s-2], f[s-2], d[s-1], false);
    shift = trailing_svd.s_min;
  }

  void compute_leading_wilkinson_shift() {
    bidiagonal_2x2_decomposition<scalar_t> leading_svd(
                d[r], f[r], d[r+1], false);
    shift = leading_svd.s_min;
  }

  /// Do one Golub-Kahan iteration chasing nonzeroes from top-left to bottom-right
  void do_implicit_shift_qr_iteration_downward(bool compute_shift=false) {
    using math::copysign;

    if (compute_shift) compute_trailing_wilkinson_shift();

    /* Givens rotation zeroing off-diagonal of leading 2x2 block
       of S = B^T B - shift I (implicit QR shift)
       For numerical stability, S is scaled by d[r].
    */
    scalar_t s00 = (std::abs(d[r]) - shift)*(copysign(1., d[r]) + shift/d[r]);
    scalar_t s01 = f[r];
    givens::rotation<scalar_t> g1(r, r+1);
    g1.zero_x1(s00, s01);

    // Apply G1 on the right to B, filling in b(r+1, r), so we now need...
    scalar_t z;
    g1.apply(d[r], f[r]);
    g1.apply_assuming_null_x0(z, d[r+1]); // z = b(r+1, r)
    q_v.multiply_by(g1);

    // ... to chase non-zeroes down the diagonal
    for (int k=r; k < s-2; ++k) {
      scalar_t t;
      givens::rotation<scalar_t> g_u, g_v;
      g_u.chase_nonzero_from_x1_to_z0(d[k], f[k]  , t,
                                      z   , d[k+1], f[k+1]);
      z = t;
      q_u.multiply_by(g_u);
      int l = k+1;
      g_v.chase_nonzero_from_x1_to_z0(f[l-1], d[l], t,
                                      z     , f[l], d[l+1]);
      z = t;
      q_v.multiply_by(g_v);
    }
    givens::rotation<scalar_t> g_u;
    g_u.chase_nonzero_from_x1_off(d[s-2], f[s-2],
                                  z     , d[s-1]);
    q_u.multiply_by(g_u);

    // Accumulate
    q_u.apply_downward_on_right(u, r);
    q_v.apply_downward_on_right(v, r);
  }

  /// Do one Golub-Kahan iteration chasing nonzeroes from bottom-right to top-left
  void do_implicit_shift_qr_iteration_upward(bool compute_shift=false) {
    using math::copysign;

    if (compute_shift) compute_leading_wilkinson_shift();

    /* Givens rotation zeroing off-diagonal of leading 2x2 block
       of S = B^T B - shift I (implicit QR shift)
       For numerical stability, S is scaled by d[r].
    */
    scalar_t s00 = (std::abs(d[s-1]) - shift)*(copysign(1., d[s-1]) + shift/d[s-1]);
    scalar_t s01 = f[s-2];
    givens::rotation<scalar_t> g1;
    g1.zero_x1(s00, s01);

    // Apply G1 on the left to B, filling in b(s-1, s-2), so we now need...
    scalar_t z;
    g1.apply(d[s-1], f[s-2]);
    g1.apply_assuming_null_x0(z, d[s-2]); // z = b(s-1, s-2)
    q_u.multiply_by(g1);

    // ... to chase non-zeroes up the diagonal
    for (int k=s-1; k >= r+2; --k) {
      scalar_t t;
      givens::rotation<scalar_t> g_v, g_u;
      g_v.chase_nonzero_from_x1_to_z0(d[k], f[k-1], t,
                                       z  , d[k-1], f[k-2]);
      z = t;
      q_v.multiply_by(g_v);
      int l = k-1;
      g_u.chase_nonzero_from_x1_to_z0(f[l], d[l]  , t,
                                      z   , f[l-1], d[l-1]);
      z = t;
      q_u.multiply_by(g_u);
    }
    givens::rotation<scalar_t> g_v;
    g_v.chase_nonzero_from_x1_off(d[r+1], f[r],
                                   z    , d[r]);
    q_v.multiply_by(g_v);

    // Accumulate
    q_u.apply_upward_on_right(u, s-1);
    q_v.apply_upward_on_right(v, s-1);
  }

  /// Do one Demmel-Kahan iteration,
  /// chasing nonzeroes from top-left to bottom-right
  void do_implicit_zero_shift_qr_iteration_downward(bool compute_shift=false) {
    if (compute_shift) shift = 0;

    // Givens rotation zeroing off-diagonal of leading 2x2 block of S=B^t B
    givens::rotation<scalar_t> g1;
    g1.zero_x1(d[r], f[r]);

    // Apple G1 on the left to B, filling in b(r+1, r), so we need...
    scalar_t z;
    g1.apply_assuming_null_x0(z, d[r+1]);
    q_v.multiply_by(g1);

    // ... to chase non-zeroes down the diagonal
    for (int k=r; k < s-2; ++k) {
      scalar_t t;
      givens::demmel_kahan_rotations<scalar_t> dk;
      dk.chase_non_zero_from_z_to_t(d[k], f[k]  ,
                                    z   , d[k+1], f[k+1],
                                          t     , d[k+2]);
      q_u.multiply_by(dk.g1);
      q_v.multiply_by(dk.g2);
      z = t;
    }
    givens::rotation<scalar_t> g_u;
    g_u.chase_nonzero_from_x1_to_y0(d[s-2], f[s-2],
                                    z     , d[s-1]);
    q_u.multiply_by(g_u);

    // Accumulate
    q_u.apply_downward_on_right(u, r);
    q_v.apply_downward_on_right(v, r);
  }

  /// Do one Demmel-Kahan iteration,
  /// chasing nonzeroes from bottom-right to top-left
  void do_implicit_zero_shift_qr_iteration_upward(bool compute_shift=false) {
    if (compute_shift) shift = 0;

    // Givens rotation zeroing off-diagonal of leading 2x2 block of S=B^t B
    givens::rotation<scalar_t> g1;
    g1.zero_x1(d[s-1], f[s-2]);

    // Apply G1 on the right to B, filling in b(s-1, s-2), so we need...
    scalar_t z;
    g1.apply_assuming_null_x0(z, d[s-2]);
    q_u.multiply_by(g1);

    // ... to chase non-zeroes down the diagonal
    for (int k=s-1; k >= r+2; --k) {
      scalar_t t;
      givens::demmel_kahan_rotations<scalar_t> dk;
      dk.chase_non_zero_from_z_to_t(d[k], f[k-1],
                                    z   , d[k-1], f[k-2],
                                          t     , d[k-2]);
      q_v.multiply_by(dk.g1);
      q_u.multiply_by(dk.g2);
      z = t;
    }
    givens::rotation<scalar_t> g_v;
    g_v.chase_nonzero_from_x1_to_y0(d[r+1], f[r],
                                    z     , d[r]);
    q_v.multiply_by(g_v);

    // Accumulate
    q_u.apply_upward_on_right(u, s-1);
    q_v.apply_upward_on_right(v, s-1);
  }

  /// Sort singular values in descending order
  /** Implementation note: we use selection sort if U or V is accumulated.

      Any sort algorithm needs to swap values which are out-of-order.
      Here it means swapping not only two singular values
      but also the associated singular vectors. Thus the cost of one swap
      is 2n ops, compared to the cost of one comparison which is 1 op.
      Since selection sort reaches the theoretical minimum number of swaps (n)
      while performing n(n+1)/2 comparisons, it is the algorithm of choice here.
  */
  void sort() {
    if (sorted) return;
    int n = d.size();
    if (!q_u.effective && !q_v.effective) {
      // No accumulation => best sorting algorithm available
      std::sort(d.begin(), d.end(), std::greater<scalar_t>());
    }
    else {
      // Accumulation => selection sort
      for (int i=0; i<n; ++i) {
        scalar_t *p = std::max_element(&d[i], d.end());
        if (p > &d[i]) {
          std::swap(*p, d[i]);
          if (q_u.effective) u.swap_columns(p - &d[0], i);
          if (q_v.effective) v.swap_columns(p - &d[0], i);
        }
      }
    }
    sorted = true;
  }

  /// Numerical rank
  /** A matrix A is said to have numerical delta-rank equal to k if

        k = min { rank(B) | ||A-B||_2 <= delta }

      After ordering the singular values in decreasing order, we can use the
      following result
        sigma[1] >= ... >= sigma[k] >= delta >= sigma[k+1] >= ... >= sigma[n]
  */
  std::size_t numerical_rank(scalar_t delta) {
    sort();
    scalar_t *p = std::upper_bound(d.begin(), d.end(), delta,
                                   std::greater_equal<scalar_t>());
    return p - d.begin();
  }
};


/// Reconstruct a matrix from a thin SVD decomposition: A = U S V^T
/** The matrices dimensions are respectively:

    - U : m x p
    - S : p x p
    - V : n x p
    - A : m x n
*/
template <typename T>
af::versa<T, af::c_grid<2> >
reconstruct(af::const_ref<T, af::mat_grid> const &u,
            af::const_ref<T, af::mat_grid> const &v,
            af::const_ref<T> const &sigma)
{
  int m=u.n_rows(), p=sigma.size(), n=v.n_rows();
  SCITBX_ASSERT(u.n_columns() == p);
  SCITBX_ASSERT(v.n_columns() == p);
  typedef af::c_grid<2> dim;
  typedef af::versa<T, dim> matrix_t;
  typedef af::ref<T, dim> matrix_ref_t;
  matrix_t result(dim(m,n));
  matrix_ref_t a = result.ref();
  for (int i=0; i<m; ++i)
  for (int j=0; j<n; ++j) {
    T a_ij = 0;
    for (int k=0; k<p; ++k) a_ij += sigma[k]*u(i,k)*v(j,k);
    a(i,j) = a_ij;
  }
  return result;
}

template <typename FloatType>
struct decompose {
  typedef af::c_grid<2> dim;
  typedef af::versa<FloatType, dim> matrix_t;
  typedef af::ref<FloatType, af::mat_grid> matrix_ref_t;
  typedef af::const_ref<FloatType, af::mat_grid> matrix_const_ref_t;

  matrix_t u, v;
  af::shared<FloatType> sigma;
  bool has_u, has_v;

  decompose(matrix_ref_t const &a,
    FloatType crossover=static_cast<FloatType>(5./3),
    bool accumulate_u=false, bool accumulate_v=false)
    : has_u(accumulate_u), has_v(accumulate_v)
  {
    using namespace scitbx::matrix;
    using namespace scitbx::matrix::householder;
    int n_cols = a.n_columns(),
      n_rows = a.n_rows();
    if (n_rows > n_cols*crossover ||
      n_cols > n_rows*crossover)
    {
      matrix_t t, q;
      bool taller = n_rows > n_cols*crossover;
      if (taller) {
        qr_decomposition<FloatType> qr(a, accumulate_u);
        t = copy_upper_triangle<FloatType>(a);
        if (accumulate_u) {
          qr.accumulate_q_in_place();
          q = af::mat_const_ref_as_versa<FloatType>(a);
        }
      }
      else {
        lq_decomposition<FloatType> lq(a, accumulate_v);
        t = matrix::copy_lower_triangle<FloatType>(a);
        if (accumulate_v) {
          lq.accumulate_q_in_place();
          q = af::mat_const_ref_as_versa<FloatType>(a);
        }
      }
      bidiagonalisation<FloatType> bdg(t.ref());
      if (accumulate_u) {
        u = bdg.u();
      }
      if (accumulate_v) {
        v = bdg.v();
      }
      std::pair<af::shared<FloatType>, af::shared<FloatType> > df =
        af::matrix_upper_bidiagonal<FloatType>(t.ref());
      bidiagonal_decomposition<FloatType> svd(
        df.first.ref(), df.second.ref(), svd::upper_bidiagonal_kind,
        u.ref(), accumulate_u, v.ref(), accumulate_v);
      svd.compute();
      SCITBX_ASSERT(svd.has_converged);
      svd.sort();
      if (taller) {
        if (accumulate_u) {
          u = af::matrix_multiply<FloatType>(q.ref(), u.ref());
        }
      }
      else if (accumulate_v) {
        v = af::matrix_transpose_multiply<FloatType>(q.ref(), v.ref());
      }
      sigma = df.first;
    }
    else {
      bidiagonalisation<FloatType> bdg(a);
      if (accumulate_u) {
        u = bdg.u();
      }
      if (accumulate_v) {
        v = bdg.v();
      }
      int m_kind;
      std::pair<af::shared<FloatType>, af::shared<FloatType> > df;
      if (n_rows >= n_cols) {
        df = af::matrix_upper_bidiagonal<FloatType>(a);
        m_kind = upper_bidiagonal_kind;
      }
      else {
        df = af::matrix_lower_bidiagonal<FloatType>(a);
        m_kind = lower_bidiagonal_kind;
      }
      bidiagonal_decomposition<FloatType> svd(
        df.first.ref(), df.second.ref(), m_kind,
        u.ref(), accumulate_u, v.ref(), accumulate_v);
      svd.compute();
      SCITBX_ASSERT(svd.has_converged);
      svd.sort();
      sigma = df.first;
    }
  }

  matrix_t getU() const {
    SCITBX_ASSERT(has_u);
    return u;
  }
  matrix_t getV() const {
    SCITBX_ASSERT(has_v);
    return v;
  }
  af::shared<FloatType> getSigma() const {
    return sigma;
  }

  matrix_t reconstruct() const {
    SCITBX_ASSERT(has_u && has_v);
    return svd::reconstruct<FloatType>(
      u.const_ref(), v.const_ref(), sigma.const_ref());
  }

  std::size_t numerical_rank(FloatType delta) {
    FloatType *p = std::upper_bound(sigma.begin(), sigma.end(), delta,
      std::greater_equal<FloatType>());
    return p - sigma.begin();
  }
};

}}} // scitbx::matrix::svd

#endif // GUARD
