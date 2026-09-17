#pragma once
#include <smtbx/error.h>
#include <fast_linalg/lapacke.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

namespace smtbx { namespace ED {

  /** @brief The working storage one eigendecomposition needs.

  Eight small arrays, and at eleven beams every one of them is a few hundred
  bytes -- so the allocation is a real fraction of the call, and there are tens
  of thousands of calls per pass. Held by the caller and reused, they are
  allocated once for a whole beam group instead.

  Sizing is done here rather than at each use so that the arrays grow to the
  largest matrix seen and then stop; resize() on a vector already big enough
  does nothing, which is the common case after the first call.
  */
  template <typename FloatType>
  struct hermitian_eigen_scratch {
    std::vector<std::complex<FloatType> > t, q, v, p, sub;
    std::vector<FloatType> d, e;
    std::vector<std::size_t> order;

    void size_for(std::size_t n) {
      if (t.size() < n*n) {
        t.resize(n*n);
        q.resize(n*n);
      }
      if (v.size() < n) {
        v.resize(n);
        p.resize(n);
        sub.resize(n);
        d.resize(n);
        e.resize(n);
        order.resize(n);
      }
    }
  };

  /// Orders indices by the value they point at; a functor for C++98's sake
  template <typename FloatType>
  struct ascending_by {
    FloatType const *values;
    ascending_by(FloatType const *values) : values(values) {}
    bool operator()(std::size_t i, std::size_t j) const {
      return values[i] < values[j];
    }
  };

  /** @brief Eigendecomposition of a small complex Hermitian matrix.

  Written for the dynamical electron diffraction code, where it is the hottest
  operation there is: one per beam group per orientation, some tens of
  thousands per pass over a dataset, on matrices of order ten beams.

  At that size a LAPACK driver is the wrong tool. `zheev` is built for matrices
  where blocking pays and O(n^3) arithmetic dominates; on an 11-beam matrix its
  workspace query, its per-call allocation and its block machinery cost several
  times the arithmetic they are organising, and the allocation becomes
  contention once the reflection loop is threaded. Measured here, `zheev` takes
  about 20 us on an 11-beam matrix where the same algorithm written directly
  for small n takes about 8.

  The algorithm is the standard one and there is nothing novel in it:

    1. Householder reflectors reduce the Hermitian matrix to tridiagonal form,
       accumulating the unitary transform as they go. Each reflector is applied
       as a rank-2 update using the Hermitian symmetry, not as two matrix
       products, so the reduction is O(n^3) overall rather than O(n^4).
    2. The tridiagonal that comes out has a real diagonal but a complex
       subdiagonal. A diagonal unitary makes the subdiagonal real and positive
       without touching the eigenvalues, leaving a real symmetric tridiagonal.
    3. Implicit-shift QL with Wilkinson shifts finds its eigenpairs, the plane
       rotations being accumulated into the same transform.
    4. Eigenpairs are sorted ascending, as every caller here expects.

  On entry @p a is a row-major n by n Hermitian matrix. On exit it holds the
  eigenvectors in its columns -- a[i*n + k] is component i of eigenvector k --
  and @p ev the eigenvalues in ascending order, @p ev having room for n of
  them. That is the layout LAPACK's `heev` leaves behind, so this is a drop-in
  for it.

  Only the upper triangle and the diagonal are read, matching the LAPACK_UPPER
  this replaces, so a caller which has filled only that triangle keeps working.
  The diagonal's imaginary part is ignored rather than trusted, it being zero
  for a Hermitian matrix.

  Accuracy is that of the algorithm, which is backward stable: the computed
  eigenpairs are exact for a matrix within rounding of the one given. It does
  not agree bit for bit with LAPACK -- no two implementations of this do -- but
  the difference is at the level of the last digit or two.
  */
  template <typename FloatType>
  void hermitian_eigen_unblocked(std::complex<FloatType> *a, std::size_t n,
                                 FloatType *ev,
                                 hermitian_eigen_scratch<FloatType> &s)
  {
    typedef std::complex<FloatType> complex_t;
    if (n == 0) {
      return;
    }
    if (n == 1) {
      ev[0] = a[0].real();
      a[0] = complex_t(1, 0);
      return;
    }
    s.size_for(n);

    // The matrix is reduced in place in `t`; `q` accumulates the transform and
    // ends up holding the eigenvectors. Reused storage, so both are filled
    // rather than constructed -- q needs its zeros written every time.
    complex_t *t = &s.t[0], *q = &s.q[0];
    std::copy(a, a + n*n, t);
    std::fill(q, q + n*n, complex_t(0, 0));
    for (std::size_t i = 0; i < n; i++) {
      q[i*n + i] = complex_t(1, 0);
      // mirror the upper triangle down, so the update below may read either
      for (std::size_t j = 0; j < i; j++) {
        t[i*n + j] = std::conj(t[j*n + i]);
      }
      t[i*n + i] = complex_t(t[i*n + i].real(), 0);
    }

    complex_t *v = &s.v[0], *p = &s.p[0], *sub = &s.sub[0];
    std::fill(sub, sub + n, complex_t(0, 0));

    // --- 1. Householder reduction to tridiagonal form ---
    for (std::size_t k = 0; k + 2 < n; k++) {
      // the column below the subdiagonal, which the reflector must annihilate
      FloatType norm_sq = 0;
      for (std::size_t i = k + 1; i < n; i++) {
        norm_sq += std::norm(t[i*n + k]);
      }
      const FloatType norm = std::sqrt(norm_sq);
      if (norm == 0) {
        sub[k] = complex_t(0, 0);
        continue;
      }
      const complex_t x0 = t[(k + 1)*n + k];
      const FloatType abs_x0 = std::abs(x0);
      // reflect onto -phase(x0)*|x|, the sign that avoids cancellation in v
      const complex_t phase = abs_x0 > 0 ? x0/abs_x0 : complex_t(1, 0);
      const complex_t alpha = -phase*norm;

      for (std::size_t i = k + 1; i < n; i++) {
        v[i] = t[i*n + k];
      }
      v[k + 1] -= alpha;
      FloatType v_norm_sq = 0;
      for (std::size_t i = k + 1; i < n; i++) {
        v_norm_sq += std::norm(v[i]);
      }
      if (v_norm_sq == 0) {           // already in the wanted form
        sub[k] = alpha;
        continue;
      }
      const FloatType beta = 2/v_norm_sq;

      /* H A H with H = I - beta v v^H comes out as a rank-2 update,
         A <- A - v w^H - w v^H with w = beta A v - (beta/2)(v^H beta A v) v,
         which is what keeps this O(n^2) per reflector instead of O(n^3). */
      complex_t vhp(0, 0);
      for (std::size_t i = k + 1; i < n; i++) {
        complex_t s(0, 0);
        for (std::size_t j = k + 1; j < n; j++) {
          s += t[i*n + j]*v[j];
        }
        p[i] = beta*s;
        vhp += std::conj(v[i])*p[i];
      }
      // real for a Hermitian matrix; taking the real part keeps it so exactly
      const FloatType kappa = beta*vhp.real()/2;
      for (std::size_t i = k + 1; i < n; i++) {
        p[i] -= kappa*v[i];
      }
      for (std::size_t i = k + 1; i < n; i++) {
        for (std::size_t j = k + 1; j < n; j++) {
          t[i*n + j] -= v[i]*std::conj(p[j]) + p[i]*std::conj(v[j]);
        }
      }
      for (std::size_t i = k + 1; i < n; i++) {
        t[i*n + i] = complex_t(t[i*n + i].real(), 0);
      }
      sub[k] = alpha;

      // q <- q H, again as a rank-1 update rather than a product
      for (std::size_t i = 0; i < n; i++) {
        complex_t s(0, 0);
        for (std::size_t j = k + 1; j < n; j++) {
          s += q[i*n + j]*v[j];
        }
        s *= beta;
        for (std::size_t j = k + 1; j < n; j++) {
          q[i*n + j] -= s*std::conj(v[j]);
        }
      }
    }
    if (n >= 2) {
      sub[n - 2] = t[(n - 1)*n + (n - 2)];
    }

    // --- 2. make the subdiagonal real, by a diagonal unitary ---
    /* Multiplying column k of q by a phase changes no eigenvalue and no
       eigenvector up to that same phase, which nothing downstream can see. */
    /* e is one longer than the subdiagonal it holds, the last entry staying
       zero: the search for a negligible subdiagonal below runs off the end
       when the block is already split, and finds the zero there. */
    FloatType *d = &s.d[0], *e = &s.e[0];
    std::fill(e, e + n, FloatType(0));
    for (std::size_t i = 0; i < n; i++) {
      d[i] = t[i*n + i].real();
    }
    /* With D(k+1) = D(k)*phase(sub[k]), element (k+1,k) of D^H T D comes out
       as conj(phase)*sub[k] = |sub[k]|, and the diagonal is untouched because
       |D(k)| is one. */
    complex_t accumulated(1, 0);
    for (std::size_t k = 0; k + 1 < n; k++) {
      const FloatType abs_sub = std::abs(sub[k]);
      e[k] = abs_sub;
      if (abs_sub > 0) {
        accumulated *= sub[k]/abs_sub;
      }
      for (std::size_t i = 0; i < n; i++) {
        q[i*n + (k + 1)] *= accumulated;
      }
    }

    // --- 3. implicit-shift QL on the real symmetric tridiagonal ---
    const FloatType eps = std::numeric_limits<FloatType>::epsilon();
    for (std::size_t l = 0; l < n; l++) {
      int iter = 0;
      std::size_t m;
      do {
        for (m = l; m + 1 < n; m++) {
          const FloatType dd = std::abs(d[m]) + std::abs(d[m + 1]);
          if (std::abs(e[m]) <= eps*dd) {
            break;
          }
        }
        if (m == l) {
          break;
        }
        SMTBX_ASSERT(iter++ < 60);   // no tridiagonal of this size should
        // Wilkinson shift, chosen as the eigenvalue of the trailing 2x2
        // nearer to d[l], which is what makes the convergence cubic
        FloatType g = (d[l + 1] - d[l])/(2*e[l]);
        FloatType r = std::sqrt(g*g + 1);
        g = d[m] - d[l] + e[l]/(g + (g >= 0 ? std::abs(r) : -std::abs(r)));
        FloatType s = 1, c = 1, pp = 0;
        std::size_t i = m;
        while (i-- > l) {
          FloatType f = s*e[i], b = c*e[i];
          r = std::sqrt(f*f + g*g);
          e[i + 1] = r;
          if (r == 0) {               // a zero subdiagonal: deflate and restart
            d[i + 1] -= pp;
            e[m] = 0;
            break;
          }
          s = f/r;
          c = g/r;
          g = d[i + 1] - pp;
          r = (d[i] - g)*s + 2*c*b;
          pp = s*r;
          d[i + 1] = g + pp;
          g = c*r - b;
          // the same rotation, applied to the two eigenvector columns
          complex_t *qi = &q[i], *qi1 = &q[i + 1];
          for (std::size_t row = 0; row < n; row++, qi += n, qi1 += n) {
            const complex_t u = *qi, w = *qi1;
            *qi1 = s*u + c*w;
            *qi = c*u - s*w;
          }
        }
        if (r == 0 && i + 1 > l) {
          continue;
        }
        d[l] -= pp;
        e[l] = g;
        e[m] = 0;
      } while (m != l);
    }

    // --- 4. ascending, which is what heev would have given ---
    std::size_t *order = &s.order[0];
    for (std::size_t i = 0; i < n; i++) {
      order[i] = i;
    }
    std::sort(order, order + n, ascending_by<FloatType>(d));
    for (std::size_t k = 0; k < n; k++) {
      ev[k] = d[order[k]];   // NOLINT: a is rewritten from q, not read here
      const std::size_t src = order[k];
      for (std::size_t i = 0; i < n; i++) {
        a[i*n + k] = q[i*n + src];
      }
    }
  }

  /** @brief Where LAPACK starts winning again.

  The advantage above is an advantage at small n only, and it has to end: the
  reduction here is unblocked, so once the working set stops fitting in cache
  LAPACK's blocking is worth more than the per-call overhead it costs. Measured
  on random Hermitian matrices, us per call, this against LAPACK zheev:

      n     200     300     400     600     800
      ratio 1.78x   1.58x   1.27x   0.92x   0.74x

  so the crossover sits between 400 and 600, and LAPACK's lead keeps growing
  above it while the lead below is capped at about 3x. The threshold therefore
  sits at the crossover rather than beyond it, the curve being shallow enough
  either side that the exact placement costs little.

  Well outside anything a beam group reaches -- `beam_n` defaults to 10 and a
  few tens is a large calculation -- so this is about not being surprising if
  someone turns it up, not about a case anyone runs today.
  */
  inline std::size_t hermitian_eigen_blocked_from() { return 512; }

  /** @brief Eigendecomposition of a Hermitian matrix, by whichever is faster.

  See hermitian_eigen_unblocked for the layout, which is LAPACK's, and for what
  the two implementations do and do not share. They agree to the last digit or
  two but not bit for bit, so a matrix either side of the threshold is solved
  to the same accuracy by a different route.

  Without fast_linalg there is no LAPACK to defer to and the unblocked path is
  used at every size. That is slower for a very large matrix and correct at
  all of them, which is the right way round: such a build previously could not
  run this code at all, every call throwing.

  Compiled in is not the same as available: fast_linalg loads its BLAS at run
  time and every entry point asserts that this has happened. So the test is on
  is_initialised() and not on the macro alone -- otherwise a process which
  never called it would throw here, having a working solver right beside it.
  */
  template <typename FloatType>
  void hermitian_eigen(std::complex<FloatType> *a, std::size_t n,
                       FloatType *ev,
                       hermitian_eigen_scratch<FloatType> &s)
  {
#if defined(USE_FAST_LINALG)
    if (n >= hermitian_eigen_blocked_from() && fast_linalg::is_initialised()) {
      lapack_int info = fast_linalg::heev(
        fast_linalg::LAPACK_ROW_MAJOR, 'V', fast_linalg::LAPACK_UPPER,
        static_cast<lapack_int>(n), a, static_cast<lapack_int>(n), ev);
      SMTBX_ASSERT(!info)(info);
      return;
    }
#endif
    hermitian_eigen_unblocked(a, n, ev, s);
  }

  /// For a caller with nowhere to keep scratch, at the cost of allocating it
  template <typename FloatType>
  void hermitian_eigen(std::complex<FloatType> *a, std::size_t n,
                       FloatType *ev)
  {
    hermitian_eigen_scratch<FloatType> s;
    hermitian_eigen(a, n, ev, s);
  }

}}
