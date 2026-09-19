#ifndef SCITBX_MATRIX_MATRIX_VECTOR_OPERATIONS_H
#define SCITBX_MATRIX_MATRIX_VECTOR_OPERATIONS_H

#include <scitbx/matrix/vector_operations.h>
#include <fast_linalg/cblas.h>
#include <fast_linalg/lapacke.h>
#include <vector>
#include <cstddef>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace scitbx { namespace matrix {

  /// Where a BLAS call starts being worth making, in elements of the matrix
  /** A gemv call costs a function pointer indirection, the library's own
      dispatch and its decision about threading -- a microsecond or so, against
      the twenty nanoseconds a 6x6 product takes in the loop below. Nearly
      every caller in cctbx passes 6x6 (ADP restraints, geometry, Hirshfeld
      test) and does so once per restraint, so sending those through a library
      would be slower and would change their arithmetic for nothing. The note
      on symmetric_packed_u_rank_1_update below records the same finding for
      spr, arrived at the same way.

      A design matrix is the case this exists for: at that size the loop is not
      remotely competitive.
   */
  static const long blas_worthwhile_size = 10000;

  /// y := alpha A x + beta y for a general m x n matrix A
  /** A shall be passed in row-major layout to a.
      x is n long and y is m long.
   */
  template <typename T>
  void matrix_vector(int m, int n,
                     T const *a, T const *x, T *y,
                     T alpha=1, T beta=0)
  {
    if (static_cast<long>(m)*n >= blas_worthwhile_size
        && fast_linalg::is_initialised()) {
      using namespace fast_linalg;
      gemv(CblasRowMajor, CblasNoTrans, m, n, alpha, a, n, x, 1, beta, y, 1);
      return;
    }
    scale_vector(m, y, beta);
    if (alpha == 0) return;
    for (int i=0; i<m; ++i) {
      T ax_i = 0;
      for (int j=0; j<n; ++j) {
        T a_ij = *a++;
        ax_i += a_ij * x[j];
      }
      y[i] += alpha * ax_i;
    }
  }

  /// y := alpha A^T x + beta y for a general m x n matrix A
  /** A shall be passed in row-major layout to a.
      x is m long and y is n long -- the other way round from matrix_vector,
      m and n always describing A itself rather than A^T.
   */
  template <typename T>
  void matrix_transposed_vector(int m, int n,
                                T const *a, T const *x, T *y,
                                T alpha=1, T beta=0)
  {
    if (static_cast<long>(m)*n >= blas_worthwhile_size
        && fast_linalg::is_initialised()) {
      using namespace fast_linalg;
      gemv(CblasRowMajor, CblasTrans, m, n, alpha, a, n, x, 1, beta, y, 1);
      return;
    }
    // n and not m: y is as long as A has columns. Every caller of this was
    // square until the design matrix of a refinement came along, which hid
    // that scaling m of them either left the tail of y unscaled or ran off
    // the end of it.
    scale_vector(n, y, beta);
    if (alpha == 0) return;
    for (int i=0; i<m; ++i) {
      T alpha_x_i = alpha * x[i];
      for (int j=0; j<n; ++j) {
        T a_ij = *a++;
        y[j] += a_ij * alpha_x_i;
      }
    }
  }

  /// @name Mixed precision products: a matrix stored narrow, accumulated wide
  /** A stored design matrix is read once per product and multiplied by once, so
      the products are bandwidth bound and not flop bound. Holding the matrix in
      float and accumulating in double therefore halves the time for a
      negligible loss of accuracy -- far better than accumulating in single as a
      plain sgemv would. See notes/gemv_precision.cpp.

      **No BLAS offers this**, which is why it is written out here: sgemv
      accumulates in single, and expressing the mix through numpy upcasts the
      matrix, allocating the memory the exercise exists to avoid. The loops are
      trivial and memory bound, so a library would do no better.

      Unlike the two above these take no alpha or beta: every caller wants
      y := A x exactly, and the general form would only invite the aliasing
      question these avoid by writing y once.
   */
  //@{
  template <typename StoreType, typename AccType>
  void matrix_vector_mixed(int m, int n,
                           StoreType const *a, AccType const *x, AccType *y)
  {
    #if defined(_OPENMP)
    #pragma omp parallel for schedule(static)
    #endif
    for (long long i = 0; i < static_cast<long long>(m); ++i) {
      StoreType const *row = a + static_cast<std::size_t>(i)*n;
      AccType s = 0;
      for (int j = 0; j < n; ++j) s += static_cast<AccType>(row[j])*x[j];
      y[i] = s;
    }
  }

  template <typename StoreType, typename AccType>
  void matrix_transposed_vector_mixed(int m, int n,
                                      StoreType const *a, AccType const *x,
                                      AccType *y)
  {
    /* Each thread accumulates a private copy of y and they are summed at the
       end. y is as long as A has columns, small enough that the private copies
       stay in cache and A is still streamed exactly once. Accumulating
       into a shared y instead would race, and locking per row would serialise
       the whole thing.
     */
    int n_threads = 1;
    #if defined(_OPENMP)
    n_threads = omp_get_max_threads();
    #endif
    std::vector<AccType> partial(
      static_cast<std::size_t>(n_threads)*n, AccType(0));
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
      int t = 0;
      #if defined(_OPENMP)
      t = omp_get_thread_num();
      #endif
      AccType *p = &partial[static_cast<std::size_t>(t)*n];
      #if defined(_OPENMP)
      #pragma omp for schedule(static)
      #endif
      for (long long i = 0; i < static_cast<long long>(m); ++i) {
        StoreType const *row = a + static_cast<std::size_t>(i)*n;
        AccType const x_i = x[i];
        for (int j = 0; j < n; ++j) p[j] += static_cast<AccType>(row[j])*x_i;
      }
    }
    for (int j = 0; j < n; ++j) {
      AccType s = 0;
      for (int t = 0; t < n_threads; ++t) {
        s += partial[static_cast<std::size_t>(t)*n + j];
      }
      y[j] = s;
    }
  }
  //@}

  /// y := alpha A x + beta y for symmetric matrix A
  /** The upper diagonal of A packed by row is passed in the range
      [ a, a + n*(n+1)/2 )
   */
  template <typename T>
  void symmetric_packed_u_vector(int n,
                                 T const *a, T const *x, T *y,
                                 T alpha=1, T beta=0)
  {
    scale_vector(n, y, beta);
    if (alpha == 0) return;
    for (unsigned i=0; i<n; ++i) {
      T alpha_x_i = alpha * x[i];
      T ax_i = 0;
      T a_ii = *a++;
      y[i] += a_ii * alpha_x_i;
      for (unsigned j=i+1; j<n; ++j) {
        T a_ij = *a++, &a_ji = a_ij;
        ax_i += a_ij * x[j];
        y[j] += a_ji * alpha_x_i;
      }
      y[i] += alpha * ax_i;
    }
  }


  /// A := alpha x x^T + A for symmetric matrix A
  /** The upper diagonal of A packed by rows shall be passed in the range
      [ a, a +  n*(n+1)/2 )
    Note that using fast_linalg is slower here
    if (fast_linalg::is_initialised()) {
      using namespace fast_linalg;
      spr(CblasRowMajor, CblasUpper, n, alpha, x, 1, a);
    }
  */
  template <typename T>
  void symmetric_packed_u_rank_1_update(int n,
                                        T *a, T const *x,
                                        T alpha=1)
  {
    if (alpha == 0.0) {
      return;
    }
    for (int i = 0; i < n; ++i) {
      int len = (n - i);
      if (x[i] == 0) {
        a += len;
        continue;
      }
      T alpha_x_i = alpha * x[i];
      for (int j = i; j < n; ++j) {
        *a++ += alpha_x_i * x[j];
      }
    }
  }

  /// x^T A x for symmetric matrix A
  /** The upper diagonal of A packed by rows shall be passed in the range
   [ a, a +  n*(n+1)/2 )
   */
  template <typename T>
  T quadratic_form_packed_u(int n, T const *a, T const *x) {
    T diag = 0, off_diag = 0;
    for (int i=0; i<n; ++i) {
      T a_ii = *a++;
      diag += x[i] * a_ii * x[i];
      T off_diag_i = 0;
      for (int j=i+1; j<n; ++j) {
        T a_ij = *a++;
        off_diag_i += a_ij * x[j];
      }
      off_diag += x[i] * off_diag_i;
    }
    return diag + 2.*off_diag;
  }

  /// x^T A x for a square matrix A
  /** A shall be passed as a in row-major layout */
  template <typename T>
  T quadratic_form(int n, T const *a, T const *x) {
    T s = 0;
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) {
      T a_ij = *a++;
      s += x[i] * a_ij * x[j];
    }
    return s;
  }

  /// x := P x where P is a permutation matrix
  /** P = E_{r-1} ... E_0
      where each E_k is the identity with row k and p[k] interchanged.
   */
  template <typename T, typename IndexType>
  void permutation_vector(int n, T *x, IndexType const *p) {
    for (int i=0; i<n; ++i) {
      if (p[i] != i) std::swap(x[i], x[ p[i] ]);
    }
  }

  /// x := P^T x where P is a permutation matrix
  /** P = E_{r-1} ... E_0
   where each E_k is the identity with row k and p[k] interchanged.
   Note: P^T = P^{-1}
   */
  template <typename T, typename IndexType>
  void permutation_transposed_vector(int n, T *x, IndexType const *p) {
    for (int i=n-1; i>=0; --i) {
      if (p[i] != i) std::swap(x[i], x[ p[i] ]);
    }
  }
}}

#endif // GUARD
