#ifndef SCITBX_MATRIX_MATRIX_VECTOR_OPERATIONS_H
#define SCITBX_MATRIX_MATRIX_VECTOR_OPERATIONS_H

#include <scitbx/matrix/vector_operations.h>

namespace scitbx { namespace matrix {

  /// y := alpha A x + beta y for a general m x n matrix A
  /** A shall be passed in row-major layout to a */
  template <typename T>
  void matrix_vector(int m, int n,
                     T const *a, T const *x, T *y,
                     T alpha=1, T beta=0)
  {
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
  /** A shall be passed in row-major layout to a */
  template <typename T>
  void matrix_transposed_vector(int m, int n,
                                T const *a, T const *x, T *y,
                                T alpha=1, T beta=0)
  {
    scale_vector(m, y, beta);
    if (alpha == 0) return;
    for (int i=0; i<m; ++i) {
      T alpha_x_i = alpha * x[i];
      for (int j=0; j<n; ++j) {
        T a_ij = *a++;
        y[j] += a_ij * alpha_x_i;
      }
    }
  }

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
   */
  template <typename T>
  void symmetric_packed_u_rank_1_update(int n,
                                        T *a, T const *x,
                                        T alpha=1)
  {
    for (int i=0; i<n; ++i) {
      T alpha_x_i = alpha * x[i];
      for (int j=i; j<n; ++j) {
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
