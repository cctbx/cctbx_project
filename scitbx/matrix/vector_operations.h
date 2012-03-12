#ifndef SCITBX_MATRIX_VECTOR_OPERATIONS_H
#define SCITBX_MATRIX_VECTOR_OPERATIONS_H

namespace scitbx { namespace matrix {

  /// x(0:n) *= alpha
  template <typename T>
  inline
  void scale_vector(int n, T *x, T alpha) {
    if (alpha == 0) std::fill_n(x, n, T(0));
    else if (alpha != 1) for (int i=0; i<n; ++i) x[i] *= alpha;
  }

  /// Euclidean scalar product of x(0:n) and y(0:n)
  template <typename T>
  inline
  T dot(int n, T const *x, T const *y) {
    T s = 0;
    for (int i=0; i<n; ++i) s += x[i] * y[i];
    return s;
  }

}}

#endif // GUARD
