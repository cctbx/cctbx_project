#ifndef SCITBX_MATRIX_DIAGONAL_H
#define SCITBX_MATRIX_DIAGONAL_H

namespace scitbx { namespace matrix {

  template <typename NumType>
  inline
  void
  diagonal(
    const NumType *a,
    std::size_t n,
    NumType *d)
  {
    std::size_t ii = 0;
    for (std::size_t i=0;i<n;i++,ii+=n+1) {
      d[i] = a[ii];
    }
  }

  template <typename NumType>
  inline
  NumType
  diagonal_sum(
    const NumType *a,
    std::size_t n)
  {
    NumType result = 0;
    std::size_t ii = 0;
    for (std::size_t i=0;i<n;i++,ii+=n+1) {
      result += a[ii];
    }
    return result;
  }

  template <typename NumType>
  inline
  NumType
  diagonal_product(
    const NumType *a,
    std::size_t n)
  {
    NumType result = 1;
    std::size_t ii = 0;
    for (std::size_t i=0;i<n;i++,ii+=n+1) {
      result *= a[ii];
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_DIAGONAL_H
