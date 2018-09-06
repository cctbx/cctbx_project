#pragma once

#if defined(CCTBX_FAST_LINALG_USES_OPENBLAS)
#define LAPACK_COMPLEX_CUSTOM
#include <complex>
#define lapack_complex_float  std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

#else
#error No LAPACKE has been configured
#endif

/// Overloads of LAPACKE functions
/**
 *   Some terms used throughout:
 *
 *   - RFP, for Rectangular Full Packed: a format for symmetric and triangular
 *     matrices that takes the same amount of space as the traditional
 *     packed format but allows BLAS level 3 algorithms to be used.
 */
namespace fast_linalg {
  /// @name Cholesky factorisation of A in-place (RFP)
  //@{
  inline lapack_int
  pftrf( int matrix_order, char transr, char uplo,
        lapack_int n, double* a )
  {
    return LAPACKE_dpftrf(matrix_order, transr, uplo, n, a);
  }

  inline lapack_int
  pftrf( int matrix_order, char transr, char uplo,
        lapack_int n, float* a )
  {
    return LAPACKE_spftrf(matrix_order, transr, uplo, n, a);
  }
  //@}

  /// @name Rank-k update (RFP)
  //@{
  inline lapack_int
  sfrk(int matrix_order, char transr, char uplo, char trans,
       lapack_int n, lapack_int k, double alpha,
       const double* a, lapack_int lda, double beta,
       double* c)
  {
    return LAPACKE_dsfrk(matrix_order, transr, uplo, trans,
                         n, k, alpha,
                         a, lda, beta,
                         c);
  }

  inline lapack_int
  sfrk(int matrix_order, char transr, char uplo, char trans,
       lapack_int n, lapack_int k, float alpha,
       const float* a, lapack_int lda, float beta,
       float* c)
  {
    return LAPACKE_ssfrk(matrix_order, transr, uplo, trans,
                         n, k, alpha,
                         a, lda, beta,
                         c);
  }
  //@}

  /// @name Copies arf in RFP format to ap in standard packed format
  //@{
  inline lapack_int
  tfttp(int matrix_order, char transr, char uplo,
        lapack_int n, const double* arf, double* ap )
  {
    return LAPACKE_dtfttp(matrix_order, transr, uplo, n, arf, ap);
  }

  inline lapack_int
  tfttp(int matrix_order, char transr, char uplo,
        lapack_int n, const float* arf, float* ap )
  {
    return LAPACKE_stfttp(matrix_order, transr, uplo, n, arf, ap);
  }
  //@}

  /// @name Copies ap in standard packed format to arf in RFP format
  //@{
  inline lapack_int
  tpttf(int matrix_order, char transr, char uplo,
        lapack_int n, const double* ap, double* arf)
  {
    return LAPACKE_dtpttf(matrix_order, transr, uplo, n, ap, arf);
  }

  inline lapack_int
  tpttf(int matrix_order, char transr, char uplo,
        lapack_int n, const float* ap, float* arf)
  {
    return LAPACKE_stpttf(matrix_order, transr, uplo, n, ap, arf);
  }
  //@}

  /// @name Inverse of Cholesky decomposed matrix (RFP)
  //@{
  inline lapack_int
  pftri(int matrix_order, char transr, char uplo,
        lapack_int n, double* a )
  {
    return LAPACKE_dpftri(matrix_order, transr, uplo, n, a);
  }

  inline lapack_int
  pftri(int matrix_order, char transr, char uplo,
        lapack_int n, float* a )
  {
    return LAPACKE_spftri(matrix_order, transr, uplo, n, a);
  }
  //@}

}
