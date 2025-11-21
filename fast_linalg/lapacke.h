#pragma once

/// Overloads of LAPACKE functions
/**
 *   Some terms used throughout:
 *
 *   - RFP, for Rectangular Full Packed: a format for symmetric and triangular
 *     matrices that takes the same amount of space as the traditional
 *     packed format but allows BLAS level 3 algorithms to be used.
 */
#include <string>
#include <scitbx/error.h>
#include <complex>
#ifndef lapack_int
 // note: this can be long
#define lapack_int int
#endif

namespace fast_linalg {
  const int LAPACK_ROW_MAJOR = 101,
    LAPACK_COL_MAJOR = 102;

  const char LAPACK_EIGENVALUES = 'N',
    LAPACK_EIGENVALUES_AND_EIGENVECTORS = 'V',
    LAPACK_UPPER = 'U',
    LAPACK_LOWER = 'L';
}

#if defined(USE_FAST_LINALG)
#include <boost/config.hpp>
#ifndef fast_linalg_api
#define fast_linalg_api BOOST_SYMBOL_EXPORT
#endif

extern "C" {
  fast_linalg_api bool is_fast_linalg_initialised();
  // throws exceptions
  fast_linalg_api void initialise_fast_linalg(const std::string &lib_name);

  /// deletes stored wrapper pointer and the shared library handle
  fast_linalg_api void finalise_fast_linalg();

  fast_linalg_api lapack_int lapack_pftrf(int matrix_order, char transr,
    char uplo, lapack_int n, double* a);
  fast_linalg_api lapack_int lapack_spftrf(int matrix_order, char transr,
    char uplo, lapack_int n, float* a);

  fast_linalg_api lapack_int lapack_sfrk(int matrix_order, char transr,
    char uplo, char trans, lapack_int n, lapack_int k, double alpha,
    const double* a, lapack_int lda, double beta,
    double* c);
  fast_linalg_api lapack_int lapack_ssfrk(int matrix_order, char transr,
    char uplo, char trans, lapack_int n, lapack_int k, float alpha,
    const float* a, lapack_int lda, float beta,
    float* c);

  fast_linalg_api lapack_int lapack_tfttp(int matrix_order, char transr,
    char uplo, lapack_int n, const double* arf, double* ap);
  fast_linalg_api lapack_int lapack_stfttp(int matrix_order, char transr,
    char uplo, lapack_int n, const float* arf, float* ap);

  fast_linalg_api lapack_int lapack_tpttf(int matrix_order, char transr,
    char uplo, lapack_int n, const double* ap, double* arf);
  fast_linalg_api lapack_int lapack_stpttf(int matrix_order, char transr,
    char uplo, lapack_int n, const float* ap, float* arf);

  fast_linalg_api lapack_int lapack_pftri(int matrix_order, char transr,
    char uplo, lapack_int n, double* a);
  fast_linalg_api lapack_int lapack_spftri(int matrix_order, char transr,
    char uplo, lapack_int n, float* a);

  fast_linalg_api lapack_int lapack_dsyev(int matrix_order, char jobz,
    char uplo, lapack_int n, double* a, lapack_int lda,
    double* w);
  fast_linalg_api lapack_int lapack_ssyev(int matrix_order, char jobz,
    char uplo, lapack_int n, float* a, lapack_int lda,
    float* w);

  fast_linalg_api lapack_int lapack_zheev(int matrix_order, char jobz,
    char uplo, lapack_int n, std::complex<double>* a, lapack_int lda,
    double* w);
  fast_linalg_api lapack_int lapack_cheev(int matrix_order, char jobz,
    char uplo, lapack_int n, std::complex<float>* a, lapack_int lda,
    float* w);

  fast_linalg_api lapack_int lapack_zgeev(int matrix_order, char jobvl,
    char jovr, lapack_int n, std::complex<double>* a, lapack_int lda,
    std::complex<double>* w, std::complex<double>* vl,
    lapack_int ldvl, std::complex<double>* vr, lapack_int ldvr);
  fast_linalg_api lapack_int lapack_cgeev(int matrix_order, char jobvl,
    char jovr,lapack_int n, std::complex<float>* a, lapack_int lda,
    std::complex<float>* w, std::complex<float>*vl, lapack_int ldvl,
    std::complex<float>*vr, lapack_int ldvr);

  fast_linalg_api lapack_int lapack_cgetrf(int matrix_order,
    lapack_int m, lapack_int n, std::complex<float>* a, lapack_int lda,
    lapack_int* ipiv);
  fast_linalg_api lapack_int lapack_zgetrf(int matrix_order,
    lapack_int m, lapack_int n, std::complex<double>* a, lapack_int lda,
    lapack_int* ipiv);

  fast_linalg_api lapack_int lapack_cgetri(int matrix_order,
    lapack_int n, std::complex<float>* a, lapack_int lda,
    lapack_int* ipiv);
  fast_linalg_api lapack_int lapack_zgetri(int matrix_order,
    lapack_int n, std::complex<double>* a, lapack_int lda,
    lapack_int* ipiv);

  fast_linalg_api void cblas_ssyr(int Order, int Uplo, int N, float Alpha,
    const float *X, int incX, float* A, int lda);
  fast_linalg_api void cblas_dsyr(int Order, int Uplo, int N, double Alpha,
    const double *X, int incX, double* A, int lda);

  fast_linalg_api void cblas_ssyr(int Order, int Uplo, int N, float alpha,
    const float *X, int incX, float *A, int lda);
  fast_linalg_api void cblas_dsyr(int Order, int Uplo, int N, double alpha,
    const double *X, int incX, double *A, int lda);

  fast_linalg_api void cblas_sspr(int Order, int Uplo, int N, float alpha,
    const float* X, int incX, float* A);
  fast_linalg_api void cblas_dspr(int Order, int Uplo, int N, double alpha,
    const double* X, int incX, double* A);

  fast_linalg_api void cblas_ssyrk(int Order, int Uplo, int Trans, int N, int K,
    float alpha, const float *A, int lda, float beta, float *C, int ldc);
  fast_linalg_api void cblas_dsyrk(int Order, int Uplo, int Trans, int N, int K,
    double alpha, const double* A, int lda, double beta, double* C, int ldc);

  fast_linalg_api void cblas_sgemm(int Order, int TransA, int TransB,
    int M, int N, int K,
    float alpha, const float* A, int lda,
    const float* B, int ldb, float beta, float* C, int ldc);
  fast_linalg_api void cblas_dgemm(int Order, int TransA, int TransB,
    int M, int N, int K,
    double alpha, const double* A, int lda,
    const double* B, int ldb, double beta, double* C, int ldc);

};

namespace fast_linalg {
  inline bool is_initialised() {
    return is_fast_linalg_initialised();
  }

  /// pass lib name like 'openblas.so' or full library path
  inline void initialise(const std::string &lib_name) {
    initialise_fast_linalg(lib_name);
  }

  /// cleans up the library
  inline void finalise() {
    finalise_fast_linalg();
  }

  /// @name Cholesky factorisation of A in-place (RFP)
  //@{
  inline lapack_int pftrf(int matrix_order, char transr, char uplo,
    lapack_int n, double* a)
  {
    return lapack_pftrf(matrix_order, transr, uplo, n, a);
  }

  inline lapack_int pftrf(int matrix_order, char transr, char uplo,
    lapack_int n, float* a)
  {
    return lapack_spftrf(matrix_order, transr, uplo, n, a);
  }
  //@}

  /// @name Rank-k update (RFP)
  //@{
  inline lapack_int sfrk(int matrix_order, char transr,
    char uplo, char trans, lapack_int n, lapack_int k, double alpha,
    const double* a, lapack_int lda, double beta,
    double* c)
  {
    return lapack_sfrk(matrix_order, transr,
      uplo, trans, n, k, alpha, a, lda, beta, c);
  }

  inline lapack_int sfrk(int matrix_order, char transr,
    char uplo, char trans, lapack_int n, lapack_int k, float alpha,
    const float* a, lapack_int lda, float beta,
    float* c)
  {
    return lapack_ssfrk(matrix_order, transr,
      uplo, trans, n, k, alpha, a, lda, beta, c);
  }
  //@}

  /// @name Copies arf in RFP format to ap in standard packed format
  //@{
  inline lapack_int tfttp(int matrix_order, char transr,
    char uplo, lapack_int n, const double* arf, double* ap)
  {
    return lapack_tfttp(matrix_order, transr, uplo, n, arf, ap);
  }

  inline lapack_int tfttp(int matrix_order, char transr,
    char uplo, lapack_int n, const float* arf, float* ap)
  {
    return lapack_stfttp(matrix_order, transr, uplo, n, arf, ap);
  }
  //@}

  /// @name Copies ap in standard packed format to arf in RFP format
  //@{
  inline lapack_int tpttf(int matrix_order, char transr, char uplo,
    lapack_int n, const double* ap, double* arf)
  {
    return lapack_tpttf(matrix_order, transr, uplo, n, ap, arf);
  }

  inline lapack_int tpttf(int matrix_order, char transr, char uplo,
    lapack_int n, const float* ap, float* arf)
  {
    return lapack_stpttf(matrix_order, transr, uplo, n, ap, arf);
  }
  //@}

  /// @name Inverse of Cholesky decomposed matrix (RFP)
  //@{
  inline lapack_int pftri(int matrix_order, char transr, char uplo,
    lapack_int n, double* a)
  {
    return lapack_pftri(matrix_order, transr, uplo, n, a);
  }

  inline lapack_int pftri(int matrix_order, char transr, char uplo,
    lapack_int n, float* a)
  {
    return lapack_spftri(matrix_order, transr, uplo, n, a);
  }
  //@}

  /// @name Performs eigenvalue decomposition of a real-symm matrix
  //@{

  inline lapack_int syev(int matrix_order, char jobz,
    char uplo, lapack_int n, float* a, lapack_int lda,
    float* w)
  {
    return lapack_ssyev(matrix_order, jobz, uplo, n, a, lda, w);
  }

  inline lapack_int syev(int matrix_order, char jobz,
    char uplo, lapack_int n, double* a, lapack_int lda,
    double* w)
  {
    return lapack_dsyev(matrix_order, jobz, uplo, n, a, lda, w);
  }
  //@}

  /// @name Performs eigenvalue decomposition of a Hermitian matrix
  //@{

  inline lapack_int heev(int matrix_order, char jobz,
    char uplo, lapack_int n, std::complex<float>* a, lapack_int lda,
    float* w)
  {
    return lapack_cheev(matrix_order, jobz, uplo, n, a, lda, w);
  }

  inline lapack_int heev(int matrix_order, char jobz,
    char uplo, lapack_int n, std::complex<double>* a, lapack_int lda,
    double* w)
  {
    return lapack_zheev(matrix_order, jobz, uplo, n, a, lda, w);
  }
  //@}

  /// @name Performs eigenvalue decomposition of a generic complex matrix
  //@{

  inline lapack_int geev(int matrix_order, char jobvl, char jobvr,
    lapack_int n, std::complex<float>* a, lapack_int lda,
    std::complex<float>* w, std::complex<float>* vl, lapack_int ldvl,
    std::complex<float>* vr, lapack_int ldvr)
  {
    return lapack_cgeev(matrix_order, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
  }

  inline lapack_int geev(int matrix_order, char jobvl, char jobvr,
    lapack_int n, std::complex<double>* a, lapack_int lda,
    std::complex<double>* w, std::complex<double>* vl, lapack_int ldvl,
    std::complex<double>* vr, lapack_int ldvr)
  {
    return lapack_zgeev(matrix_order, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
  }
  //@}

  /// @name Performs LU decomosition of a complex matrix
  //@{
  inline lapack_int getrf(int matrix_order,
    lapack_int m, lapack_int n, std::complex<float>* a, lapack_int lda,
    lapack_int* ipiv)
  {
    return lapack_cgetrf(matrix_order, m, n, a, lda, ipiv);
  }
  inline lapack_int getrf(int matrix_order,
    lapack_int m, lapack_int n, std::complex<double>* a, lapack_int lda,
    lapack_int* ipiv)
  {
    return lapack_zgetrf(matrix_order, m, n, a, lda, ipiv);
  }
  //@}

  /// @name Performs complex matrix inversion after a call to getrf
  //@{
  inline lapack_int lapack_getri(int matrix_order,
    lapack_int n, std::complex<float>* a, lapack_int lda,
    lapack_int* ipiv)
  {
    return lapack_cgetri(matrix_order, n, a, lda, ipiv);
  }

  inline lapack_int getri(int matrix_order,
    lapack_int n, std::complex<double>* a, lapack_int lda,
    lapack_int* ipiv)
  {
    return lapack_zgetri(matrix_order, n, a, lda, ipiv);
  }
  //@}

  /// @name Performs a symmetric rank-1 update (upacked matrix!)
  //@{
  inline void syr(int Order, int Uplo, int N, float Alpha,
    const float *X, int incX, float* A, int lda)
  {
    cblas_ssyr(Order, Uplo, N, Alpha, X, incX, A, lda);
  }

  inline void syr(int Order, int Uplo, int N, double Alpha,
    const double *X, int incX, double* A, int lda)
  {
    cblas_dsyr(Order, Uplo, N, Alpha, X, incX, A, lda);
  }
  //@}

  /// @name Performs a symmetric rank-1 update (packed matrix)
  //@{
  inline void spr(int Order, int Uplo, int N, float Alpha,
    const float* X, int incX, float* A)
  {
    cblas_sspr(Order, Uplo, N, Alpha, X, incX, A);
  }

  inline void spr(int Order, int Uplo, int N, double Alpha,
    const double* X, int incX, double* A)
  {
    cblas_dspr(Order, Uplo, N, Alpha, X, incX, A);
  }
  //@}

  /// @name Performs a symmetric rank-k update
  //@{
  inline void syrk(int Order, int Uplo, int Trans, int N, int K,
    float alpha, const float *A, int lda, float beta, float *C, int ldc)
  {
    cblas_ssyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
  }

  inline void syrk(int Order, int Uplo, int Trans, int N, int K,
    double alpha, const double* A, int lda, double beta, double* C, int ldc)
  {
    cblas_dsyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
  }
  //@}

  /// @name Performs a generic matrix multiplication
  //@{
  inline void gemm(int Order, int TransA, int TransB,
    int M, int N, int K,
    float alpha, const float* A, int lda,
    const float* B, int ldb, float beta, float* C, int ldc)
  {
    cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  inline void gemm(int Order, int TransA, int TransB,
    int M, int N, int K,
    double alpha, const double* A, int lda,
    const double* B, int ldb, double beta, double* C, int ldc)
  {
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
  //@}

}
#else
namespace fast_linalg {

  inline bool is_initialised() {
    return false;
  }

  inline void initialise(const std::string &lib_name) {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  inline void finalise() {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int pftrf(int matrix_order, char transr, char uplo,
    lapack_int n, FloatType* a)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int sfrk(int matrix_order, char transr,
    char uplo, char trans, lapack_int n, lapack_int k, FloatType alpha,
    const FloatType* a, lapack_int lda, FloatType beta,
    FloatType* c)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int tfttp(int matrix_order, char transr,
    char uplo, lapack_int n, const FloatType* arf, FloatType* ap)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int tpttf(int matrix_order, char transr, char uplo,
    lapack_int n, const FloatType* ap, FloatType* arf)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int pftri(int matrix_order, char transr, char uplo,
    lapack_int n, FloatType* a)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  void spr(int, int, int, FloatType, const FloatType *, int,
    FloatType *)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  inline lapack_int syev(int, char, char, lapack_int, FloatType*,
    lapack_int, FloatType*)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  inline lapack_int heev(int, char, char, lapack_int, std::complex<FloatType>*,
    lapack_int, FloatType*)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  inline lapack_int geev(int, char, char, lapack_int,
    std::complex<FloatType>*, lapack_int,
    std::complex<FloatType>*, std::complex<FloatType>*,
    lapack_int, std::complex<FloatType>*, lapack_int)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int getrf(int, lapack_int, lapack_int,
    std::complex<FloatType>* , lapack_int,
    std::complex<FloatType>*)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  lapack_int lapack_getri(int, lapack_int, std::complex<FloatType>*,
    lapack_int, std::complex<FloatType>*)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  void syr(int, int, int, FloatType, const FloatType*, int,
    FloatType*, int)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  void syrk(int, int, int, int, int, FloatType, const FloatType *, int,
    FloatType, FloatType *, int)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }

  template <typename FloatType>
  void gemm(int, int, int, int, int, int, FloatType, const FloatType*, int,
    const FloatType*, int, FloatType, FloatType*, int)
  {
    throw SCITBX_NOT_IMPLEMENTED();
  }
}

#endif // USE_FAST_LINALG
