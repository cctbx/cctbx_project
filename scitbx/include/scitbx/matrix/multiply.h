#ifndef SCITBX_MATRIX_MULTIPLY_H
#define SCITBX_MATRIX_MULTIPLY_H

namespace scitbx { namespace matrix {

  //! Generic matrix multiplication function.
  /*! AB[ma, nb] = A[ma, na] * B[na, nb]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  inline
  void
  multiply(
    const NumTypeA *A,
    const NumTypeB *B,
    std::size_t ma,
    std::size_t na,
    std::size_t nb,
    NumTypeAB *AB)
  {
    for (std::size_t i=0;i<ma;i++) {
      for (std::size_t k=0;k<nb;k++) {
        *AB = NumTypeAB(0);
        for (std::size_t j=0;j<na;j++) {
          *AB += NumTypeAB(A[i*na+j] * B[j*nb+k]);
        }
        AB++;
      }
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_MULTIPLY_H
