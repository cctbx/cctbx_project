/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragment from cctbx/basic/matrix_lite.h (rwgk)
 */

#ifndef SCITBX_MATRIX_MULTIPLY_H
#define SCITBX_MATRIX_MULTIPLY_H

namespace scitbx {

  //! Generic matrix multiplication function.
  /*! AB[ma, nb] = A[ma, na] * B[na, nb]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  inline
  void
  matrix_multiply(
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

} // namespace scitbx

#endif // SCITBX_MATRIX_MULTIPLY_H
