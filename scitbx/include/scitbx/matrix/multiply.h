#ifndef SCITBX_MATRIX_MULTIPLY_H
#define SCITBX_MATRIX_MULTIPLY_H

namespace scitbx { namespace matrix {

  //! Generic matrix multiplication function: a * b
  /*! ab[ma, nb] = a[ma, na] * b[na, nb]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  inline
  void
  multiply(
    const NumTypeA *a,
    const NumTypeB *b,
    unsigned ma,
    unsigned na,
    unsigned nb,
    NumTypeAB *ab)
  {
    for (unsigned i=0;i<ma;i++) {
      for (unsigned k=0;k<nb;k++) {
        *ab = static_cast<NumTypeAB>(0);
        for (unsigned j=0;j<na;j++) {
          *ab += NumTypeAB(a[i*na+j] * b[j*nb+k]);
        }
        ab++;
      }
    }
  }

  //! Generic matrix multiplication function: a.transpose() * a
  /*! ata[na, na] = a.transpose()[ma, na] * a[ma, na]
   */
  template <typename NumTypeA,
            typename NumTypeAB>
  inline
  void
  transpose_multiply_as_packed_u(
    const NumTypeA *a,
    unsigned ma,
    unsigned na,
    NumTypeAB *ata)
  {
    if (ma == 0) {
      std::fill(ata, ata+(na*(na+1)/2), static_cast<NumTypeAB>(0));
      return;
    }
    std::size_t ik = 0;
    for(unsigned i=0;i<na;i++) {
      for(unsigned k=i;k<na;k++) {
        ata[ik++] = a[i] * a[k];
      }
    }
    unsigned jna = na;
    for(unsigned j=1;j<ma;j++,jna+=na) {
      ik = 0;
      for(unsigned i=0;i<na;i++) {
        for(unsigned k=i;k<na;k++) {
          ata[ik++] += a[jna + i] * a[jna + k];
        }
      }
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_MULTIPLY_H
