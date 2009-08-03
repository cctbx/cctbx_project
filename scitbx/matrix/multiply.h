#ifndef SCITBX_MATRIX_MULTIPLY_H
#define SCITBX_MATRIX_MULTIPLY_H

namespace scitbx { namespace matrix {

  //! Generic matrix multiplication function: a * b
  /*! ab[ar, bc] = a[ar, ac] * b[ac, bc]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  void
  multiply(
    const NumTypeA *a,
    const NumTypeB *b,
    unsigned ar,
    unsigned ac,
    unsigned bc,
    NumTypeAB *ab)
  {
    unsigned i0 = 0;
    for (unsigned i=0; i<ar; i++, i0+=ac) {
      for (unsigned k=0;k<bc;k++) {
        NumTypeAB s = 0;
        unsigned ij = i0;
        unsigned jk = k;
        for (unsigned j=0;j<ac;j++,jk+=bc) {
          s += a[ij++] * b[jk];
        }
        *ab++ = s;
      }
    }
  }

  //! Generic matrix multiplication function: a.transpose() * b
  /*! ab[ac, bc] = a[ar, ac].transpose() * b[ar, bc]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAtB>
  void
  transpose_multiply(
    const NumTypeA *a,
    const NumTypeB *b,
    unsigned ar,
    unsigned ac,
    unsigned bc,
    NumTypeAtB *atb)
  {
    unsigned arac = ar * ac;
    for (unsigned i=0;i<ac;i++) {
      for (unsigned k=0;k<bc;k++) {
        NumTypeAtB s = 0;
        unsigned jk = k;
        for (unsigned ji=i;ji<arac;ji+=ac,jk+=bc) {
          s += a[ji] * b[jk];
        }
        *atb++ = s;
      }
    }
  }

  //! Generic matrix multiplication function: a * b.transpose()
  /*! ab[ar, br] = a[ar, ac] * b[br, ac].transpose()
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  void
  multiply_transpose(
    const NumTypeA *a,
    const NumTypeB *b,
    unsigned ar,
    unsigned ac,
    unsigned br,
    NumTypeAB *ab)
  {
    unsigned i0 = 0;
    for (unsigned i=0; i<ar; i++, i0+=ac) {
      unsigned kj = 0;
      for (unsigned k=0;k<br;k++) {
        NumTypeAB s = 0;
        unsigned ij = i0;
        for (unsigned j=0;j<ac;j++) {
          s += a[ij++] * b[kj++];
        }
        *ab++ = s;
      }
    }
  }

  //! Generic matrix multiplication function: a * b, with b in packed_u format.
  /*! ab[ar, ac] = a[ar, ac] * b[ac, ac]
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB>
  void
  multiply_packed_u(
    const NumTypeA* a,
    const NumTypeB* b,
    unsigned ar,
    unsigned ac,
    NumTypeAB *ab)
  {
    unsigned i0 = 0;
    for(unsigned i=0;i<ar;i++,i0+=ac) {
      for(unsigned k=0;k<ac;k++) {
        NumTypeAB s = 0;
        unsigned ij = i0;
        unsigned jk = k;
        for(unsigned j=0;j<k;jk+=ac-(++j)) {
          s += a[ij++] * b[jk];
        }
        for(unsigned j=k;j<ac;j++) {
          s += a[ij++] * b[jk++];
        }
        *ab++ = s;
      }
    }
  }

  /*! \brief Generic matrix multiplication function:
      a * b * a.transpose(), with b in packed_u format.
   */
  /*! size of b is ac*(ac+1)/2;
      size of ab is ar*ac;
      size of abat is ar*(ar+1)/2
   */
  template <typename NumTypeA,
            typename NumTypeB,
            typename NumTypeAB,
            typename NumTypeABAt>
  void
  multiply_packed_u_multiply_lhs_transpose(
    const NumTypeA* a,
    const NumTypeB* b,
    unsigned ar,
    unsigned ac,
    NumTypeAB *ab,
    NumTypeABAt *abat)
  {
    multiply_packed_u(a, b, ar, ac, ab);
    unsigned i0 = 0;
    for(unsigned i=0;i<ar;i++,i0+=ac) {
      unsigned k0 = i0;
      for(unsigned k=i;k<ar;k++,k0+=ac) {
        NumTypeABAt s = 0;
        unsigned ij = i0;
        unsigned kj = k0;
        for(unsigned j=0;j<ac;j++) {
          s += ab[ij++] * a[kj++];
        }
        *abat++ = s;
      }
    }
  }

  //! Generic matrix multiplication function: a.transpose() * a
  /*! ata[ac, ac] = a[ar, ac].transpose() * a[ar, ac]
   */
  template <typename NumTypeA,
            typename NumTypeAB>
  void
  transpose_multiply_as_packed_u(
    const NumTypeA *a,
    unsigned ar,
    unsigned ac,
    NumTypeAB *ata)
  {
    if (ar == 0) {
      std::fill(ata, ata+(ac*(ac+1)/2), static_cast<NumTypeAB>(0));
      return;
    }
    std::size_t ik = 0;
    for(unsigned i=0;i<ac;i++) {
      for(unsigned k=i;k<ac;k++) {
        ata[ik++] = a[i] * a[k];
      }
    }
    unsigned jac = ac;
    for(unsigned j=1;j<ar;j++,jac+=ac) {
      ik = 0;
      for(unsigned i=0;i<ac;i++) {
        for(unsigned k=i;k<ac;k++) {
          ata[ik++] += a[jac + i] * a[jac + k];
        }
      }
    }
  }

  //! Generic matrix multiplication function: a.transpose() * d * a
  /*! a = n x n square matrix
      d = n diagonal elements of diagonal matrix
   */
  template <typename NumTypeA,
            typename NumTypeD,
            typename NumTypeAD>
  void
  transpose_multiply_diagonal_multiply_as_packed_u(
    const NumTypeA *a,
    const NumTypeD *diagonal_elements,
    unsigned n,
    NumTypeAD *atda)
  {
    std::size_t ik = 0;
    for(unsigned i=0;i<n;i++) {
      NumTypeAD ad = a[i] * diagonal_elements[0];
      for(unsigned k=i;k<n;k++) {
        atda[ik++] = ad * a[k];
      }
    }
    unsigned jac = n;
    for(unsigned j=1;j<n;j++,jac+=n) {
      ik = 0;
      for(unsigned i=0;i<n;i++) {
        NumTypeAD ad = a[jac + i] * diagonal_elements[j];
        for(unsigned k=i;k<n;k++) {
          atda[ik++] += ad * a[jac + k];
        }
      }
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_MULTIPLY_H
