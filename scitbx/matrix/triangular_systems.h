#ifndef SCITBX_MATRIX_TRIANGULAR_SYSTEMS_H
#define SCITBX_MATRIX_TRIANGULAR_SYSTEMS_H

namespace scitbx { namespace matrix {

  /// Forward substitution
  /** Solve L x = b for a n x n lower diagonal matrix L.
      L is stored packed by rows, i.e. if L is
        [ L00 0   0   0   ]
        [ L10 L11 0   0   ]
        [ L20 L21 L22 0   ]
        [ L30 L31 L32 L33 ]
      then l points at an array
        [ L00 L10 L11 L20 L21 L22 L30 L31 L32 L33 ]
      This is refered to as the packed L storage throughout the scitbx
      b is overwritten by the solution x.
      unit_diag specifies whether the diagonal of L contains only zeroes.

   Reference: Algorithm 3.1.1 in Golub & Van Loan
  */
  template <typename FloatType>
  void forward_substitution(int n, FloatType const *l, FloatType *b,
                            bool unit_diag=false)
  {
    for (int i=0; i<n; ++i, ++l) {
      for (int j=0; j<i; ++j, ++l) {
        FloatType l_ij = *l;
        b[i] -= l_ij * b[j];
      }
      if (!unit_diag) {
        FloatType l_ii = *l;
        b[i] /= l_ii;
      }
    }
  }

  /// Backward substitution given the transpose
  /** Solve L^T x = b for otherwise the same arguments as
      matrix::forward_substitution. This is the same as
      U x = b with U upper diagonal stored packed by columns.

   Reference: Algorithm 3.1.4 in Golub & Van Loan
  */
  template <typename FloatType>
  void back_substitution_given_transpose(int n, FloatType const *l, FloatType *b,
                                         bool unit_diag=false)
  {
    FloatType const *u = l + n*(n+1)/2 - 1;
    for (int j=n-1; j>=0; --j) {
      if (!unit_diag) {
        FloatType u_jj = *u;
        b[j] /= u_jj;
      }
      --u;
      for (int i=j-1; i>=0; --i, --u) {
        FloatType u_ij = *u;
        b[i] -= b[j] * u_ij;
      }
    }
  }

  /// Back substitution
  /** Solve U x = b for a n x n upper diagonal matrix U.
      U is stored packed by rows, i.e. if U is
        [ R00 R01 R02 R03 ]
        [  0  R11 R12 R13 ]
        [  0   0  R22 R23 ]
        [  0   0   0  R33 ]
      Then u points at an array
        [ R00 R01 R02 R03 R11 R12 R13 R22 R23 R33 ]
      This is refered to as the packed U storage throughout the scitbx
      b is overwritten by the solution x.
      unit_diag specifies whether the diagonal of U contains only zeroes.

   Reference: Algorithm 3.1.2 in Golub & Van Loan
  */
  template <typename FloatType>
  void back_substitution(int n, FloatType const *u, FloatType *b,
                         bool unit_diag=false)
  {
    u += n*(n+1)/2 - 1;
    for (int i=n-1; i>=0; --i, --u) {
      for (int j=n-1; j>=i+1; --j, --u) {
        FloatType u_ij = *u;
        b[i] -= u_ij * b[j];
      }
      if (!unit_diag) {
        FloatType u_ii = *u;
        b[i] /= u_ii;
      }
    }
  }

  /// Forward substitution given the transpose
  /** Solve U^T x = b for otherwise the same arguments as
   matrix::back_substitution.  This is the same as
   L x = b with L lower diagonal stored packed by columns.

   Reference: Algorithm 3.1.3 in Golub & Van Loan
  */
  template <typename FloatType>
  void forward_substitution_given_transpose(int n,
                                            FloatType const *u, FloatType *b,
                                            bool unit_diag=false)
  {
    FloatType const *l = u;
    for (int j=0; j<n; ++j) {
      if (!unit_diag) {
        FloatType l_jj = *l;
        b[j] /= l_jj;
      }
      ++l;
      for (int i=j+1; i<n; ++i, ++l) {
        FloatType l_ij = *l;
        b[i] -= b[j] * l_ij;
      }
    }
  }

}}


#endif // GUARD
