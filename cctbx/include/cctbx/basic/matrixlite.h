/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Oct: Moved tensor transformations from adptbx (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 May: added: identidy, isDiagonal, transpose
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_MATRIXLITE_H
#define CCTBX_BASIC_MATRIXLITE_H

#include <vector>
#include <cctbx/array_family/type_holder.h>
#include <cctbx/array_family/tiny.h>

namespace cctbx {
  namespace MatrixLite {

    template <class NumType>
    void
    identity(NumType *M, std::size_t n, NumType const& diagonal = 1)
    {
      std::size_t i;
      for(i=0;i<n*n;i++) M[i] = 0;
      for(i=0;i<n*n;i+=n+1) M[i] = diagonal;
    }

    template <class NumType>
    bool
    isDiagonal(const NumType *M, const std::size_t nr, const std::size_t nc)
    {
      if (nr != nc) return false;
      for (std::size_t ir = 0; ir < nr; ir++)
        for (std::size_t ic = 0; ic < nc; ic++)
          if (ir != ic && M[ir * nc + ic]) return false;
      return true;
    }

    template <class NumType>
    void
    transpose(const NumType *M, const std::size_t nr, const std::size_t nc,
              NumType *Mt)
    {
      for (std::size_t ir = 0; ir < nr; ir++)
        for (std::size_t ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
    }

    template <class NumType>
    void
    transpose(NumType *M, const std::size_t nr, const std::size_t nc)
    {
      std::vector<NumType> Mt(nr * nc);
      for (std::size_t ir = 0; ir < nr; ir++)
        for (std::size_t ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
      for (std::size_t i = 0; i < nr * nc; i++) M[i] = Mt[i];
    }

    template <typename NumTypeA,
              typename NumTypeB,
              typename NumTypeAB>
    void
    multiply(const NumTypeA *A, const NumTypeB *B,
             const std::size_t ma,
             const std::size_t na, const std::size_t nb,
             NumTypeAB *AB) {
      // AB[ma, nb] = A[ma, na] * B[na, nb]
      for (std::size_t i = 0; i < ma; i++) {
        for (std::size_t k = 0; k < nb; k++) {
          *AB = NumTypeAB(0);
          for (std::size_t j = 0; j < na; j++) {
            *AB += NumTypeAB(A[i * na + j] * B[j * nb + k]);
          }
          AB++;
        }
      }
    }

    template <class NumType>
    af::tiny<NumType, 3>
    DiagonalElements(af::tiny_plain<NumType, 9> const& M)
    {
      af::tiny<NumType, 3> result;
      for(std::size_t i=0;i<3;i++) {
        result[i] = M[i * (3 + 1)];
      }
      return result;
    }

    template <class NumType>
    inline NumType
    Determinant(af::tiny_plain<NumType, 9> const& M) {
      return   M[0] * (M[4] * M[8] - M[5] * M[7])
             - M[1] * (M[3] * M[8] - M[5] * M[6])
             + M[2] * (M[3] * M[7] - M[4] * M[6]);
    }

    template <class NumType>
    af::tiny<NumType, 9>
    CoFactorMxTp(af::tiny_plain<NumType, 9> const& M)
    {
      af::tiny<NumType, 9> result;
      result[0] =  M[4] * M[8] - M[5] * M[7];
      result[1] = -M[1] * M[8] + M[2] * M[7];
      result[2] =  M[1] * M[5] - M[2] * M[4];
      result[3] = -M[3] * M[8] + M[5] * M[6];
      result[4] =  M[0] * M[8] - M[2] * M[6];
      result[5] = -M[0] * M[5] + M[2] * M[3];
      result[6] =  M[3] * M[7] - M[4] * M[6];
      result[7] = -M[0] * M[7] + M[1] * M[6];
      result[8] =  M[0] * M[4] - M[1] * M[3];
      return result;
    }

    template <class NumType>
    af::tiny<NumType, 3>
    cross_product(af::tiny_plain<NumType, 3> const& a,
                  af::tiny_plain<NumType, 3> const& b)
    {
      af::tiny<NumType, 3> result;
      result[0] = a[1] * b[2] - b[1] * a[2];
      result[1] = a[2] * b[0] - b[2] * a[0];
      result[2] = a[0] * b[1] - b[0] * a[1];
      return result;
    }

    template <class FloatType6, class FloatType33>
    af::tiny<FloatType33, 3*3>
    CondensedSymMx33_as_FullSymMx33(
      af::tiny<FloatType6, 6> const& Mcond,
      af::type_holder<FloatType33>)
    {
      af::tiny<FloatType33, 3*3> Mfull;
      Mfull[0] = Mcond[0];
      Mfull[1] = Mcond[3];
      Mfull[2] = Mcond[4];
      Mfull[3] = Mcond[3];
      Mfull[4] = Mcond[1];
      Mfull[5] = Mcond[5];
      Mfull[6] = Mcond[4];
      Mfull[7] = Mcond[5];
      Mfull[8] = Mcond[2];
      return Mfull;
    }

    template <class FloatType33, class FloatType6>
    inline af::tiny<FloatType6, 6>
    FullSymMx33_as_CondensedSymMx33(
      af::tiny<FloatType33, 3*3> const& Mfull,
      af::type_holder<FloatType6>)
    {
      af::tiny<FloatType6, 6> Mcond;
      Mcond[0] = Mfull[0];
      Mcond[1] = Mfull[4];
      Mcond[2] = Mfull[8];
      Mcond[3] = Mfull[1];
      Mcond[4] = Mfull[2];
      Mcond[5] = Mfull[5];
      return Mcond;
    }

    template <class FloatType>
    af::tiny<FloatType, 9>
    FullTensorTransformation(af::tiny<FloatType, 9> const& C,
                             af::tiny<FloatType, 9> const& T)
    {
      af::tiny<FloatType, 9> CT;
      multiply<FloatType>(C.begin(), T.begin(), 3, 3, 3, CT.begin());
      af::tiny<FloatType, 9> Ct;
      transpose<FloatType>(C.begin(), 3, 3, Ct.begin());
      af::tiny<FloatType, 9> CTCt;
      multiply<FloatType>(CT.begin(), Ct.begin(), 3, 3, 3, CTCt.begin());
      return CTCt;
    }

    template <class FloatTypeC, class FloatTypeT>
    inline af::tiny<FloatTypeT, 6>
    CondensedTensorTransformation(af::tiny<FloatTypeC, 9> const& C,
                                  af::tiny<FloatTypeT, 6> const& Tcond)
    {
      return
        FullSymMx33_as_CondensedSymMx33(
          FullTensorTransformation(C,
            CondensedSymMx33_as_FullSymMx33(Tcond,
              af::type_holder<FloatTypeT>())),
              af::type_holder<FloatTypeT>());
    }

    template <typename NumTypeM,
              typename NumTypeV>
    inline
    af::tiny<NumTypeV, 3>
    matrix_mul_vector(
      af::tiny<NumTypeM, 9> const& m,
      af::tiny<NumTypeV, 3> const& v)
    {
      return af::tiny<NumTypeV, 3>(
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2]);
    }

    template <typename NumTypeV,
              typename NumTypeM>
    inline
    af::tiny<NumTypeV, 3>
    vector_mul_matrix(
      af::tiny<NumTypeV, 3> const& v,
      af::tiny<NumTypeM, 9> const& m)
    {
      return af::tiny<NumTypeV, 3>(
        m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
        m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
        m[2] * v[0] + m[5] * v[1] + m[8] * v[2]);
    }

  } // namespace MatrixLite
} // namespace cctbx

#endif // CCTBX_BASIC_MATRIXLITE_H
