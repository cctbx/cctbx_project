// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 16: Moved tensor transformations from adptbx (rwgk)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 May 07 added: identidy, isDiagonal, transpose
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_MATRIXLITE_H
#define CCTBX_BASIC_MATRIXLITE_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <boost/array.hpp>
#include <cctbx/utils.h>

namespace boost {

    template<class AnyType, std::size_t N>
    bool operator== (const array<AnyType,N>& x, const AnyType& value) {
        for (std::size_t i = 0; i < x.size(); i++)
            if (x[i] != value) return false;
        return true;
    }

    template<class AnyType, std::size_t N>
    bool operator!= (const array<AnyType,N>& x, const AnyType& value) {
        for (std::size_t i = 0; i < x.size(); i++)
            if (x[i] != value) return true;
        return false;
    }

    template<class NumType, std::size_t N>
    NumType operator* (const array<NumType,N>& lhs,
                       const array<NumType,N>& rhs) {
        NumType result = 0;
        for (std::size_t i = 0; i < lhs.size(); i++) result += lhs[i] * rhs[i];
        return result;
    }

    template<class NumType, std::size_t N>
    const array<NumType,N> operator/ (const array<NumType,N>& lhs,
                                      const NumType& rhs) {
        array<NumType,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) result[i] = lhs[i] / rhs;
        return result;
    }

    template<class NumType, std::size_t N>
    const array<NumType,N> operator* (const NumType& lhs,
                                      const array<NumType,N>& rhs) {
        array<NumType,N> result;
        for (std::size_t i = 0; i < rhs.size(); i++) result[i] = lhs * rhs[i];
        return result;
    }

    template<class NumType, std::size_t N>
    array<NumType,N> operator+ (const array<NumType,N>& lhs,
                                const array<NumType,N>& rhs) {
        array<NumType,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) {
            result[i] = lhs[i] + rhs[i];
        }
        return result;
    }

    template<class NumType, std::size_t N>
    array<NumType,N> operator- (const array<NumType,N>& lhs,
                                const array<NumType,N>& rhs) {
        array<NumType,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) {
            result[i] = lhs[i] - rhs[i];
        }
        return result;
    }

    template<class NumType, std::size_t N>
    array<NumType,N> operator- (const array<NumType,N>& rhs) {
        array<NumType,N> result;
        for (std::size_t i = 0; i < rhs.size(); i++) {
            result[i] = -rhs[i];
        }
        return result;
    }

    template<class AnyType, std::size_t N>
    std::ostream& operator<<(std::ostream& os, const array<AnyType,N>& x) {
        os << "(";
        if (x.size() > 0) {
            for (std::size_t i = 0;;) {
                os << x[i];
                i++;
                if (i == x.size()) break;
                os << ",";
            }
        }
        os << ")";
        return os;
    }

    template <class AnyType, std::size_t N>
    std::size_t
    array_min_index(const boost::array<AnyType, N>& a) {
      return std::min_element(a.begin(), a.end()) - a.begin();
    }

    template <class AnyType, std::size_t N>
    std::size_t
    array_max_index(const boost::array<AnyType, N>& a) {
      return std::max_element(a.begin(), a.end()) - a.begin();
    }

    template <class NumType, std::size_t N>
    boost::array<NumType, N>
    array_abs(const boost::array<NumType, N>& a)
    {
      boost::array<NumType, N> result;
      for (std::size_t i = 0; i < N; i++) {
        if (a[i] < 0) result[i] = -a[i];
        else          result[i] =  a[i];
      }
      return result;
    }

} // namespace boost

namespace cctbx {
  namespace MatrixLite {

    namespace itype {

      //! Helper class for passing coordinate vectors.
      typedef boost::array<int, 3> Vec3;

      //! Helper class for passing 3x3 matrices.
      typedef boost::array<int, 3 * 3> Mx33;
    }

    namespace dtype {

      //! Helper class for passing coordinate vectors.
      typedef boost::array<double, 3> Vec3;

      //! Helper class for passing 3x3 matrices.
      typedef boost::array<double, 3 * 3> Mx33;
    }

    template <class NumType>
    void
    identity(NumType *M, const std::size_t n, const NumType& diagonal = 1)
    {
      int i;
      for(i=0;i<n*n;i++) M[i] = 0;
      for(i=0;i<n*n;i+=n+1) M[i] = diagonal;
    }

    template <class NumType>
    bool
    isDiagonal(const NumType *M, const std::size_t nr, const std::size_t nc)
    {
      if (nr != nc) return false;
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          if (ir != ic && M[ir * nc + ic]) return false;
      return true;
    }

    template <class NumType>
    void
    transpose(const NumType *M, const std::size_t nr, const std::size_t nc,
              NumType *Mt)
    {
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
    }

    template <class NumType>
    void
    transpose(NumType *M, const std::size_t nr, const std::size_t nc)
    {
      std::vector<NumType> Mt(nr * nc);
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
      for (int i = 0; i < nr * nc; i++) M[i] = Mt[i];
    }

    template <class NumType>
    void
    multiply(const NumType *A, const NumType *B,
             const std::size_t ma,
             const std::size_t na, const std::size_t nb,
             NumType *AB) {
      // AB[ma, nb] = A[ma, na] * B[na, nb]
      for (int i = 0; i < ma; i++) {
        for (int k = 0; k < nb; k++) {
          *AB = 0;
          for (int j = 0; j < na; j++) *AB += A[i * na + j] * B[j * nb + k];
          AB++;
        }
      }
    }

    template <class NumType>
    boost::array<NumType, 3>
    DiagonalElements(const boost::array<NumType, 9>& M)
    {
      boost::array<NumType, 3> result;
      for(std::size_t i=0;i<3;i++) {
        result[i] = M[i * (3 + 1)];
      }
      return result;
    }

    template <class NumType>
    inline NumType
    Determinant(const boost::array<NumType, 9>& M) {
      return   M[0] * (M[4] * M[8] - M[5] * M[7])
             - M[1] * (M[3] * M[8] - M[5] * M[6])
             + M[2] * (M[3] * M[7] - M[4] * M[6]);
    }

    template <class NumType>
    boost::array<NumType, 9>
    CoFactorMxTp(const boost::array<NumType, 9>& M)
    {
      boost::array<NumType, 9> result;
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
    boost::array<NumType, 3>
    cross_product(const boost::array<NumType, 3>& a,
                  const boost::array<NumType, 3>& b)
    {
      boost::array<NumType, 3> result;
      result[0] = a[1] * b[2] - b[1] * a[2];
      result[1] = a[2] * b[0] - b[2] * a[0];
      result[2] = a[0] * b[1] - b[0] * a[1];
      return result;
    }

    template <class FloatType, std::size_t N>
    bool
    approx_equal_scaled(const boost::array<FloatType, N>& a,
                        const boost::array<FloatType, N>& b,
                        FloatType scaled_tolerance) {
      for (std::size_t i = 0; i < N; i++) {
        if (!cctbx::approx_equal_scaled(a[i], b[i], scaled_tolerance)) {
          return false;
        }
      }
      return true;
    }

    template <class FloatType>
    struct return_type {};

    template <class FloatType6, class FloatType33>
    boost::array<FloatType33, 3*3>
    CondensedSymMx33_as_FullSymMx33(
      const boost::array<FloatType6, 6>& Mcond,
      return_type<FloatType33>)
    {
      boost::array<FloatType33, 3*3> Mfull;
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
    inline boost::array<FloatType6, 6>
    FullSymMx33_as_CondensedSymMx33(
      const boost::array<FloatType33, 3*3>& Mfull,
      return_type<FloatType6>)
    {
      boost::array<FloatType6, 6> Mcond;
      Mcond[0] = Mfull[0];
      Mcond[1] = Mfull[4];
      Mcond[2] = Mfull[8];
      Mcond[3] = Mfull[1];
      Mcond[4] = Mfull[2];
      Mcond[5] = Mfull[5];
      return Mcond;
    }

    template <class FloatType>
    boost::array<FloatType, 9>
    FullTensorTransformation(const boost::array<FloatType, 9>& C,
                             const boost::array<FloatType, 9>& T)
    {
      boost::array<FloatType, 9> CT;
      multiply<FloatType>(C.elems, T.elems, 3, 3, 3, CT.elems);
      boost::array<FloatType, 9> Ct;
      transpose<FloatType>(C.elems, 3, 3, Ct.elems);
      boost::array<FloatType, 9> CTCt;
      multiply<FloatType>(CT.elems, Ct.elems, 3, 3, 3, CTCt.elems);
      return CTCt;
    }

    template <class FloatTypeC, class FloatTypeT>
    inline boost::array<FloatTypeT, 6>
    CondensedTensorTransformation(const boost::array<FloatTypeC, 9>& C,
                                  const boost::array<FloatTypeT, 6>& Tcond)
    {
      return
        FullSymMx33_as_CondensedSymMx33(
          FullTensorTransformation(C,
            CondensedSymMx33_as_FullSymMx33(Tcond, return_type<FloatTypeT>())),
            return_type<FloatTypeT>());
    }

  } // namespace MatrixLite
} // namespace cctbx

#endif // CCTBX_BASIC_MATRIXLITE_H
