// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
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

namespace boost {

    template<class T, std::size_t N>
    bool operator== (const array<T,N>& x, const T& value) {
        for (std::size_t i = 0; i < x.size(); i++)
            if (x[i] != value) return false;
        return true;
    }

    template<class T, std::size_t N>
    bool operator!= (const array<T,N>& x, const T& value) {
        for (std::size_t i = 0; i < x.size(); i++)
            if (x[i] != value) return true;
        return false;
    }

    template<class T, std::size_t N>
    T operator* (const array<T,N>& lhs, const array<T,N>& rhs) {
        T result = 0;
        for (std::size_t i = 0; i < lhs.size(); i++) result += lhs[i] * rhs[i];
        return result;
    }

    template<class T, std::size_t N>
    const array<T,N> operator/ (const array<T,N>& lhs, const T& rhs) {
        array<T,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) result[i] = lhs[i] / rhs;
        return result;
    }

    template<class T, std::size_t N>
    const array<T,N> operator* (const T& lhs, const array<T,N>& rhs) {
        array<T,N> result;
        for (std::size_t i = 0; i < rhs.size(); i++) result[i] = lhs * rhs[i];
        return result;
    }

    template<class T, std::size_t N>
    array<T,N> operator+ (const array<T,N>& lhs, const array<T,N>& rhs) {
        array<T,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) {
            result[i] = lhs[i] + rhs[i];
        }
        return result;
    }

    template<class T, std::size_t N>
    array<T,N> operator- (const array<T,N>& lhs, const array<T,N>& rhs) {
        array<T,N> result;
        for (std::size_t i = 0; i < lhs.size(); i++) {
            result[i] = lhs[i] - rhs[i];
        }
        return result;
    }

    template<class T, std::size_t N>
    array<T,N> operator- (const array<T,N>& rhs) {
        array<T,N> result;
        for (std::size_t i = 0; i < rhs.size(); i++) {
            result[i] = -rhs[i];
        }
        return result;
    }

    template<class T, std::size_t N>
    std::ostream& operator<<(std::ostream& os, const array<T,N>& x) {
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

    template <class T>
    void identity(T *M, const std::size_t n, const T& diagonal = 1)
    {
      int i;
      for(i=0;i<n*n;i++) M[i] = 0;
      for(i=0;i<n*n;i+=n+1) M[i] = diagonal;
    }

    template <class T>
    bool isDiagonal(const T *M, const std::size_t nr, const std::size_t nc)
    {
      if (nr != nc) return false;
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          if (ir != ic && M[ir * nc + ic]) return false;
      return true;
    }

    template <class T>
    void transpose(const T *M, const std::size_t nr, const std::size_t nc,
                   T *Mt)
    {
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
    }

    template <class T>
    void transpose(T *M, const std::size_t nr, const std::size_t nc)
    {
      std::vector<T> Mt(nr * nc);
      for (int ir = 0; ir < nr; ir++)
        for (int ic = 0; ic < nc; ic++)
          Mt[ic * nr + ir] = M[ir * nc + ic];
      for (int i = 0; i < nr * nc; i++) M[i] = Mt[i];
    }

    template <class T>
    void multiply(const T *A, const T *B,
                  const std::size_t ma,
                  const std::size_t na, const std::size_t nb,
                  T *AB) {
      // AB[ma, nb] = A[ma, na] * B[na, nb]
      for (int i = 0; i < ma; i++) {
        for (int k = 0; k < nb; k++) {
          *AB = 0;
          for (int j = 0; j < na; j++) *AB += A[i * na + j] * B[j * nb + k];
          AB++;
        }
      }
    }

  } // namespace MatrixLite
} // namespace cctbx

#endif // CCTBX_BASIC_MATRIXLITE_H
