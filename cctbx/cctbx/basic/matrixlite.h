// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_MATRIXLITE_H
#define CCTBX_BASIC_MATRIXLITE_H

#include <cstddef>
#include <boost/array.hpp>

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
  void transpose(const T *M, const std::size_t nr, const std::size_t nc, T *Mt)
  {
    for (int ir = 0; ir < nr; ir++)
      for (int ic = 0; ic < nc; ic++)
        Mt[ic * nr + ir] = M[ir * nc + ic];
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
}

#endif // CCTBX_BASIC_MATRIXLITE_H
