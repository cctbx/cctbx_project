// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_ADAPTORS_H
#define CCTBX_FFTBX_ADAPTORS_H

#include <cctbx/ndim.h>
#include <cctbx/fftbx/triple.h>

namespace cctbx { namespace fftbx {

  //! Adaptors for 3-dimensional transforms.
  namespace adaptors {

    //! Base class for 1-dimensional vector -> 3-dimensional array %adaptors.
    template <class VectorType>
    class base_3d
    {
      public:
        //! Constructor given an iterator and a triple of array dimensions.
        base_3d(typename VectorType::iterator Start, const triple& M)
          : m_Start(Start), m_M(M[0], M[1], 2 * M[2]) {}
        //! Iterator corresponding to the index triple.
        typename VectorType::iterator
        operator[](const triple& I) const {
          return m_Start + c_index_1d<3>()(m_M, I);
        }
        //! Reference to the array element corresponding to the index triple.
        typename VectorType::value_type&
        operator()(std::size_t ix, std::size_t iy, std::size_t iz) {
          return *(operator[](triple(ix, iy, iz)));
        }
      protected:
        typename VectorType::iterator m_Start;
        triple m_M;
    };

  } // namespace adaptors
}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_ADAPTORS_H
