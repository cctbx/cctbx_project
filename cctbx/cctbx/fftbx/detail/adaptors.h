// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_DETAIL_ADAPTORS_H
#define CCTBX_FFTBX_DETAIL_ADAPTORS_H

#include <boost/array.hpp>
#include <cctbx/ndim.h>

namespace cctbx { namespace fftbx {
  namespace detail {

    template <class VectorType, std::size_t D>
    class array_tp
    {
      public:
        template <typename SizeTypeX, typename SizeTypeY>
        array_tp(typename VectorType::iterator Start,
                 const SizeTypeX& Nx,
                 const SizeTypeY& Ny)
          : m_Start(Start) {
          boost::array<std::size_t, D> N = { Nx, Ny };
          m_N = N;
        }
        template <typename SizeTypeX, typename SizeTypeY, typename SizeTypeZ>
        array_tp(typename VectorType::iterator Start,
                 const SizeTypeX& Nx,
                 const SizeTypeY& Ny,
                 const SizeTypeZ& Nz)
          : m_Start(Start) {
          boost::array<std::size_t, D> N = { Nx, Ny, Nz };
          m_N = N;
        }
        template <typename SizeTypeX, typename SizeTypeY>
        typename VectorType::value_type&
        operator()(const SizeTypeX& ix,
                   const SizeTypeY& iy) {
          boost::array<std::size_t, D> I = { ix, iy };
          return m_Start[fortran_index_1d<D>()(m_N, I)];
        }
        template <typename SizeTypeX, typename SizeTypeY, typename SizeTypeZ>
        typename VectorType::value_type&
        operator()(const SizeTypeX& ix,
                   const SizeTypeY& iy,
                   const SizeTypeZ& iz) {
          boost::array<std::size_t, D> I = { ix, iy, iz };
          return m_Start[fortran_index_1d<D>()(m_N, I)];
        }
      private:
        typename VectorType::iterator m_Start;
        boost::array<std::size_t, D> m_N;
    };

  } // namespace detail
}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_DETAIL_ADAPTORS_H
