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

namespace cctbx { namespace fftbx {
  namespace detail {

    template <class VectorType>
    class array_2d_tp
    {
      public:
        array_2d_tp(typename VectorType::iterator Start,
                    std::size_t Nx,
                    std::size_t Ny)
          : m_Start(Start),
            m_Nx(Nx),
            m_Ny(Ny) {}
        typename VectorType::value_type&
        operator()(std::size_t ix, std::size_t iy) {
          return m_Start[iy * m_Nx + ix];
        }
      private:
        typename VectorType::iterator m_Start;
        std::size_t m_Nx;
        std::size_t m_Ny;
    };

    template <class VectorType>
    class array_3d_tp
    {
      public:
        array_3d_tp(typename VectorType::iterator Start,
                    std::size_t Nx,
                    std::size_t Ny,
                    std::size_t Nz)
          : m_Start(Start),
            m_Nx(Nx),
            m_Ny(Ny),
            m_Nz(Nz) {}
        typename VectorType::value_type&
        operator()(std::size_t ix, std::size_t iy, std::size_t iz) {
          return m_Start[(iz * m_Ny + iy) * m_Nx + ix];
        }
      private:
        typename VectorType::iterator m_Start;
        std::size_t m_Nx;
        std::size_t m_Ny;
        std::size_t m_Nz;
    };

  } // namespace detail
}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_DETAIL_ADAPTORS_H
