/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy from cctbx/fftbx (rwgk)
     2002 Apr: resurrected (was adaptors.h) to restore performance (rwgk)
     2001 Nov: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef SCITBX_FFTPACK_DETAIL_REF_H
#define SCITBX_FFTPACK_DETAIL_REF_H

namespace scitbx { namespace fftpack { namespace detail {

  template <class ElementType>
  class ref_2d_tp
  {
    public:
      ref_2d_tp(ElementType* Start,
                std::size_t Nx,
                std::size_t)
        : m_Start(Start),
          m_Nx(Nx) {}
      ElementType&
      operator()(std::size_t ix, std::size_t iy) {
        return m_Start[iy * m_Nx + ix];
      }
    private:
      ElementType* m_Start;
      std::size_t m_Nx;
  };

  template <class ElementType>
  class ref_3d_tp
  {
    public:
      ref_3d_tp(ElementType* Start,
                std::size_t Nx,
                std::size_t Ny,
                std::size_t)
        : m_Start(Start),
          m_Nx(Nx),
          m_Ny(Ny) {}
      ElementType&
      operator()(std::size_t ix, std::size_t iy, std::size_t iz) {
        return m_Start[(iz * m_Ny + iy) * m_Nx + ix];
      }
    private:
      ElementType* m_Start;
      std::size_t m_Nx;
      std::size_t m_Ny;
  };

}}} // namespace scitbx::fftpack::detail

#endif // SCITBX_FFTPACK_DETAIL_REF_H
