// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
     XXX update all
     XXX const NdimAccessor&
 */

#ifndef CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H
#define CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H

#include <cctbx/array.h>
#include <cctbx/fftbx/complex_to_complex.h>

namespace cctbx { namespace fftbx {

  //! 3-dimensional complex-to-complex Fast Fourier Transformation.
  template <typename RealType,
            typename ComplexType = std::complex<RealType> >
  class complex_to_complex_3d
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef RealType real_type;
      typedef ComplexType complex_type;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Convenience typedef.
      typedef cctbx::array<std::size_t, 3> triple;

      //! Default constructor.
      complex_to_complex_3d() {}
      //! Initialization for transforms of lengths N.
      /*! See also: Constructor of complex_to_complex.
       */
      complex_to_complex_3d(const boost::array<std::size_t, 3>& N);
      //! Initialization for transforms of lengths N0, N1, N2.
      /*! See also: Constructor of complex_to_complex.
       */
      complex_to_complex_3d(std::size_t N0, std::size_t N1, std::size_t N2);
      //! Access the N (or N0, N1, N2) that was passed to the constructor.
      triple N() const {
        return triple(m_fft1d[0].N(), m_fft1d[1].N(), m_fft1d[2].N());
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex
       */
      template <typename NdimAccessor>
      void forward(NdimAccessor Map) {
        transform(select_sign<forward_tag>(), Map);
      }
      //! In-place "backward" Fourier transformation.
      /*! See also: complex_to_complex
       */
      template <typename NdimAccessor>
      void backward(NdimAccessor Map) {
        transform(select_sign<backward_tag>(), Map);
      }
      //! Generic in-place Fourier transformation.
      template <typename Tag, typename NdimAccessor>
      void transform(select_sign<Tag> tag, NdimAccessor Map) {
  // FUTURE: move out of class body
  {
    complex_type* Seq = m_Seq.begin();
    for (std::size_t iz = 0; iz < m_fft1d[2].N(); iz++) {
      for (std::size_t iy = 0; iy < m_fft1d[1].N(); iy++) {
        std::size_t ix;
        for (ix = 0; ix < m_fft1d[0].N(); ix++) {
          Seq[ix] = Map[triple(ix, iy, iz)];
        }
        // Transform along x (slow direction)
        m_fft1d[0].transform(tag, Seq);
        for (ix = 0; ix < m_fft1d[0].N(); ix++) {
          Map[triple(ix, iy, iz)] = Seq[ix];
        }
      }
      for (std::size_t ix = 0; ix < m_fft1d[0].N(); ix++) {
        std::size_t iy;
        for (iy = 0; iy < m_fft1d[1].N(); iy++) {
          Seq[iy] = Map[triple(ix, iy, iz)];
        }
        // Transform along y (medium direction)
        m_fft1d[1].transform(tag, Seq);
        for (iy = 0; iy < m_fft1d[1].N(); iy++) {
          Map[triple(ix, iy, iz)] = Seq[iy];
        }
      }
    }
    for (std::size_t ix = 0; ix < m_fft1d[0].N(); ix++) {
      for (std::size_t iy = 0; iy < m_fft1d[1].N(); iy++) {
        // Transform along z (fast direction)
        m_fft1d[2].transform(tag, &Map[triple(ix, iy, 0)]);
      }
    }
  }
      }
    private:
      boost::array<complex_to_complex<real_type, complex_type>, 3> m_fft1d;
      std::vector<complex_type> m_Seq;
  };

  template <typename RealType, typename ComplexType>
  complex_to_complex_3d<RealType, ComplexType
    >::complex_to_complex_3d(const boost::array<std::size_t, 3>& N)
    : m_Seq(cctbx::vector::max(triple(N)))
  {
    for(std::size_t i=0;i<3;i++) {
      m_fft1d[i] = complex_to_complex<real_type, complex_type>(N[i]);
    }
  }

  template <typename RealType, typename ComplexType>
  complex_to_complex_3d<RealType, ComplexType
    >::complex_to_complex_3d(std::size_t N0, std::size_t N1, std::size_t N2)
    : m_Seq(cctbx::vector::max(triple(N0, N1, N2)))
  {
    m_fft1d[0] = complex_to_complex<real_type, complex_type>(N0);
    m_fft1d[1] = complex_to_complex<real_type, complex_type>(N1);
    m_fft1d[2] = complex_to_complex<real_type, complex_type>(N2);
  }

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H
