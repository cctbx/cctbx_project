// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
     XXX avoid z, z+1 by casting to complex
 */

#ifndef CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H
#define CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H

#include <cctbx/array.h>
#include <cctbx/fftbx/complex_to_complex.h>
#include <cctbx/fftbx/real_to_complex.h>

namespace cctbx { namespace fftbx {

  /*! \brief Dimensions of 3-dimensional real-to-complex array
      given dimensions of real array.
   */
  /*! The complex array contains product(Ncomplex) complex
      values, i.e. product(2*Ncomplex) real values.
   */
  template <typename IntegerType, std::size_t D>
  inline boost::array<std::size_t, D>
  Ncomplex_from_Nreal(const boost::array<IntegerType, D>& Nreal) {
    boost::array<std::size_t, D> result = Nreal;
    result[D-1] = Ncomplex_from_Nreal(result[D-1]);
    return result;
  }

  //! 3-dimensional real-to-complex Fast Fourier Transformation.
  /*! The real-to-complex Fourier transform of a real array
      is Hermitian. I.e., Map(i,j,k) is the conjugate complex
      of Map(-i,-j,-k). Explointing this symmetry leads to
      reduced memory usage and faster Fourier transformations.
      <p>
      In this implementation, the Hermitian symmetry is exploited
      by omitting the negative half-space in the third dimension.
      I.e., the real-to-complex transformed array contains
      only Nreal/2+1 (Nreal_to_Ncomplex()) complex values
      in the third dimension.
      <p>
      Note that sligthly more than half the data are present
      in the real-to-complex transformed array: both
      Map(i,j,0) and Map(-i,-j,0) are present. It would be
      impractical to remove this remaining symmetry.
      <b>For the backward transform, it is important to
      provide both Map(i,j,0) and Map(-i,-j,0)</b>.
   */
  template <typename RealType,
            typename ComplexType = std::complex<RealType> >
  class real_to_complex_3d
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef RealType real_type;
      typedef ComplexType complex_type;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Convenience typedef.
      typedef cctbx::array<std::size_t, 3> triple;

      //! Default constructor.
      real_to_complex_3d() {}
      //! Initialization for transforms of lengths Nreal.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(const boost::array<std::size_t, 3>& Nreal)
        : m_Nreal(Nreal) {
        init();
      }
      //! Initialization for transforms of lengths N0, N1, N2.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(std::size_t N0, std::size_t N1, std::size_t N2)
        : m_Nreal(N0, N1, N2) {
        init();
      }
      //! XXX
      triple Nreal() const { return m_Nreal; }
      //! XXX
      triple Mreal() const {
        triple result = Ncomplex_from_Nreal(m_Nreal);
        result[2] *= 2;
        return result;
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex, real_to_complex
       */
      template <typename NdimAccessor>
      void forward(NdimAccessor Map)
  // FUTURE: move out of class body
  {
    real_type* Seq = m_Seq.begin();
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.forward(&Map[triple(ix, iy, 0)]);
      }
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Seq[2*iy] = Map[triple(ix, iy, 2*iz)];
          Seq[2*iy+1] = Map[triple(ix, iy, 2*iz+1)];
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<forward_tag>(), Seq);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map[triple(ix, iy, 2*iz)] = Seq[2*iy];
          Map[triple(ix, iy, 2*iz+1)] = Seq[2*iy+1];
        }
      }
    }
    for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Seq[2*ix] = Map[triple(ix, iy, 2*iz)];
          Seq[2*ix+1] = Map[triple(ix, iy, 2*iz+1)];
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<forward_tag>(), Seq);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map[triple(ix, iy, 2*iz)] = Seq[2*ix];
          Map[triple(ix, iy, 2*iz+1)] = Seq[2*ix+1];
        }
      }
    }
  }
      //! In-place "backward" Fourier transformation.
      /*! <b>It is important to provide both
          Map(i,j,0) and Map(-i,-j,0)</b>. See class details.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      template <typename NdimAccessor>
      void backward(NdimAccessor Map)
  // FUTURE: move out of class body
  {
    real_type* Seq = m_Seq.begin();
    for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Seq[2*ix] = Map[triple(ix, iy, 2*iz)];
          Seq[2*ix+1] = Map[triple(ix, iy, 2*iz+1)];
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<backward_tag>(), Seq);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map[triple(ix, iy, 2*iz)] = Seq[2*ix];
          Map[triple(ix, iy, 2*iz+1)] = Seq[2*ix+1];
        }
      }
      for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Seq[2*iy] = Map[triple(ix, iy, 2*iz)];
          Seq[2*iy+1] = Map[triple(ix, iy, 2*iz+1)];
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<backward_tag>(), Seq);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map[triple(ix, iy, 2*iz)] = Seq[2*iy];
          Map[triple(ix, iy, 2*iz+1)] = Seq[2*iy+1];
        }
      }
    }
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.backward(&Map[triple(ix, iy, 0)]);
      }
    }
  }
    private:
      void init();
    private:
      triple m_Nreal;
      complex_to_complex<real_type, complex_type> m_fft1d_x;
      complex_to_complex<real_type, complex_type> m_fft1d_y;
      real_to_complex<real_type, complex_type>    m_fft1d_z;
      std::vector<real_type> m_Seq;
  };

  template <typename RealType, typename ComplexType>
  void real_to_complex_3d<RealType, ComplexType>::init()
  {
    m_fft1d_x = complex_to_complex<real_type, complex_type>(m_Nreal[0]);
    m_fft1d_y = complex_to_complex<real_type, complex_type>(m_Nreal[1]);
    m_fft1d_z = real_to_complex<real_type, complex_type>(m_Nreal[2]);
    m_Seq.resize(2 * cctbx::vector::max(Ncomplex_from_Nreal(m_Nreal)));
  }

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H
