// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H
#define CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H

#include <algorithm>
#include <cctbx/fftbx/complex_to_complex.h>
#include <cctbx/fftbx/real_to_complex.h>
#include <cctbx/fftbx/adaptors.h>

namespace cctbx { namespace fftbx {

  /*! \brief Dimensions of 3-dimensional real-to-complex array
      given dimensions of real array.
   */
  /*! The complex array contains Ncomplex.product() complex
      values, i.e. 2*Ncomplex.product() real values.
   */
  inline triple Ncomplex_from_Nreal(const triple& Nreal) {
    return triple(Nreal[0], Nreal[1], Ncomplex_from_Nreal(Nreal[2]));
  }

  namespace adaptors {

    //! 1-dimensional vector -> 3-dimensional array of real numbers.
    /*! Currently, the only difference between cc_3d and rc_3d is
        that the dimensions for rc_3d are adjusted using
        Ncomplex_from_Nreal().
     */
    template <class VectorType>
    class rc_3d : public base_3d<VectorType>
    {
      public:
        //! Constructor given an iterator and a triple of array dimensions.
        rc_3d(typename VectorType::iterator Start, const triple& Nreal)
          : base_3d<VectorType>(Start, Ncomplex_from_Nreal(Nreal)) {}
    };

  } // namespace adaptors

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
  template <class VectorType>
  class real_to_complex_3d
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef typename VectorType::iterator iterator_type;
      typedef typename VectorType::value_type value_type;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Default constructor.
      real_to_complex_3d() {}
      //! Initialization for transforms of lengths Nreal.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(const boost::array<std::size_t, 3>& Nreal)
        : m_Nreal(Nreal) {
        init();
      }
      //! Initialization for transforms of lengths Nx, Ny, Nz.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(std::size_t Nx, std::size_t Ny, std::size_t Nz)
        : m_Nreal(Nx, Ny, Nz) {
        init();
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex, real_to_complex
       */
      void forward(adaptors::rc_3d<VectorType>& Map);
      //! In-place "backward" Fourier transformation.
      /*! <b>It is important to provide both
          Map(i,j,0) and Map(-i,-j,0)</b>. See class details.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      void backward(adaptors::rc_3d<VectorType>& Map);
      //! In-place "forward" Fourier transformation.
      /*! It is assumed that the dimensions of the Map are equal to
          Nreal.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      void forward(VectorType& Map) {
        adaptors::rc_3d<VectorType> a(Map.begin(), m_Nreal);
        forward(a);
      }
      //! In-place "backward" Fourier transformation.
      /*! It is assumed that the dimensions of the Map are equal to
          Ncomplex_from_Nreal().
          <p>
          <b>It is important to provide both
          Map(i,j,0) and Map(-i,-j,0)</b>. See class details.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      void backward(VectorType& Map) {
        adaptors::rc_3d<VectorType> a(Map.begin(), m_Nreal);
        backward(a);
      }
    private:
      void init();
    private:
      triple m_Nreal;
      complex_to_complex<VectorType> m_fft1d_x;
      complex_to_complex<VectorType> m_fft1d_y;
      real_to_complex<VectorType>   m_fft1d_z;
      VectorType m_Seq;
  };

  template <class VectorType>
  void real_to_complex_3d<VectorType>::init()
  {
    m_fft1d_x = complex_to_complex<VectorType>(m_Nreal[0]);
    m_fft1d_y = complex_to_complex<VectorType>(m_Nreal[1]);
    m_fft1d_z = real_to_complex<VectorType>(m_Nreal[2]);
    m_Seq.resize(2 * Ncomplex_from_Nreal(m_Nreal).max());
  }

  template <class VectorType>
  void
  real_to_complex_3d<VectorType>::forward(adaptors::rc_3d<VectorType>& Map)
  {
    iterator_type Seq_begin = m_Seq.begin();
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.forward(Map[triple(ix, iy, 0)]);
      }
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          m_Seq[2*iy] = Map(ix, iy, 2*iz); // Seq_begin
          m_Seq[2*iy+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<forward_tag>(), Seq_begin);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map(ix, iy, 2*iz) = m_Seq[2*iy];
          Map(ix, iy, 2*iz+1) = m_Seq[2*iy+1];
        }
      }
    }
    for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          m_Seq[2*ix] = Map(ix, iy, 2*iz);
          m_Seq[2*ix+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<forward_tag>(), Seq_begin);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map(ix, iy, 2*iz) = m_Seq[2*ix];
          Map(ix, iy, 2*iz+1) = m_Seq[2*ix+1];
        }
      }
    }
  }

  template <class VectorType>
  void
  real_to_complex_3d<VectorType>::backward(adaptors::rc_3d<VectorType>& Map)
  {
    iterator_type Seq_begin = m_Seq.begin();
    for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          m_Seq[2*ix] = Map(ix, iy, 2*iz);
          m_Seq[2*ix+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<backward_tag>(), Seq_begin);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map(ix, iy, 2*iz) = m_Seq[2*ix];
          Map(ix, iy, 2*iz+1) = m_Seq[2*ix+1];
        }
      }
      for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          m_Seq[2*iy] = Map(ix, iy, 2*iz);
          m_Seq[2*iy+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<backward_tag>(), Seq_begin);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map(ix, iy, 2*iz) = m_Seq[2*iy];
          Map(ix, iy, 2*iz+1) = m_Seq[2*iy+1];
        }
      }
    }
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.backward(Map[triple(ix, iy, 0)]);
      }
    }
  }

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H
