// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H
#define CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H

#include <boost/array.hpp>
#include <cctbx/fftbx/complex_to_complex.h>
#include <cctbx/fftbx/adaptors.h>

namespace cctbx { namespace fftbx {

  namespace adaptors {

    //! 1-dimensional vector -> 3-dimensional array of complex numbers.
    /*! Currently the complex numbers are represented by two
        consecutive real numbers (real part and imaginary part).
        It is planned that this will change in the future
        (i.e. std::complex<> will be used).
     */
    template <class VectorType>
    class cc_3d : public base_3d<VectorType>
    {
      public:
        //! Constructor given an iterator and a triple of array dimensions.
        cc_3d(typename VectorType::iterator Start, const triple& N)
          : base_3d<VectorType>(Start, N) {}
    };

  } // namespace adaptors

  //! 3-dimensional complex-to-complex Fast Fourier Transformation.
  template <class VectorType>
  class complex_to_complex_3d
  {
    public:
      //! Default constructor.
      complex_to_complex_3d() {}
      //! Initialization for transforms of lengths N.
      /*! See also: Constructor of complex_to_complex.
       */
      complex_to_complex_3d(const boost::array<std::size_t, 3>& N);
      //! Initialization for transforms of lengths Nx, Ny, Nz.
      /*! See also: Constructor of complex_to_complex.
       */
      complex_to_complex_3d(std::size_t Nx, std::size_t Ny, std::size_t Nz);
      //! Access the N (or Nx, Ny, Nz) that was passed to the constructor.
      triple N() const {
        return triple(m_fft1d[0].N(), m_fft1d[1].N(), m_fft1d[2].N());
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex
       */
      void forward(adaptors::cc_3d<VectorType>& Map);
      //! In-place "backward" Fourier transformation.
      /*! See also: complex_to_complex
       */
      void backward(adaptors::cc_3d<VectorType>& Map);
      //! In-place "forward" Fourier transformation.
      /*! It is assumed that the dimensions of the Map are equal to N().
          See also: complex_to_complex
       */
      void forward(VectorType& Map) {
        adaptors::cc_3d<VectorType> a(
          Map.begin(),
          triple(m_fft1d[0].N(), m_fft1d[1].N(), m_fft1d[2].N()));
        transform(select_sign<forward_tag>(), a);
      }
      //! In-place "backward" Fourier transformation.
      /*! It is assumed that the dimensions of the Map are equal to N().
          See also: complex_to_complex
       */
      void backward(VectorType& Map) {
        adaptors::cc_3d<VectorType> a(
          Map.begin(),
          triple(m_fft1d[0].N(), m_fft1d[1].N(), m_fft1d[2].N()));
        transform(select_sign<backward_tag>(), a);
      }
      //! Generic in-place Fourier transformation.
      template <class Tag>
      void transform(select_sign<Tag>,
                     adaptors::cc_3d<VectorType>& Map) {
  // FUTURE: move out of class body
  {
    typename VectorType::iterator Seq_begin = m_Seq.begin();
    for (std::size_t iz = 0; iz < m_fft1d[2].N(); iz++) {
      for (std::size_t iy = 0; iy < m_fft1d[1].N(); iy++) {
        std::size_t ix;
        for (ix = 0; ix < m_fft1d[0].N(); ix++) {
          m_Seq[2*ix] = Map(ix, iy, 2*iz);
          m_Seq[2*ix+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along x (slow direction)
        m_fft1d[0].transform(select_sign<Tag>(), Seq_begin);
        for (ix = 0; ix < m_fft1d[0].N(); ix++) {
          Map(ix, iy, 2*iz) = m_Seq[2*ix];
          Map(ix, iy, 2*iz+1) = m_Seq[2*ix+1];
        }
      }
      for (std::size_t ix = 0; ix < m_fft1d[0].N(); ix++) {
        std::size_t iy;
        for (iy = 0; iy < m_fft1d[1].N(); iy++) {
          m_Seq[2*iy] = Map(ix, iy, 2*iz);
          m_Seq[2*iy+1] = Map(ix, iy, 2*iz+1);
        }
        // Transform along y (medium direction)
        m_fft1d[1].transform(select_sign<Tag>(), Seq_begin);
        for (iy = 0; iy < m_fft1d[1].N(); iy++) {
          Map(ix, iy, 2*iz) = m_Seq[2*iy];
          Map(ix, iy, 2*iz+1) = m_Seq[2*iy+1];
        }
      }
    }
    for (std::size_t ix = 0; ix < m_fft1d[0].N(); ix++) {
      for (std::size_t iy = 0; iy < m_fft1d[1].N(); iy++) {
        // Transform along z (fast direction)
        m_fft1d[2].transform(select_sign<Tag>(), Map[triple(ix, iy, 0)]);
      }
    }
  }
      }
    private:
      boost::array<complex_to_complex<VectorType>, 3> m_fft1d;
      VectorType m_Seq;
  };

  template <class VectorType>
  complex_to_complex_3d<VectorType>::complex_to_complex_3d(
    const boost::array<std::size_t, 3>& N)
  {
    for(std::size_t i=0;i<3;i++) {
      m_fft1d[i] = complex_to_complex<VectorType>(N[i]);
    }
    m_Seq.resize(2 * triple(N).max());
  }

  template <class VectorType>
  complex_to_complex_3d<VectorType>::complex_to_complex_3d(
    std::size_t Nx, std::size_t Ny, std::size_t Nz)
  {
    m_fft1d[0] = complex_to_complex<VectorType>(Nx);
    m_fft1d[1] = complex_to_complex<VectorType>(Ny);
    m_fft1d[2] = complex_to_complex<VectorType>(Nz);
    m_Seq.resize(2 * triple(Nx, Ny, Nz).max());
  }


}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_COMPLEX_TO_COMPLEX_3D_H
