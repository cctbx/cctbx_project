// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Dec 21: iterator-based interface (rwgk)
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H
#define CCTBX_FFTBX_REAL_TO_COMPLEX_3D_H

#include <cctbx/carray.h>
#include <cctbx/fftbx/complex_to_complex.h>
#include <cctbx/fftbx/real_to_complex.h>

namespace cctbx { namespace fftbx {

  /*! \brief Physical dimensions of 3-dimensional real-to-complex array
      as complex array, given generic dimensions of real array.
   */
  /*! The real-to-complex array contains product(Ncomplex) complex
      values, i.e. product(2*Ncomplex) real values.
      <p>
      See also: Mreal_from_Nreal()
   */
  template <typename IntegerType, std::size_t D>
  inline carray<IntegerType, D>
  Ncomplex_from_Nreal(const carray<IntegerType, D>& Nreal) {
    carray<IntegerType, D> result = Nreal;
    result[D-1] = Ncomplex_from_Nreal(result[D-1]);
    return result;
  }

  /*! \brief Physical dimensions of 3-dimensional real-to-complex array
      as real array, given generic dimensions of real array.
   */
  /*! The real-to-complex array contains product(Ncomplex) complex
      values, i.e. product(2*Ncomplex) real values.
      <p>
      See also: Ncomplex_from_Nreal()
   */
  template <typename IntegerType, std::size_t D>
  inline carray<IntegerType, D>
  Mreal_from_Nreal(const carray<IntegerType, D>& Nreal) {
    carray<IntegerType, D> result = Nreal;
    result[D-1] = Mreal_from_Nreal(result[D-1]);
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

      //! Default constructor.
      real_to_complex_3d() {}
      //! Initialization for transforms of lengths Nreal.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(const int3& Nreal)
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
      //! Generic dimensions of real array.
      int3 Nreal() const { return m_Nreal; }
      //! Physical dimensions of real-to-complex array as complex array.
      /*! See also: Mreal(), Ncomplex_from_Nreal()
       */
      int3 Ncomplex() const {
        return Ncomplex_from_Nreal(m_Nreal);
      }
      //! Physical dimensions of real-to-complex array as real array.
      /*! See also: Ncomplex(), Mreal_from_Nreal()
       */
      int3 Mreal() const {
        return Mreal_from_Nreal(m_Nreal);
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex, real_to_complex
       */
      template <typename RealOrComplexVecRefNd>
      void forward(RealOrComplexVecRefNd Map) {
        typedef
        typename RealOrComplexVecRefNd::value_type real_or_complex_type;
        forward(Map, real_or_complex_type());
      }
      //! In-place "backward" Fourier transformation.
      /*! <b>It is important to provide both
          Map(i,j,0) and Map(-i,-j,0)</b>. See class details.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      template <typename RealOrComplexVecRefNd>
      void backward(RealOrComplexVecRefNd Map) {
        typedef
        typename RealOrComplexVecRefNd::value_type real_or_complex_type;
        backward(Map, real_or_complex_type());
      }
    private:
      void init();
      // Cast map of complex to map of real.
      template <typename VecRefNd>
      void forward(VecRefNd Map, complex_type) {
        typedef typename VecRefNd::dimension_type dim_type;
        dim_type dim(Map.dim()[0], Map.dim()[1], Map.dim()[2] * 2);
        vecrefnd<real_type, dim_type> rmap(Map.handle(), dim);
        forward(rmap, real_type());
      }
      // Core routine always works on real maps.
      template <typename VecRefNd>
      void forward(VecRefNd Map, real_type)
  // FUTURE: move out of class body
  {
    // TODO: avoid i, i+1 by casting to complex
    real_type* Seq = &(*(m_Seq.begin()));
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.forward(&Map(int3(ix, iy, 0)));
      }
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Seq[2*iy] = Map(int3(ix, iy, 2*iz));
          Seq[2*iy+1] = Map(int3(ix, iy, 2*iz+1));
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<forward_tag>(), Seq);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map(int3(ix, iy, 2*iz)) = Seq[2*iy];
          Map(int3(ix, iy, 2*iz+1)) = Seq[2*iy+1];
        }
      }
    }
    for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
      for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Seq[2*ix] = Map(int3(ix, iy, 2*iz));
          Seq[2*ix+1] = Map(int3(ix, iy, 2*iz+1));
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<forward_tag>(), Seq);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map(int3(ix, iy, 2*iz)) = Seq[2*ix];
          Map(int3(ix, iy, 2*iz+1)) = Seq[2*ix+1];
        }
      }
    }
  }
      // Cast map of complex to map of real.
      template <typename VecRefNd>
      void backward(VecRefNd Map, complex_type) {
        typedef typename VecRefNd::dimension_type dim_type;
        dim_type dim(Map.dim()[0], Map.dim()[1], Map.dim()[2] * 2);
        vecrefnd<real_type, dim_type> rmap(Map.handle(), dim);
        backward(rmap, real_type());
      }
      // Core routine always works on real maps.
      template <typename VecRefNd>
      void backward(VecRefNd Map, real_type)
  // FUTURE: move out of class body
  {
    // TODO: avoid i, i+1 by casting to complex
    real_type* Seq = &(*(m_Seq.begin()));
    for (std::size_t iz = 0; iz < m_fft1d_z.Ncomplex(); iz++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        std::size_t ix;
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Seq[2*ix] = Map(int3(ix, iy, 2*iz));
          Seq[2*ix+1] = Map(int3(ix, iy, 2*iz+1));
        }
        // Transform along x (slow direction)
        m_fft1d_x.transform(select_sign<backward_tag>(), Seq);
        for (ix = 0; ix < m_Nreal[0]; ix++) {
          Map(int3(ix, iy, 2*iz)) = Seq[2*ix];
          Map(int3(ix, iy, 2*iz+1)) = Seq[2*ix+1];
        }
      }
      for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
        std::size_t iy;
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Seq[2*iy] = Map(int3(ix, iy, 2*iz));
          Seq[2*iy+1] = Map(int3(ix, iy, 2*iz+1));
        }
        // Transform along y (medium direction)
        m_fft1d_y.transform(select_sign<backward_tag>(), Seq);
        for (iy = 0; iy < m_Nreal[1]; iy++) {
          Map(int3(ix, iy, 2*iz)) = Seq[2*iy];
          Map(int3(ix, iy, 2*iz+1)) = Seq[2*iy+1];
        }
      }
    }
    for (std::size_t ix = 0; ix < m_Nreal[0]; ix++) {
      for (std::size_t iy = 0; iy < m_Nreal[1]; iy++) {
        // Transform along z (fast direction)
        m_fft1d_z.backward(&Map(int3(ix, iy, 0)));
      }
    }
  }
    private:
      int3 m_Nreal;
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
