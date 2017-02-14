#ifndef SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_H
#define SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_H

#include <scitbx/fftpack/complex_to_complex.h>
#include <scitbx/fftpack/real_to_complex.h>
#include <omptbx/omp_or_stubs.h>

#define SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP

namespace scitbx { namespace fftpack {

  /*! \brief Physical dimensions of 3-dimensional real-to-complex array
      as complex array, given generic dimensions of real array.
   */
  /*! The real-to-complex array contains product(n_complex) complex
      values, i.e. product(2*n_complex) real values.
      <p>
      See also: m_real_from_n_real()
   */
  template <typename IntegerType, std::size_t D>
  inline af::tiny<IntegerType, D>
  n_complex_from_n_real(const af::tiny<IntegerType, D>& n_real)
  {
    af::tiny<IntegerType, D> result = n_real;
    result[D-1] = n_complex_from_n_real(result[D-1]);
    return result;
  }

  /*! \brief Physical dimensions of 3-dimensional real-to-complex array
      as real array, given generic dimensions of complex array.
   */
  /*! The real-to-complex array contains product(n_complex) complex
      values, i.e. product(2*n_complex) real values.
      <p>
      See also: n_complex_from_n_real()
   */
  template <typename IntegerType, std::size_t D>
  inline af::tiny<IntegerType, D>
  n_real_from_n_complex(const af::tiny<IntegerType, D>& n_complex)
  {
    af::tiny<IntegerType, D> result = n_complex;
    result[D-1] *= 2;
    return result;
  }

  /*! \brief Physical dimensions of 3-dimensional real-to-complex array
      as real array, given generic dimensions of real array.
   */
  /*! The real-to-complex array contains product(n_complex) complex
      values, i.e. product(2*n_complex) real values.
      <p>
      See also: n_complex_from_n_real()
   */
  template <typename IntegerType, std::size_t D>
  inline af::tiny<IntegerType, D>
  m_real_from_n_real(const af::tiny<IntegerType, D>& n_real)
  {
    af::tiny<IntegerType, D> result = n_real;
    result[D-1] = m_real_from_n_real(result[D-1]);
    return result;
  }

  //! 3-dimensional real-to-complex Fast Fourier Transformation.
  /*! The real-to-complex Fourier transform of a real array
      is Hermitian. I.e., map(i,j,k) is the conjugate complex
      of map(-i,-j,-k). Exploiting this symmetry leads to
      reduced memory usage and faster Fourier transformations.
      <p>
      In this implementation, the Hermitian symmetry is exploited
      by omitting the negative half-space in the third dimension.
      I.e., the real-to-complex transformed array contains
      only n_real/2+1 (n_real_to_n_complex()) complex values
      in the third dimension.
      <p>
      Note that sligthly more than half the data are present
      in the real-to-complex transformed array: both
      map(i,j,0) and map(-i,-j,0) are present. It would be
      impractical to remove this remaining symmetry.
      <b>For the backward transform, it is important to
      provide both map(i,j,0) and map(-i,-j,0)</b>.
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
      //! Initialization for transforms of lengths n_real.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(const af::int3& n_real)
        : n_real_(n_real)
      {
        init();
      }
      //! Initialization for transforms of lengths n0, n1, n2.
      /*! See also: Constructors of complex_to_complex and real_to_complex.
       */
      real_to_complex_3d(std::size_t n0, std::size_t n1, std::size_t n2)
        : n_real_(n0, n1, n2)
      {
        init();
      }
      //! Generic dimensions of real array.
      af::int3 n_real() const { return n_real_; }
      //! Physical dimensions of real-to-complex array as complex array.
      /*! See also: m_real(), n_complex_from_n_real()
       */
      af::int3 n_complex() const
      {
        return n_complex_from_n_real(n_real_);
      }
      //! Physical dimensions of real-to-complex array as real array.
      /*! See also: n_complex(), m_real_from_n_real()
       */
      af::int3 m_real() const
      {
        return m_real_from_n_real(n_real_);
      }
      //! In-place "forward" Fourier transformation.
      /*! See also: complex_to_complex, real_to_complex
       */
      template <typename RealOrComplexMapType>
      void forward(RealOrComplexMapType map)
      {
        typedef typename RealOrComplexMapType::value_type
          real_or_complex_type;
        forward(map, real_or_complex_type());
      }
      //! In-place "backward" Fourier transformation.
      /*! <b>It is important to provide both
          map(i,j,0) and map(-i,-j,0)</b>. See class details.
          <p>
          See also: complex_to_complex, real_to_complex
       */
      template <typename RealOrComplexMapType>
      void backward(RealOrComplexMapType map)
      {
        typedef typename RealOrComplexMapType::value_type
          real_or_complex_type;
        backward(map, real_or_complex_type());
      }
    private:
      void init();
      // Cast map of complex to map of real.
      template <typename MapType>
      void forward(MapType map, complex_type)
      {
        typedef typename MapType::accessor_type accessor_type;
        af::ref<real_type, accessor_type> rmap(
          reinterpret_cast<real_type*>(map.begin()),
          n_real_from_n_complex(map.accessor()));
        forward(rmap, real_type());
      }
      // Core routine always works on real maps.
      template <typename MapType>
      void forward(MapType map, real_type)
  // FUTURE: move out of class body
  {
    // TODO: avoid i, i+1 by casting to complex
    int nx = n_real_[0];
    int ny = n_real_[1];
    int nzc = fft1d_z_.n_complex();
    int seq_size = 2 * std::max(std::max(nx, ny), nzc);
    scitbx::auto_array<real_type> seq_and_scratch;
    if (omp_in_parallel() == 0) omp_set_dynamic(0);
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
    #pragma omp parallel
#endif
    {
      int num_threads = omp_get_num_threads();
      int i_thread = omp_get_thread_num();
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp single
#endif
      {
        seq_and_scratch = scitbx::auto_array<real_type>(
          new real_type[2 * seq_size * num_threads]);
      }
      real_type* seq = seq_and_scratch.get() + 2 * seq_size * i_thread;
      real_type* scratch = seq + seq_size;
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp for
#endif
      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          // Transform along z (fast direction)
          fft1d_z_.forward(&map(ix, iy, 0), scratch);
        }
        for (int iz = 0; iz < nzc; iz++) {
          for (int iy = 0; iy < ny; iy++) {
            seq[2*iy] = map(ix, iy, 2*iz);
            seq[2*iy+1] = map(ix, iy, 2*iz+1);
          }
          // Transform along y (medium direction)
          fft1d_y_.transform(select_sign<forward_tag>(), seq, scratch);
          for (int iy = 0; iy < ny; iy++) {
            map(ix, iy, 2*iz) = seq[2*iy];
            map(ix, iy, 2*iz+1) = seq[2*iy+1];
          }
        }
      }
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp for
#endif
      for (int iy = 0; iy < ny; iy++) {
        for (int iz = 0; iz < nzc; iz++) {
          for (int ix = 0; ix < nx; ix++) {
            seq[2*ix] = map(ix, iy, 2*iz);
            seq[2*ix+1] = map(ix, iy, 2*iz+1);
          }
          // Transform along x (slow direction)
          fft1d_x_.transform(select_sign<forward_tag>(), seq, scratch);
          for (int ix = 0; ix < nx; ix++) {
            map(ix, iy, 2*iz) = seq[2*ix];
            map(ix, iy, 2*iz+1) = seq[2*ix+1];
          }
        }
      }
    }
  }
      // Cast map of complex to map of real.
      template <typename MapType>
      void backward(MapType map, complex_type)
      {
        typedef typename MapType::accessor_type accessor_type;
        af::ref<real_type, accessor_type> rmap(
          reinterpret_cast<real_type*>(map.begin()),
          n_real_from_n_complex(map.accessor()));
        backward(rmap, real_type());
      }
      // Core routine always works on real maps.
      template <typename MapType>
      void backward(MapType map, real_type)
  // FUTURE: move out of class body
  {
    // TODO: avoid i, i+1 by casting to complex
    int nx = n_real_[0];
    int ny = n_real_[1];
    int nzc = fft1d_z_.n_complex();
    int seq_size = 2 * std::max(std::max(nx, ny), nzc);
    scitbx::auto_array<real_type> seq_and_scratch;
    if (omp_in_parallel() == 0) omp_set_dynamic(0);
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
    #pragma omp parallel
#endif
    {
      int num_threads = omp_get_num_threads();
      int i_thread = omp_get_thread_num();
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp single
#endif
      {
        seq_and_scratch = scitbx::auto_array<real_type>(
          new real_type[2 * seq_size * num_threads]);
      }
      real_type* seq = seq_and_scratch.get() + 2 * seq_size * i_thread;
      real_type* scratch = seq + seq_size;
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp for
#endif
      for (int iz = 0; iz < nzc; iz++) {
        for (int iy = 0; iy < ny; iy++) {
          for (int ix = 0; ix < nx; ix++) {
            seq[2*ix] = map(ix, iy, 2*iz);
            seq[2*ix+1] = map(ix, iy, 2*iz+1);
          }
          // Transform along x (slow direction)
          fft1d_x_.transform(select_sign<backward_tag>(), seq, scratch);
          for (int ix = 0; ix < nx; ix++) {
            map(ix, iy, 2*iz) = seq[2*ix];
            map(ix, iy, 2*iz+1) = seq[2*ix+1];
          }
        }
        for (int ix = 0; ix < nx; ix++) {
          for (int iy = 0; iy < ny; iy++) {
            seq[2*iy] = map(ix, iy, 2*iz);
            seq[2*iy+1] = map(ix, iy, 2*iz+1);
          }
          // Transform along y (medium direction)
          fft1d_y_.transform(select_sign<backward_tag>(), seq, scratch);
          for (int iy = 0; iy < ny; iy++) {
            map(ix, iy, 2*iz) = seq[2*iy];
            map(ix, iy, 2*iz+1) = seq[2*iy+1];
          }
        }
      }
#if !defined(SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_NO_PRAGMA_OMP)
      #pragma omp for
#endif
      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          // Transform along z (fast direction)
          fft1d_z_.backward(&map(ix, iy, 0), scratch);
        }
      }
    }
  }
    private:
      af::int3 n_real_;
      complex_to_complex<real_type, complex_type> fft1d_x_;
      complex_to_complex<real_type, complex_type> fft1d_y_;
      real_to_complex<real_type, complex_type>    fft1d_z_;
  };

  template <typename RealType, typename ComplexType>
  void real_to_complex_3d<RealType, ComplexType>::init()
  {
    fft1d_x_ = complex_to_complex<real_type, complex_type>(n_real_[0]);
    fft1d_y_ = complex_to_complex<real_type, complex_type>(n_real_[1]);
    fft1d_z_ = real_to_complex<real_type, complex_type>(n_real_[2]);
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_REAL_TO_COMPLEX_3D_H
