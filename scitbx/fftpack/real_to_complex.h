#ifndef SCITBX_FFTPACK_REAL_TO_COMPLEX_H
#define SCITBX_FFTPACK_REAL_TO_COMPLEX_H

#include <scitbx/fftpack/factorization.h>
#include <scitbx/fftpack/detail/ref.h>
#include <scitbx/array_family/shared.h>
#include <boost/scoped_array.hpp>
#include <complex>
#include <cmath>

namespace scitbx { namespace fftpack {

  /*! \brief Length of a real-to-complex sequence given
      length of the real sequence.
   */
  /*! n_complex = n_real / 2 + 1. I.e., for even n_real, the
      real-to-complex transformed sequence viewed as
      consecutive real numbers (real and imaginary parts)
      contains two additional slots for the given
      floating-point type. For odd n_real, the complex
      sequence contains one additional slot.
   */
  inline std::size_t n_complex_from_n_real(std::size_t n_real)
  {
    return n_real/2+1;
  }
  /*! \brief Minimum array dimension for a real-to-complex sequence
      as number of real values, given length of the real sequence.
   */
  /*! Equivalent to 2 * n_complex_from_n_real(n_real).
   */
  inline std::size_t m_real_from_n_real(std::size_t n_real)
  {
    return 2 * n_complex_from_n_real(n_real);
  }

  //! Real-to-complex Fast Fourier Transformation.
  template <typename RealType,
            typename ComplexType = std::complex<RealType> >
  class real_to_complex : public factorization
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef RealType real_type;
      typedef ComplexType complex_type;
      typedef detail::ref_2d_tp<real_type> dim2;
      typedef detail::ref_3d_tp<real_type> dim3;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Default constructor.
      real_to_complex() : factorization() {}
      //! Initialization for transforms of length n_real.
      /*! This constructor determines the factorization of n_real,
          pre-computes some constants, determines the
          "twiddle factors" needed in the transformation
          (n_real RealType values), and allocates scratch space
          (n_real RealType values).
       */
      real_to_complex(std::size_t n_real);
      //! Length of real sequence.
      /*! See also: n_complex_from_n_real()
       */
      std::size_t n_real() { return n_; }
      //! Length of complex sequence.
      /*! See also: n_complex_from_n_real()
       */
      std::size_t n_complex() { return n_complex_; }
      /*! \brief Minimum array dimension for a real-to-complex sequence
          as number of real values.
       */
      /*! See also: m_real_from_n_real()
       */
      std::size_t m_real() { return 2 * n_complex_; }
      //! Access to the n_real() pre-computed "twiddle factors."
      af::shared<real_type> wa() const
      {
        return wa_;
      }

      /*! \brief In-place "forward" Fourier transformation of a
          sequence of n_real() real numbers to n_complex()
          complex numbers.
       */
      /*! For the forward transformation, the sign in the exponent
          is "-". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void
      forward(
        ComplexOrRealIterOrPtrType seq_begin,
        real_type* scratch=0)
      {
        forward_adaptor(&(*seq_begin), scratch);
      }

      /*! \brief In-place "backward" Fourier transformation of a sequence
          of n_complex() complex numbers to n_real() real numbers.
       */
      /*! For the backward transform, the sign in the exponent
          is "+". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void
      backward(
        ComplexOrRealIterOrPtrType seq_begin,
        real_type* scratch=0)
      {
        backward_adaptor(&(*seq_begin), scratch);
      }

    private:
      std::size_t n_complex_;
      af::shared<real_type> wa_;
      void forward_adaptor(complex_type* seq_begin, real_type* scratch);
      void forward_adaptor(real_type* seq_begin, real_type* scratch);
      void backward_adaptor(complex_type* seq_begin, real_type* scratch);
      void backward_adaptor(real_type* seq_begin, real_type* scratch);
      void forward_compressed(real_type* seq_begin, real_type* scratch);
      void backward_compressed(real_type* seq_begin, real_type* scratch);
      void passf2(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1);
      void passf3(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2);
      void passf4(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2,
                  const real_type* wa3);
      void passf5(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2,
                  const real_type* wa3,
                  const real_type* wa4);
      void passfg(std::size_t ido,
                  std::size_t ip,
                  std::size_t l1,
                  std::size_t idl1,
                  real_type* cc_start,
                  real_type* c1_start,
                  real_type* c2_start,
                  real_type* ch_start,
                  real_type* ch2_start,
                  const real_type* wa);
      void passb2(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1);
      void passb3(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2);
      void passb4(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2,
                  const real_type* wa3);
      void passb5(std::size_t ido,
                  std::size_t l1,
                  real_type* cc_start,
                  real_type* ch_start,
                  const real_type* wa1,
                  const real_type* wa2,
                  const real_type* wa3,
                  const real_type* wa4);
      void passbg(std::size_t ido,
                  std::size_t ip,
                  std::size_t l1,
                  std::size_t idl1,
                  real_type* cc_start,
                  real_type* c1_start,
                  real_type* c2_start,
                  real_type* ch_start,
                  real_type* ch2_start,
                  const real_type* wa);
  };

  template <typename RealType, typename ComplexType>
  real_to_complex<RealType,
                  ComplexType>::real_to_complex(std::size_t n_real)
    : factorization(n_real, true),
      n_complex_(n_complex_from_n_real(n_real)),
      wa_(n_real)
  {
    // Computation of the sin and cos terms.
    // Based on the second part of fftpack41/rffti1.f.
    if (n_ < 2) return;
    real_type* wa = &(*(wa_.begin()));
    const real_type tpi = real_type(8) * std::atan(real_type(1));
    real_type argh = tpi / real_type(n_);
    std::size_t is = 0;
    std::size_t nfm1 = factors_.size()-1;
    std::size_t l1 = 1;
    if (nfm1 == 0) return;
    for (std::size_t k1 = 0; k1 < nfm1; k1++) {
      std::size_t ip = factors_[k1];
      std::size_t ld = 0;
      std::size_t l2 = l1*ip;
      std::size_t ido = n_/l2;
      std::size_t ipm = ip-1;
      for (std::size_t j = 1; j <= ipm; j++) {
        ld = ld+l1;
        std::size_t i = is;
        real_type argld = real_type(ld)*argh;
        real_type fi = 0.;
        for (std::size_t ii = 3; ii <= ido; ii += 2) {
          i = i+2;
          fi = fi+1.;
          real_type arg = fi*argld;
          wa[i-1-1] = std::cos(arg);
          wa[i-1] = std::sin(arg);
        }
        is = is+ido;
      }
      l1 = l2;
    }
  }

  /* In the core transforms, the 1-dimensional real-to-complex
     transformed sequences are represented according to
     "compressed" FFTPACK convention which is based on the
     observation that the imaginary part of the first
     complex number in the sequence is always zero, and
     the imaginary part of the last complex number in the
     sequence is also zero if n_real() is even. FFTPACK
     removes these zeros, and therefore the number of RealType
     slots in the real and the corresonding transformed
     (complex) sequence are identical.

     For consistency with the 3-dimensional transforms,
     ease of use and ease of understanding, the
     FFTPACK compression is not exposed to the user,
     at the cost of a small runtime penalty.
   */

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::forward_adaptor(
    complex_type* seq_begin,
    real_type* scratch)
  {
    forward_adaptor(reinterpret_cast<real_type*>(seq_begin), scratch);
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::forward_adaptor(
    real_type* seq_begin,
    real_type* scratch)
  {
    if (scratch == 0) {
      boost::scoped_array<real_type> buffer(new real_type[n_]);
      scratch = buffer.get();
      forward_compressed(seq_begin, scratch);
    }
    else {
      forward_compressed(seq_begin, scratch);
    }
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is undone.
    std::copy_backward(seq_begin + 1, seq_begin + n_, seq_begin + n_ + 1);
      // Note regarding runtime overhead: in the context of 3D transforms
      // the copy was found to have no significant runtime penalty
      // (dimensions 250,300,270 used for the timings, Xeon, gcc 4.1 and 4.4).
    // Insert the trivial imaginary part.
    seq_begin[1] = real_type(0);
    // If the transform length is even, the imaginary part of the
    // last complex number in the sequence is also always zero.
    // FFTPACK does not set this imaginary part. It is done here
    // instead.
    if ((n_ & 1) == 0) {
      seq_begin[n_ + 1] = real_type(0);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::backward_adaptor(
    complex_type* seq_begin,
    real_type* scratch)
  {
    backward_adaptor(reinterpret_cast<real_type*>(seq_begin), scratch);
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType, ComplexType>::backward_adaptor(
    real_type* seq_begin,
    real_type* scratch)
  {
    if ((n_ & 1) == 0) {
      // If the transform length is even, the imaginary part of the
      // last complex number in the sequence is always zero.
      // fftpack produces the correct result even if this is not the
      // case for the input sequence, but the corresponding value is
      // used in floating-point operations. Resetting here for safety.
      seq_begin[n_ + 1] = real_type(0);
    }
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is applied before calling the core transform.
    std::copy(seq_begin + 2, seq_begin + 2 * n_complex_, seq_begin + 1);
      // See "Note regarding runtime overhead" above.
    if (scratch == 0) {
      boost::scoped_array<real_type> buffer(new real_type[n_]);
      scratch = buffer.get();
      backward_compressed(seq_begin, scratch);
    }
    else {
      backward_compressed(seq_begin, scratch);
    }
  }

}} // namespace scitbx::fftpack

#include <scitbx/fftpack/detail/real_to_complex_forward.h>
#include <scitbx/fftpack/detail/real_to_complex_backward.h>

#endif // SCITBX_FFTPACK_REAL_TO_COMPLEX_H
