/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy from cctbx/fftbx (rwgk)
     2001 Dec: iterator-based interface (rwgk)
     2001 Nov: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef SCITBX_FFTPACK_REAL_TO_COMPLEX_H
#define SCITBX_FFTPACK_REAL_TO_COMPLEX_H

#include <complex>
#include <cmath>
#include <scitbx/array_family/shared.h>
#include <scitbx/fftpack/factorization.h>
#include <scitbx/fftpack/detail/ref.h>

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
        return WA_;
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
      void forward(ComplexOrRealIterOrPtrType seq_begin)
      {
        forward_adaptor(&(*seq_begin));
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
      void backward(ComplexOrRealIterOrPtrType seq_begin)
      {
        backward_adaptor(&(*seq_begin));
      }
    private:
      std::size_t n_complex_;
      af::shared<real_type> WA_;
      af::shared<real_type> CH_;
      void forward_adaptor(complex_type* seq_begin);
      void forward_adaptor(real_type* seq_begin);
      void backward_adaptor(complex_type* seq_begin);
      void backward_adaptor(real_type* seq_begin);
      void forward_compressed(real_type* seq_begin);
      void backward_compressed(real_type* seq_begin);
      void passf2(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1);
      void passf3(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2);
      void passf4(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2,
                  const real_type* WA3);
      void passf5(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2,
                  const real_type* WA3,
                  const real_type* WA4);
      void passfg(std::size_t IDO,
                  std::size_t IP,
                  std::size_t L1,
                  std::size_t IDL1,
                  real_type* CC_start,
                  real_type* C1_start,
                  real_type* C2_start,
                  real_type* CH_start,
                  real_type* CH2_start,
                  const real_type* WA);
      void passb2(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1);
      void passb3(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2);
      void passb4(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2,
                  const real_type* WA3);
      void passb5(std::size_t IDO,
                  std::size_t L1,
                  real_type* CC_start,
                  real_type* CH_start,
                  const real_type* WA1,
                  const real_type* WA2,
                  const real_type* WA3,
                  const real_type* WA4);
      void passbg(std::size_t IDO,
                  std::size_t IP,
                  std::size_t L1,
                  std::size_t IDL1,
                  real_type* CC_start,
                  real_type* C1_start,
                  real_type* C2_start,
                  real_type* CH_start,
                  real_type* CH2_start,
                  const real_type* WA);
  };

  template <typename RealType, typename ComplexType>
  real_to_complex<RealType,
                  ComplexType>::real_to_complex(std::size_t n_real)
    : factorization(n_real, true),
      n_complex_(n_complex_from_n_real(n_real)),
      WA_(n_real),
      CH_(n_real)
  {
    // Computation of the sin and cos terms.
    // Based on the second part of fftpack41/rffti1.f.
    if (n_ < 2) return;
    real_type* WA = &(*(WA_.begin()));
    const real_type TPI = real_type(8) * std::atan(real_type(1));
    real_type ARGH = TPI / real_type(n_);
    std::size_t IS = 0;
    std::size_t NFM1 = factors_.size()-1;
    std::size_t L1 = 1;
    if (NFM1 == 0) return;
    for (std::size_t K1 = 0; K1 < NFM1; K1++) {
      std::size_t IP = factors_[K1];
      std::size_t LD = 0;
      std::size_t L2 = L1*IP;
      std::size_t IDO = n_/L2;
      std::size_t IPM = IP-1;
      for (std::size_t J = 1; J <= IPM; J++) {
        LD = LD+L1;
        std::size_t I = IS;
        real_type ARGLD = real_type(LD)*ARGH;
        real_type FI = 0.;
        for (std::size_t II = 3; II <= IDO; II += 2) {
          I = I+2;
          FI = FI+1.;
          real_type ARG = FI*ARGLD;
          WA[I-1-1] = std::cos(ARG);
          WA[I-1] = std::sin(ARG);
        }
        IS = IS+IDO;
      }
      L1 = L2;
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
  real_to_complex<RealType,
                  ComplexType>::forward_adaptor(complex_type* seq_begin)
  {
    forward_adaptor(reinterpret_cast<real_type*>(seq_begin));
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::forward_adaptor(real_type* seq_begin)
  {
    forward_compressed(seq_begin);
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is undone.
    std::copy_backward(seq_begin + 1, seq_begin + n_, seq_begin + n_ + 1);
    // Insert the trivial imaginary part.
    seq_begin[1] = real_type(0);
    // If the transform length is even, the imaginary part of the
    // last complex number in the sequence is also always zero.
    // FFTPACK does not set this imaginary part. It is done here
    // instead.
    if (n_ % 2 == 0) {
      seq_begin[n_ + 1] = real_type(0);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::backward_adaptor(complex_type* seq_begin)
  {
    backward_adaptor(reinterpret_cast<real_type*>(seq_begin));
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::backward_adaptor(real_type* seq_begin)
  {
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is applied before calling the core transform.
    std::copy(seq_begin + 2, seq_begin + 2 * n_complex_, seq_begin + 1);
    backward_compressed(seq_begin);
  }

}} // namespace scitbx::fftpack

#include <scitbx/fftpack/detail/real_to_complex_forward.h>
#include <scitbx/fftpack/detail/real_to_complex_backward.h>

#endif // SCITBX_FFTPACK_REAL_TO_COMPLEX_H
