// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Dec 21: iterator-based interface (rwgk)
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_REAL_TO_COMPLEX_H
#define CCTBX_FFTBX_REAL_TO_COMPLEX_H

#include <complex>
#include <cctbx/fixes/cmath>
#include <cctbx/array_family/shared.h>
#include <cctbx/fftbx/factorization.h>
#include <cctbx/fftbx/detail/ref.h>

namespace cctbx { namespace fftbx {

  /*! \brief Length of a real-to-complex sequence given
      length of the real sequence.
   */
  /*! Ncomplex = Nreal / 2 + 1. I.e., for even Nreal, the
      real-to-complex transformed sequence viewed as
      consecutive real numbers (real and imaginary parts)
      contains two additional slots for the given
      floating-point type. For odd Nreal, the complex
      sequence contains one additional slot.
   */
  inline std::size_t Ncomplex_from_Nreal(std::size_t Nreal) {
    return Nreal/2+1;
  }
  /*! \brief Minimum array dimension for a real-to-complex sequence
      as number of real values, given length of the real sequence.
   */
  /*! Equivalent to 2 * Ncomplex_from_Nreal(Nreal).
   */
  inline std::size_t Mreal_from_Nreal(std::size_t Nreal) {
    return 2 * Ncomplex_from_Nreal(Nreal);
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
      //! Initialization for transforms of length Nreal.
      /*! This constructor determines the factorization of Nreal,
          pre-computes some constants, determines the
          "twiddle factors" needed in the transformation
          (Nreal RealType values), and allocates scratch space
          (Nreal RealType values).
       */
      real_to_complex(std::size_t Nreal);
      //! Length of real sequence.
      /*! See also: Ncomplex_from_Nreal()
       */
      std::size_t Nreal() { return m_N; }
      //! Length of complex sequence.
      /*! See also: Ncomplex_from_Nreal()
       */
      std::size_t Ncomplex() { return m_Ncomplex; }
      /*! \brief Minimum array dimension for a real-to-complex sequence
          as number of real values.
       */
      /*! See also: Mreal_from_Nreal()
       */
      std::size_t Mreal() { return 2 * m_Ncomplex; }
      //! Access to the Nreal() pre-computed "twiddle factors."
      af::shared<real_type> WA() const {
        return m_WA;
      }

      /*! \brief In-place "forward" Fourier transformation of a
          sequence of Nreal() real numbers to Ncomplex()
          complex numbers.
       */
      /*! For the forward transformation, the sign in the exponent
          is "-". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void forward(ComplexOrRealIterOrPtrType Seq_begin) {
        forward_adaptor(&(*Seq_begin));
      }
      /*! \brief In-place "backward" Fourier transformation of a sequence
          of Ncomplex() complex numbers to Nreal() real numbers.
       */
      /*! For the backward transform, the sign in the exponent
          is "+". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void backward(ComplexOrRealIterOrPtrType Seq_begin) {
        backward_adaptor(&(*Seq_begin));
      }
    private:
      std::size_t m_Ncomplex;
      af::shared<real_type> m_WA;
      af::shared<real_type> m_CH;
      void forward_adaptor(complex_type* Seq_begin);
      void forward_adaptor(real_type* Seq_begin);
      void backward_adaptor(complex_type* Seq_begin);
      void backward_adaptor(real_type* Seq_begin);
      void forward_compressed(real_type* Seq_begin);
      void backward_compressed(real_type* Seq_begin);
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
                  ComplexType>::real_to_complex(std::size_t Nreal)
    : factorization(Nreal, true), m_WA(Nreal), m_CH(Nreal),
      m_Ncomplex(Ncomplex_from_Nreal(Nreal))
  {
    // Computation of the sin and cos terms.
    // Based on the second part of fftpack41/rffti1.f.
    if (m_N < 2) return;
    real_type* WA = &(*(m_WA.begin()));
    const real_type TPI = real_type(8) * std::atan(real_type(1));
    real_type ARGH = TPI / real_type(m_N);
    std::size_t IS = 0;
    std::size_t NFM1 = m_Factors.size()-1;
    std::size_t L1 = 1;
    if (NFM1 == 0) return;
    for (std::size_t K1 = 0; K1 < NFM1; K1++) {
      std::size_t IP = m_Factors[K1];
      std::size_t LD = 0;
      std::size_t L2 = L1*IP;
      std::size_t IDO = m_N/L2;
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
     sequence is also zero if Nreal() is even. FFTPACK
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
                  ComplexType>::forward_adaptor(complex_type* Seq_begin) {
    forward_adaptor(reinterpret_cast<real_type*>(Seq_begin));
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::forward_adaptor(real_type* Seq_begin)
  {
    forward_compressed(Seq_begin);
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is undone.
    std::copy_backward(Seq_begin + 1, Seq_begin + m_N, Seq_begin + m_N + 1);
    // Insert the trivial imaginary part.
    Seq_begin[1] = real_type(0);
    // If the transform length is even, the imaginary part of the
    // last complex number in the sequence is also always zero.
    // FFTPACK does not set this imaginary part. It is done here
    // instead.
    if (m_N % 2 == 0) {
      Seq_begin[m_N + 1] = real_type(0);
    }
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::backward_adaptor(complex_type* Seq_begin) {
    backward_adaptor(reinterpret_cast<real_type*>(Seq_begin));
  }

  template <typename RealType, typename ComplexType>
  void
  real_to_complex<RealType,
                  ComplexType>::backward_adaptor(real_type* Seq_begin)
  {
    // The imaginary part of the first coefficient is always zero.
    // FFTPACK uses this knowledge to conserve space: the sequence
    // of floating point numbers is shifted down one real-sized slot.
    // Here the shift is applied before calling the core transform.
    std::copy(Seq_begin + 2, Seq_begin + 2 * m_Ncomplex, Seq_begin + 1);
    backward_compressed(Seq_begin);
  }

}} // namespace cctbx::fftbx

#include <cctbx/fftbx/detail/real_to_complex_forward.h>
#include <cctbx/fftbx/detail/real_to_complex_backward.h>

#endif // CCTBX_FFTBX_REAL_TO_COMPLEX_H
