// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_REAL_TO_COMPLEX_H
#define CCTBX_FFTBX_REAL_TO_COMPLEX_H

#include <cctbx/fixes/cmath>
#include <cctbx/fftbx/error.h>
#include <cctbx/fftbx/factorization.h>

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
      <p>
      Note that 1-dimensional real-to-complex arrays are
      "compressed" (see real_to_complex class details),
      but 3-dimensional real-to-complex arrays are not
      compressed.
   */
  inline std::size_t Ncomplex_from_Nreal(std::size_t Nreal) {
    return Nreal/2+1;
  }

  //! Real-to-complex Fast Fourier Transformation.
  /*! TODO: Currently, the 1-dimensional real-to-complex
      transformed sequences are represented according to
      "compressed" FFTPACK convention which is based on the
      observation that the imaginary part of the first
      complex number in the sequence is always zero, and
      the imaginary part of the last complex number in the
      sequence is also zero if N() is even. FFTPACK
      removes these zeros, and therefore the number of
      slots in the real and the corresonding transformed
      (complex) sequence are identical.
      <p>
      For consistency with the 3-dimensional transforms,
      ease of use and ease of understanding, the
      1-dimensional real-to-complex transforms should
      be modified to avoid the compression.
   */
  template <class VectorType>
  class real_to_complex : public factorization
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef typename VectorType::value_type value_type;
      typedef typename VectorType::iterator iterator_type;
      typedef typename VectorType::const_iterator const_iterator_type;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Default constructor.
      real_to_complex() : factorization() {}
      //! Initialization for transforms of length N.
      /*! This constructor determines the factorization of N,
          pre-computes some constants, determines the
          "twiddle factors" needed in the transformation
          (N complex values), and allocates scratch space
          (N complex values).
       */
      real_to_complex(std::size_t N);
      //! Access to the pre-computed "twiddle factors."
      const VectorType& WA() const {
        return m_WA;
      }

      /*! \brief In-place "forward" Fourier transformation of a
          sequence of N() real numbers to Ncomplex_from_Nreal()
          complex numbers.
       */
      /*! For the forward transformation, the sign in the exponent
          is "-". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      void forward(VectorType& Seq) {
        if (Seq.size() < m_N) {
          throw error("Input sequence is too short.");
        }
        forward(Seq.begin());
      }
      /*! \brief In-place "backward" Fourier transformation of a
          sequence of Ncomplex_from_Nreal() complex
          numbers to N() real numbers.
       */
      /*! For the backward transform, the sign in the exponent
          is "+". See also class complex_to_complex.
          <p>
          See also: class details.
       */
      void backward(VectorType& Seq) {
        if (Seq.size() < m_N) {
          throw error("Input sequence is too short.");
        }
        forward(Seq.begin());
      }
      /*! \brief In-place "forward" Fourier transformation of a
          sequence of N() real numbers to Ncomplex_from_Nreal()
          complex numbers.
       */
      void forward(iterator_type Seq_begin);
      /*! \brief In-place "backward" Fourier transformation of a
          sequence of Ncomplex_from_Nreal() complex
          numbers to N() real numbers.
       */
      void backward(iterator_type Seq_begin);
    private:
      VectorType m_WA;
      VectorType m_CH;
      void passf2(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1);
      void passf3(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2);
      void passf4(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2,
                  const_iterator_type WA3);
      void passf5(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2,
                  const_iterator_type WA3,
                  const_iterator_type WA4);
      void passfg(std::size_t IDO,
                  std::size_t IP,
                  std::size_t L1,
                  std::size_t IDL1,
                  iterator_type CC_start,
                  iterator_type C1_start,
                  iterator_type C2_start,
                  iterator_type CH_start,
                  iterator_type CH2_start,
                  const_iterator_type WA);
      void passb2(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1);
      void passb3(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2);
      void passb4(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2,
                  const_iterator_type WA3);
      void passb5(std::size_t IDO,
                  std::size_t L1,
                  iterator_type CC_start,
                  iterator_type CH_start,
                  const_iterator_type WA1,
                  const_iterator_type WA2,
                  const_iterator_type WA3,
                  const_iterator_type WA4);
      void passbg(std::size_t IDO,
                  std::size_t IP,
                  std::size_t L1,
                  std::size_t IDL1,
                  iterator_type CC_start,
                  iterator_type C1_start,
                  iterator_type C2_start,
                  iterator_type CH_start,
                  iterator_type CH2_start,
                  const_iterator_type WA);
  };

  template <class VectorType>
  real_to_complex<VectorType>::real_to_complex(std::size_t N)
    : factorization(N, true), m_WA(N), m_CH(N)
  {
    // Computation of the sin and cos terms.
    // Based on the second part of fftpack41/rffti1.f.
    if (m_N < 2) return;
    const value_type TPI = value_type(8) * std::atan(value_type(1));
    value_type ARGH = TPI / value_type(m_N);
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
        value_type ARGLD = value_type(LD)*ARGH;
        value_type FI = 0.;
        for (std::size_t II = 3; II <= IDO; II += 2) {
          I = I+2;
          FI = FI+1.;
          value_type ARG = FI*ARGLD;
          m_WA[I-1-1] = std::cos(ARG);
          m_WA[I-1] = std::sin(ARG);
        }
        IS = IS+IDO;
      }
      L1 = L2;
    }
  }

}} // namespace cctbx::fftbx

#include <cctbx/fftbx/detail/real_to_complex_forward.h>
#include <cctbx/fftbx/detail/real_to_complex_backward.h>

#endif // CCTBX_FFTBX_REAL_TO_COMPLEX_H
