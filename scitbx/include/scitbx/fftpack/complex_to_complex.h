/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy from cctbx/fftbx (rwgk)
     2001 Dec: iterator-based interface (rwgk)
     2001 Nov: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H
#define SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H

#include <algorithm>
#include <complex>
#include <cmath>
#include <scitbx/array_family/shared.h>
#include <scitbx/fftpack/factorization.h>
#include <scitbx/fftpack/detail/ref.h>

namespace scitbx {

//! Fast Fourier Transformations based on FFTPACK.
/*! The Fast Fourier Transform toolbox provides
    1-dimensional and 3-dimensional in-place
    complex-to-complex and real-to-complex transforms.
    There are no restrictions on the lengths of the
    sequences to be transformed.  However, the transforms
    are most efficient when the sequence length is a
    product of 2, 3 and 5.
    <p>
    scitbx/fftpack is a pure C++ adaption of the
    complex-to-complex and real-to-complex transforms of
    FFTPACK Version 4.1 (Nov 1988) by Paul Swarztrauber
    (see http://www.scd.ucar.edu/softlib/FFTPACK.html
    and the references below).
    <p>
    The core of scitbx/fftpack is implemented exclusively in C++
    header files. Therefore it is not necessary to link
    with a support library. Only the optional Python
    bindings that are provided for convenience need to be
    compiled.
    <p>
    scitbx/fftpack is reasonably fast on the machines where it
    was tested (Tru64 Unix, Linux, Windows with Visual C++
    6), but is not as fast as FFTW (http://www.fftw.org/).
    Depending on the machine, compiler settings, transform
    lengths and the dimensionality (1D or 3D), scitbx/fftpack is
    about 30%-300% slower than FFTW for transform lengths
    that are a multiple of 2, 3 and 5. This is no surprise,
    since FFTW is a very sophisticated library that adjusts
    itself to a particular hardware, and was developed
    specifically for modern, cache-based architectures.
    In particular, the performance on Athlon Linux machines is
    very impressive. Unfortunately, the FFTW license does
    not permit the inclusion in an Open Source project
    (the GNU General Public License is too restrictive).
    Therefore the public domain FFTPACK was chosen as
    the basis for scitbx/fftpack.
    <p>
    It is hoped that scitbx/fftpack facilitates the development
    of an elegant, generic C++ interface for Fast Fourier
    Transforms. Once such an interface is worked out, it
    should be easy to provide implementations with
    different core transform algorithms, such as FFTW.
    <p>
    <i>References</i>:
    <ul>
    <li>"Vectorizing the Fast Fourier Transforms", by Paul Swarztrauber,
        Parallel Computations, G. Rodrigue, ed., Academic Press,
        New York 1982.
    <li>"Fast Fourier Transforms Algorithms for Vector Computers", by
        Paul Swarztrauber, Parallel Computing, (1984) pp.45-63.
    </ul>
 */

namespace fftpack {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  struct forward_tag {};
  struct backward_tag {};

  template <class Tag> struct select_sign {};

  // FUTURE: use partial specialiation
  template <>
  struct select_sign<forward_tag>
  {
    template <class FloatType>
    static FloatType plusminus(const FloatType&a, const FloatType& b)
    {
      return a + b;
    }
    template <class FloatType>
    static FloatType minusplus(const FloatType&a, const FloatType& b)
    {
      return a - b;
    }
    template <class FloatType>
    static FloatType unaryplusminus(const FloatType&a)
    {
      return a;
    }
    template <class FloatType>
    static FloatType unaryminusplus(const FloatType&a)
    {
      return -a;
    }
  };

  // FUTURE: use partial specialiation
  template <>
  struct select_sign<backward_tag>
  {
    template <class FloatType>
    static FloatType plusminus(const FloatType&a, const FloatType& b)
    {
      return a - b;
    }
    template <class FloatType>
    static FloatType minusplus(const FloatType&a, const FloatType& b)
    {
      return a + b;
    }
    template <class FloatType>
    static FloatType unaryplusminus(const FloatType&a)
    {
      return -a;
    }
    template <class FloatType>
    static FloatType unaryminusplus(const FloatType&a)
    {
      return a;
    }
  };

#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! Complex-to-complex Fast Fourier Transformation.
  template <typename RealType,
            typename ComplexType = std::complex<RealType> >
  class complex_to_complex : public factorization
  {
    public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      typedef RealType real_type;
      typedef ComplexType complex_type;
      typedef detail::ref_2d_tp<real_type> dim2;
      typedef detail::ref_3d_tp<real_type> dim3;
#endif // DOXYGEN_SHOULD_SKIP_THIS

      //! Default constructor.
      complex_to_complex() : factorization() {}
      //! Initialization for transforms of length n.
      /*! This constructor determines the factorization of n,
          pre-computes some constants, determines the
          "twiddle factors" needed in the transformation
          (n complex values), and allocates scratch space
          (n complex values).
       */
      complex_to_complex(std::size_t n);
      //! Access to the pre-computed "twiddle factors."
      af::shared<real_type> wa() const
      {
        return WA_;
      }

      /*! \brief In-place "forward" Fourier transformation of a
          sequence of length n().
       */
      /*! Equivalent to, but much faster than the following Python code:<pre>
          for i in range(n):
            Sum = 0
            for k in range(n):
              Sum += Seq_in[k] * exp(-2*pi*1j*i*k/n)
            Seq_out.append(Sum)
          </pre>
          Note that 1j is the imaginary number (sqrt(-1)) in Python.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void forward(ComplexOrRealIterOrPtrType seq_begin)
      {
        transform(select_sign<forward_tag>(), &(*seq_begin));
      }
      /*! \brief In-place "backward" Fourier transformation of a
          sequence of length n().
       */
      /*! Equivalent to, but much faster than the following Python code:<pre>
          for i in range(n):
            Sum = 0
            for k in range(n):
              Sum += Seq_in[k] * exp(+2*pi*1j*i*k/n)
            Seq_out.append(Sum)
          </pre>
          Note that 1j is the imaginary number (sqrt(-1)) in Python.
       */
      template <typename ComplexOrRealIterOrPtrType>
      void backward(ComplexOrRealIterOrPtrType seq_begin)
      {
        transform(select_sign<backward_tag>(), &(*seq_begin));
      }

      //! Generic in-place Fourier transformation of a contiguous sequence.
      template <typename Tag>
      void transform(select_sign<Tag> tag, complex_type* seq_begin)
      {
        transform(tag, reinterpret_cast<real_type*>(seq_begin));
      }

      //! Generic in-place Fourier transformation of a contiguous sequence.
      template <typename Tag>
      void transform(select_sign<Tag> tag, real_type* seq_begin)
  // FUTURE: move out of class body
  {
    if (n_ < 2) return;
    real_type* C = seq_begin;
    real_type* CH = &(*(CH_.begin()));
    const real_type* WA = &(*(WA_.begin()));
    bool NA = false;
    std::size_t L1 = 1;
    std::size_t IW = 0;
    for (std::size_t K1 = 0; K1 < factors_.size(); K1++) {
      std::size_t IP = factors_[K1];
      std::size_t L2 = IP*L1;
      std::size_t IDO = n_/L2;
      std::size_t IDOT = IDO+IDO;
      std::size_t IDL1 = IDOT*L1;
      if (IP == 4) {
        std::size_t IX2 = IW+IDOT;
        std::size_t IX3 = IX2+IDOT;
        if (!NA) {
          pass4(tag, IDOT,L1,C,CH,WA+IW,WA+IX2,WA+IX3);
        }
        else {
          pass4(tag, IDOT,L1,CH,C,WA+IW,WA+IX2,WA+IX3);
        }
        NA = !NA;
      }
      else if (IP == 2) {
        if (!NA) {
          pass2(tag, IDOT,L1,C,CH,WA+IW);
        }
        else {
          pass2(tag, IDOT,L1,CH,C,WA+IW);
        }
        NA = !NA;
      }
      else if (IP == 3) {
        std::size_t IX2 = IW+IDOT;
        if (!NA) {
          pass3(tag, IDOT,L1,C,CH,WA+IW,WA+IX2);
        }
        else {
          pass3(tag, IDOT,L1,CH,C,WA+IW,WA+IX2);
        }
        NA = !NA;
      }
      else if (IP == 5) {
        std::size_t IX2 = IW+IDOT;
        std::size_t IX3 = IX2+IDOT;
        std::size_t IX4 = IX3+IDOT;
        if (!NA) {
          pass5(tag, IDOT,L1,C,CH,WA+IW,WA+IX2,WA+IX3,WA+IX4);
        }
        else {
          pass5(tag, IDOT,L1,CH,C,WA+IW,WA+IX2,WA+IX3,WA+IX4);
        }
        NA = !NA;
      }
      else {
        bool NAC;
        if (!NA) {
          passg(tag, NAC, IDOT,IP,L1,IDL1,IW,C,CH,WA);
        }
        else {
          passg(tag, NAC, IDOT,IP,L1,IDL1,IW,CH,C,WA);
        }
        if (NAC) NA = !NA;
      }
      L1 = L2;
      IW = IW+(IP-1)*IDOT;
    }
    if (NA) {
      std::copy(CH, CH + 2 * n_, seq_begin);
    }
  }
    private:
      // Constants.
      real_type two_pi_;
      real_type one_half_;
      real_type cos30_;
      real_type sin18_;
      real_type cos18_;
      real_type sin36_;
      real_type cos36_;

      static real_type deg_as_rad(const real_type& phi)
      {
        return phi * std::atan(real_type(1)) / real_type(45);
      }

      // Scratch space.
      af::shared<real_type> WA_;
      af::shared<real_type> CH_;

      // Codelets for prime factors 2,3,4,5 and a general transform.

      template <class Tag>
      void pass2(select_sign<Tag>,
                 std::size_t IDO,
                 std::size_t L1,
                 real_type* CC_begin,
                 real_type* CH_begin,
                 const real_type* WA1)
  // FUTURE: move out of class body
  {
    dim3 CC(CC_begin, IDO, 2, L1);
    dim3 CH(CH_begin, IDO, L1, 2);
    if (IDO == 2) {
      for (std::size_t K = 0; K < L1; K++) {
        CH(0,K,0) = CC(0,0,K) + CC(0,1,K);
        CH(1,K,0) = CC(1,0,K) + CC(1,1,K);
        CH(0,K,1) = CC(0,0,K) - CC(0,1,K);
        CH(1,K,1) = CC(1,0,K) - CC(1,1,K);
      }
    }
    else {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 0; I0 < IDO; I0 += 2) {
          std::size_t I1 = I0 + 1;
          CH(I0,K,0) = CC(I0,0,K) + CC(I0,1,K);
          CH(I1,K,0) = CC(I1,0,K) + CC(I1,1,K);
          real_type TR2 = CC(I0,0,K) - CC(I0,1,K);
          real_type TI2 = CC(I1,0,K) - CC(I1,1,K);
          CH(I0,K,1) = select_sign<Tag>::plusminus(WA1[I0]*TR2,WA1[I1]*TI2);
          CH(I1,K,1) = select_sign<Tag>::minusplus(WA1[I0]*TI2,WA1[I1]*TR2);
        }
      }
    }
  }

      template <class Tag>
      void pass3(select_sign<Tag>,
                 std::size_t IDO,
                 std::size_t L1,
                 real_type* CC_begin,
                 real_type* CH_begin,
                 const real_type* WA1,
                 const real_type* WA2)
  // FUTURE: move out of class body
  {
    dim3 CC(CC_begin, IDO, 3, L1);
    dim3 CH(CH_begin, IDO, L1, 3);
    real_type TAUI = select_sign<Tag>::unaryminusplus(cos30_);
    if (IDO == 2) {
      for (std::size_t K = 0; K < L1; K++) {
        real_type TR2 = CC(0,1,K) + CC(0,2,K);
        real_type TI2 = CC(1,1,K) + CC(1,2,K);
        real_type CR2 = CC(0,0,K) - one_half_ * TR2;
        real_type CI2 = CC(1,0,K) - one_half_ * TI2;
        real_type CR3 = TAUI * (CC(0,1,K) - CC(0,2,K));
        real_type CI3 = TAUI * (CC(1,1,K) - CC(1,2,K));
        CH(0,K,0) = CC(0,0,K) + TR2;
        CH(1,K,0) = CC(1,0,K) + TI2;
        CH(0,K,1) = CR2 - CI3;
        CH(1,K,1) = CI2 + CR3;
        CH(0,K,2) = CR2 + CI3;
        CH(1,K,2) = CI2 - CR3;
      }
    }
    else {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 0; I0 < IDO; I0 += 2) {
          std::size_t I1 = I0 + 1;
          real_type TR2 = CC(I0,1,K) + CC(I0,2,K);
          real_type TI2 = CC(I1,1,K) + CC(I1,2,K);
          real_type CR2 = CC(I0,0,K) - one_half_ * TR2;
          real_type CI2 = CC(I1,0,K) - one_half_ * TI2;
          real_type CR3 = TAUI * (CC(I0,1,K) - CC(I0,2,K));
          real_type CI3 = TAUI * (CC(I1,1,K) - CC(I1,2,K));
          real_type DR2 = CR2 - CI3;
          real_type DI2 = CI2 + CR3;
          real_type DR3 = CR2 + CI3;
          real_type DI3 = CI2 - CR3;
          CH(I0,K,0) = CC(I0,0,K) + TR2;
          CH(I1,K,0) = CC(I1,0,K) + TI2;
          CH(I0,K,1) = select_sign<Tag>::plusminus(WA1[I0]*DR2,WA1[I1]*DI2);
          CH(I1,K,1) = select_sign<Tag>::minusplus(WA1[I0]*DI2,WA1[I1]*DR2);
          CH(I0,K,2) = select_sign<Tag>::plusminus(WA2[I0]*DR3,WA2[I1]*DI3);
          CH(I1,K,2) = select_sign<Tag>::minusplus(WA2[I0]*DI3,WA2[I1]*DR3);
        }
      }
    }
  }

      template <class Tag>
      void pass4(select_sign<Tag>,
                 std::size_t IDO,
                 std::size_t L1,
                 real_type* CC_begin,
                 real_type* CH_begin,
                 const real_type* WA1,
                 const real_type* WA2,
                 const real_type* WA3)
  // FUTURE: move out of class body
  {
    dim3 CC(CC_begin, IDO, 4, L1);
    dim3 CH(CH_begin, IDO, L1, 4);
    if (IDO == 2) {
      for (std::size_t K = 0; K < L1; K++) {
        real_type TR1 = CC(0,0,K)-CC(0,2,K);
        real_type TI1 = CC(1,0,K)-CC(1,2,K);
        real_type TR2 = CC(0,0,K)+CC(0,2,K);
        real_type TI2 = CC(1,0,K)+CC(1,2,K);
        real_type TR3 = CC(0,1,K)+CC(0,3,K);
        real_type TI3 = CC(1,1,K)+CC(1,3,K);
        real_type
        TR4 = select_sign<Tag>::unaryplusminus(CC(1,1,K)-CC(1,3,K));
        real_type
        TI4 = select_sign<Tag>::unaryplusminus(CC(0,3,K)-CC(0,1,K));
        CH(0,K,0) = TR2+TR3;
        CH(1,K,0) = TI2+TI3;
        CH(0,K,1) = TR1+TR4;
        CH(1,K,1) = TI1+TI4;
        CH(0,K,2) = TR2-TR3;
        CH(1,K,2) = TI2-TI3;
        CH(0,K,3) = TR1-TR4;
        CH(1,K,3) = TI1-TI4;
      }
    }
    else {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 0; I0 < IDO; I0 += 2) {
          std::size_t I1 = I0 + 1;
          real_type TR1 = CC(I0,0,K)-CC(I0,2,K);
          real_type TI1 = CC(I1,0,K)-CC(I1,2,K);
          real_type TR2 = CC(I0,0,K)+CC(I0,2,K);
          real_type TI2 = CC(I1,0,K)+CC(I1,2,K);
          real_type TR3 = CC(I0,1,K)+CC(I0,3,K);
          real_type TI3 = CC(I1,1,K)+CC(I1,3,K);
          real_type
          TR4 = select_sign<Tag>::unaryplusminus(CC(I1,1,K)-CC(I1,3,K));
          real_type
          TI4 = select_sign<Tag>::unaryplusminus(CC(I0,3,K)-CC(I0,1,K));
          real_type CR3 = TR2-TR3;
          real_type CI3 = TI2-TI3;
          real_type CR2 = TR1+TR4;
          real_type CR4 = TR1-TR4;
          real_type CI2 = TI1+TI4;
          real_type CI4 = TI1-TI4;
          CH(I0,K,0) = TR2+TR3;
          CH(I1,K,0) = TI2+TI3;
          CH(I0,K,1) = select_sign<Tag>::plusminus(WA1[I0]*CR2,WA1[I1]*CI2);
          CH(I1,K,1) = select_sign<Tag>::minusplus(WA1[I0]*CI2,WA1[I1]*CR2);
          CH(I0,K,2) = select_sign<Tag>::plusminus(WA2[I0]*CR3,WA2[I1]*CI3);
          CH(I1,K,2) = select_sign<Tag>::minusplus(WA2[I0]*CI3,WA2[I1]*CR3);
          CH(I0,K,3) = select_sign<Tag>::plusminus(WA3[I0]*CR4,WA3[I1]*CI4);
          CH(I1,K,3) = select_sign<Tag>::minusplus(WA3[I0]*CI4,WA3[I1]*CR4);
        }
      }
    }
  }

      template <class Tag>
      void pass5(select_sign<Tag>,
                 std::size_t IDO,
                 std::size_t L1,
                 real_type* CC_begin,
                 real_type* CH_begin,
                 const real_type* WA1,
                 const real_type* WA2,
                 const real_type* WA3,
                 const real_type* WA4)
  // FUTURE: move out of class body
  {
    dim3 CC(CC_begin, IDO, 5, L1);
    dim3 CH(CH_begin, IDO, L1, 5);
    real_type TI11 = select_sign<Tag>::unaryminusplus(cos18_);
    real_type TI12 = select_sign<Tag>::unaryminusplus(sin36_);
    if (IDO == 2) {
      for (std::size_t K = 0; K < L1; K++) {
        real_type TR2 = CC(0,1,K)+CC(0,4,K);
        real_type TI2 = CC(1,1,K)+CC(1,4,K);
        real_type TR3 = CC(0,2,K)+CC(0,3,K);
        real_type TI3 = CC(1,2,K)+CC(1,3,K);
        real_type TR4 = CC(0,2,K)-CC(0,3,K);
        real_type TI4 = CC(1,2,K)-CC(1,3,K);
        real_type TR5 = CC(0,1,K)-CC(0,4,K);
        real_type TI5 = CC(1,1,K)-CC(1,4,K);
        real_type CR2 = CC(0,0,K)+sin18_*TR2-cos36_*TR3;
        real_type CI2 = CC(1,0,K)+sin18_*TI2-cos36_*TI3;
        real_type CR3 = CC(0,0,K)-cos36_*TR2+sin18_*TR3;
        real_type CI3 = CC(1,0,K)-cos36_*TI2+sin18_*TI3;
        real_type CR5 = TI11*TR5+TI12*TR4;
        real_type CI5 = TI11*TI5+TI12*TI4;
        real_type CR4 = TI12*TR5-TI11*TR4;
        real_type CI4 = TI12*TI5-TI11*TI4;
        CH(0,K,0) = CC(0,0,K)+TR2+TR3;
        CH(1,K,0) = CC(1,0,K)+TI2+TI3;
        CH(0,K,1) = CR2-CI5;
        CH(1,K,1) = CI2+CR5;
        CH(0,K,2) = CR3-CI4;
        CH(1,K,2) = CI3+CR4;
        CH(0,K,3) = CR3+CI4;
        CH(1,K,3) = CI3-CR4;
        CH(0,K,4) = CR2+CI5;
        CH(1,K,4) = CI2-CR5;
      }
    }
    else {
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I0 = 0; I0 < IDO; I0 += 2) {
          std::size_t I1 = I0 + 1;
          real_type TR2 = CC(I0,1,K)+CC(I0,4,K);
          real_type TI2 = CC(I1,1,K)+CC(I1,4,K);
          real_type TR3 = CC(I0,2,K)+CC(I0,3,K);
          real_type TI3 = CC(I1,2,K)+CC(I1,3,K);
          real_type TR4 = CC(I0,2,K)-CC(I0,3,K);
          real_type TI4 = CC(I1,2,K)-CC(I1,3,K);
          real_type TR5 = CC(I0,1,K)-CC(I0,4,K);
          real_type TI5 = CC(I1,1,K)-CC(I1,4,K);
          real_type CR2 = CC(I0,0,K)+sin18_*TR2-cos36_*TR3;
          real_type CI2 = CC(I1,0,K)+sin18_*TI2-cos36_*TI3;
          real_type CR3 = CC(I0,0,K)-cos36_*TR2+sin18_*TR3;
          real_type CI3 = CC(I1,0,K)-cos36_*TI2+sin18_*TI3;
          real_type CR4 = TI12*TR5-TI11*TR4;
          real_type CI4 = TI12*TI5-TI11*TI4;
          real_type CR5 = TI11*TR5+TI12*TR4;
          real_type CI5 = TI11*TI5+TI12*TI4;
          real_type DR3 = CR3-CI4;
          real_type DR4 = CR3+CI4;
          real_type DI3 = CI3+CR4;
          real_type DI4 = CI3-CR4;
          real_type DR5 = CR2+CI5;
          real_type DR2 = CR2-CI5;
          real_type DI5 = CI2-CR5;
          real_type DI2 = CI2+CR5;
          CH(I0,K,0) = CC(I0,0,K)+TR2+TR3;
          CH(I1,K,0) = CC(I1,0,K)+TI2+TI3;
          CH(I0,K,1) = select_sign<Tag>::plusminus(WA1[I0]*DR2,WA1[I1]*DI2);
          CH(I1,K,1) = select_sign<Tag>::minusplus(WA1[I0]*DI2,WA1[I1]*DR2);
          CH(I0,K,2) = select_sign<Tag>::plusminus(WA2[I0]*DR3,WA2[I1]*DI3);
          CH(I1,K,2) = select_sign<Tag>::minusplus(WA2[I0]*DI3,WA2[I1]*DR3);
          CH(I0,K,3) = select_sign<Tag>::plusminus(WA3[I0]*DR4,WA3[I1]*DI4);
          CH(I1,K,3) = select_sign<Tag>::minusplus(WA3[I0]*DI4,WA3[I1]*DR4);
          CH(I0,K,4) = select_sign<Tag>::plusminus(WA4[I0]*DR5,WA4[I1]*DI5);
          CH(I1,K,4) = select_sign<Tag>::minusplus(WA4[I0]*DI5,WA4[I1]*DR5);
        }
      }
    }
  }

      template <class Tag>
      void passg(select_sign<Tag>,
                 bool& NAC,
                 std::size_t IDO,
                 std::size_t IP,
                 std::size_t L1,
                 std::size_t IDL1,
                 std::size_t IW,
                 real_type* CC_begin,
                 real_type* CH_begin,
                 const real_type* WA)
  // FUTURE: move out of class body
  {
    dim3 CC(CC_begin, IDO, IP, L1);
    dim3 C1(CC_begin, IDO, L1, IP);
    dim2 C2(CC_begin, IDL1, IP);
    dim3 CH(CH_begin, IDO, L1, IP);
    dim2 CH2(CH_begin, IDL1, IP);
    std::size_t IDOT = IDO/2;
    std::size_t IPPH = (IP+1)/2;
    std::size_t IDP = IP*IDO;
    if (IDO >= L1) {
      for (std::size_t J = 1; J < IPPH; J++) {
        std::size_t JC = IP-J;
        for (std::size_t K = 0; K < L1; K++) {
          for (std::size_t I = 0; I < IDO; I++) {
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
          }
        }
      }
      for (std::size_t K = 0; K < L1; K++) {
        for (std::size_t I = 0; I < IDO; I++) {
          CH(I,K,0) = CC(I,0,K);
        }
      }
    }
    else {
      for (std::size_t J = 1; J < IPPH; J++) {
        std::size_t JC = IP-J;
        for (std::size_t I = 0; I < IDO; I++) {
          for (std::size_t K = 0; K < L1; K++) {
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
          }
        }
      }
      for (std::size_t I = 0; I < IDO; I++) {
        for (std::size_t K = 0; K < L1; K++) {
          CH(I,K,0) = CC(I,0,K);
        }
      }
    }
    std::size_t IDL = 0;
    for (std::size_t L = 1; L < IPPH; L++) {
      std::size_t LC = IP-L;
      for (std::size_t IK = 0; IK < IDL1; IK++) {
        C2(IK,L) = CH2(IK,0)+WA[IW+IDL]*CH2(IK,1);
        C2(IK,LC) = select_sign<Tag>::unaryminusplus(
                                         WA[IW+IDL+1]*CH2(IK,IP-1));
      }
      std::size_t IDLJ = IDL;
      IDL += IDO;
      for (std::size_t J = 2; J < IPPH; J++) {
        std::size_t JC = IP-J;
        IDLJ += IDL;
        if (IDLJ >= IDP) IDLJ = IDLJ-IDP;
        real_type WAR = WA[IW+IDLJ];
        real_type WAI = WA[IW+IDLJ+1];
        for (std::size_t IK = 0; IK < IDL1; IK++) {
          C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J);
          C2(IK,LC) = select_sign<Tag>::minusplus(C2(IK,LC),WAI*CH2(IK,JC));
        }
      }
    }
    std::size_t J;
    for (J = 1; J < IPPH; J++) {
      for (std::size_t IK = 0; IK < IDL1; IK++) {
         CH2(IK,0) += CH2(IK,J);
      }
    }
    for (J = 1; J < IPPH; J++) {
      std::size_t JC = IP-J;
      for (std::size_t IK0 = 0; IK0 < IDL1; IK0 += 2) {
        std::size_t IK1 = IK0 + 1;
        CH2(IK0,J) = C2(IK0,J)-C2(IK1,JC);
        CH2(IK0,JC) = C2(IK0,J)+C2(IK1,JC);
        CH2(IK1,J) = C2(IK1,J)+C2(IK0,JC);
        CH2(IK1,JC) = C2(IK1,J)-C2(IK0,JC);
      }
    }
    if (IDO == 2) {
      NAC = true;
      return;
    }
    std::copy(CH_begin, CH_begin + IDL1, CC_begin);
    for (J = 1; J < IP; J++) {
      for (std::size_t K = 0; K < L1; K++) {
        C1(0,K,J) = CH(0,K,J);
        C1(1,K,J) = CH(1,K,J);
      }
    }
    if (IDOT <= L1) {
      std::size_t IDIJ = 2;
      for (std::size_t J = 1; J < IP; J++) {
        for (std::size_t I0 = 2; I0 < IDO; I0 += 2) {
          std::size_t I1 = I0 + 1;
          for (std::size_t K = 0; K < L1; K++) {
            C1(I0,K,J) = select_sign<Tag>::plusminus(
              WA[IW+IDIJ]*CH(I0,K,J),WA[IW+IDIJ+1]*CH(I1,K,J));
            C1(I1,K,J) = select_sign<Tag>::minusplus(
              WA[IW+IDIJ]*CH(I1,K,J),WA[IW+IDIJ+1]*CH(I0,K,J));
          }
          IDIJ += 2;
        }
        IDIJ += 2;
      }
    }
    else {
      std::size_t IDJ = 2;
      for (std::size_t J = 1; J < IP; J++) {
        for (std::size_t K = 0; K < L1; K++) {
          std::size_t IDIJ = IDJ;
          for (std::size_t I0 = 2; I0 < IDO; I0 += 2) {
            std::size_t I1 = I0 + 1;
            C1(I0,K,J) = select_sign<Tag>::plusminus(
              WA[IW+IDIJ]*CH(I0,K,J),WA[IW+IDIJ+1]*CH(I1,K,J));
            C1(I1,K,J) = select_sign<Tag>::minusplus(
              WA[IW+IDIJ]*CH(I1,K,J),WA[IW+IDIJ+1]*CH(I0,K,J));
            IDIJ += 2;
          }
        }
        IDJ += IDO;
      }
    }
    NAC = false;
    return;
  }
  };

  template <typename RealType,
            typename ComplexType>
  complex_to_complex<RealType, ComplexType>::complex_to_complex(std::size_t n)
    : factorization(n, false), WA_(2 * n), CH_(2 * n)
  {
    if (n_ < 2) return;
    // Precompute constants for real_type.
    two_pi_ = real_type(8) * std::atan(real_type(1));
    one_half_ = real_type(1) / real_type(2);
    cos30_ = std::cos(deg_as_rad(30));
    sin18_ = std::sin(deg_as_rad(18));
    cos18_ = std::cos(deg_as_rad(18));
    sin36_ = std::sin(deg_as_rad(36));
    cos36_ = std::cos(deg_as_rad(36));
    // Computation of the sin and cos terms.
    // Based on the second part of fftpack41/cffti1.f.
    real_type* WA = &(*(WA_.begin()));
    real_type ARGH = two_pi_ / real_type(n_);
    std::size_t I = 0;
    std::size_t L1 = 1;
    for (std::size_t K1 = 0; K1 < factors_.size(); K1++) {
      std::size_t IP = factors_[K1];
      std::size_t LD = 0;
      std::size_t L2 = L1 * IP;
      std::size_t IDOT = 2 * (n_ / L2) + 2;
      for (std::size_t J = 0; J < IP - 1; J++) {
        std::size_t I1 = I;
        WA[I  ] = real_type(1);
        WA[I+1] = real_type(0);
        LD += L1;
        std::size_t FI = 0;
        real_type ARGLD = real_type(LD) * ARGH;
        for (std::size_t II = 4; II <= IDOT; II += 2) {
          I += 2;
          FI++;
          real_type ARG = real_type(FI) * ARGLD;
          WA[I  ] = std::cos(ARG);
          WA[I+1] = std::sin(ARG);
        }
        if (IP > 5) {
          WA[I1  ] = WA[I  ];
          WA[I1+1] = WA[I+1];
        }
      }
      L1 = L2;
    }
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H
