#ifndef SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H
#define SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H

#include <scitbx/fftpack/factorization.h>
#include <scitbx/fftpack/detail/ref.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/auto_array.h>
#include <algorithm>
#include <complex>
#include <cmath>

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
        return wa_;
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
      void
      forward(
        ComplexOrRealIterOrPtrType seq_begin,
        real_type* scratch=0)
      {
        transform(select_sign<forward_tag>(), &(*seq_begin), scratch);
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
      void
      backward(
        ComplexOrRealIterOrPtrType seq_begin,
        real_type* scratch=0)
      {
        transform(select_sign<backward_tag>(), &(*seq_begin), scratch);
      }

      //! Generic in-place Fourier transformation of a contiguous sequence.
      template <typename Tag>
      void
      transform(
        select_sign<Tag> tag,
        complex_type* seq_begin,
        real_type* scratch=0)
      {
        transform(tag, reinterpret_cast<real_type*>(seq_begin), scratch);
      }

      //! Generic in-place Fourier transformation of a contiguous sequence.
      template <typename Tag>
      void
      transform(
        select_sign<Tag> tag,
        real_type* c, /* seq_begin */
        real_type* ch=0 /* scratch */)
  // FUTURE: move out of class body
  {
    if (n_ < 2) return;
    scitbx::auto_array<real_type> scratch;
    if (ch == 0) {
      scratch = auto_array<real_type>(new real_type[n_ * 2]);
      ch = scratch.get();
    }
    const real_type* wa = &(*(wa_.begin()));
    bool na = false;
    std::size_t l1 = 1;
    std::size_t iw = 0;
    for (std::size_t k1 = 0; k1 < factors_.size(); k1++) {
      std::size_t ip = factors_[k1];
      std::size_t l2 = ip*l1;
      std::size_t ido = n_/l2;
      std::size_t idot = ido+ido;
      std::size_t idl1 = idot*l1;
      if (ip == 4) {
        std::size_t ix2 = iw+idot;
        std::size_t ix3 = ix2+idot;
        if (!na) {
          pass4(tag, idot,l1,c,ch,wa+iw,wa+ix2,wa+ix3);
        }
        else {
          pass4(tag, idot,l1,ch,c,wa+iw,wa+ix2,wa+ix3);
        }
        na = !na;
      }
      else if (ip == 2) {
        if (!na) {
          pass2(tag, idot,l1,c,ch,wa+iw);
        }
        else {
          pass2(tag, idot,l1,ch,c,wa+iw);
        }
        na = !na;
      }
      else if (ip == 3) {
        std::size_t ix2 = iw+idot;
        if (!na) {
          pass3(tag, idot,l1,c,ch,wa+iw,wa+ix2);
        }
        else {
          pass3(tag, idot,l1,ch,c,wa+iw,wa+ix2);
        }
        na = !na;
      }
      else if (ip == 5) {
        std::size_t ix2 = iw+idot;
        std::size_t ix3 = ix2+idot;
        std::size_t ix4 = ix3+idot;
        if (!na) {
          pass5(tag, idot,l1,c,ch,wa+iw,wa+ix2,wa+ix3,wa+ix4);
        }
        else {
          pass5(tag, idot,l1,ch,c,wa+iw,wa+ix2,wa+ix3,wa+ix4);
        }
        na = !na;
      }
      else {
        bool nac;
        if (!na) {
          passg(tag, nac, idot,ip,l1,idl1,iw,c,ch,wa);
        }
        else {
          passg(tag, nac, idot,ip,l1,idl1,iw,ch,c,wa);
        }
        if (nac) na = !na;
      }
      l1 = l2;
      iw = iw+(ip-1)*idot;
    }
    if (na) {
      std::copy(ch, ch + 2 * n_, c);
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
      af::shared<real_type> wa_;

      // Codelets for prime factors 2,3,4,5 and a general transform.

      template <class Tag>
      void pass2(select_sign<Tag>,
                 std::size_t ido,
                 std::size_t l1,
                 real_type* cc_begin,
                 real_type* ch_begin,
                 const real_type* wa1)
  // FUTURE: move out of class body
  {
    dim3 cc(cc_begin, ido, 2, l1);
    dim3 ch(ch_begin, ido, l1, 2);
    if (ido == 2) {
      for (std::size_t k = 0; k < l1; k++) {
        ch(0,k,0) = cc(0,0,k) + cc(0,1,k);
        ch(1,k,0) = cc(1,0,k) + cc(1,1,k);
        ch(0,k,1) = cc(0,0,k) - cc(0,1,k);
        ch(1,k,1) = cc(1,0,k) - cc(1,1,k);
      }
    }
    else {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 0; i0 < ido; i0 += 2) {
          std::size_t i1 = i0 + 1;
          ch(i0,k,0) = cc(i0,0,k) + cc(i0,1,k);
          ch(i1,k,0) = cc(i1,0,k) + cc(i1,1,k);
          real_type tr2 = cc(i0,0,k) - cc(i0,1,k);
          real_type ti2 = cc(i1,0,k) - cc(i1,1,k);
          ch(i0,k,1) = select_sign<Tag>::plusminus(wa1[i0]*tr2,wa1[i1]*ti2);
          ch(i1,k,1) = select_sign<Tag>::minusplus(wa1[i0]*ti2,wa1[i1]*tr2);
        }
      }
    }
  }

      template <class Tag>
      void pass3(select_sign<Tag>,
                 std::size_t ido,
                 std::size_t l1,
                 real_type* cc_begin,
                 real_type* ch_begin,
                 const real_type* wa1,
                 const real_type* wa2)
  // FUTURE: move out of class body
  {
    dim3 cc(cc_begin, ido, 3, l1);
    dim3 ch(ch_begin, ido, l1, 3);
    real_type taui = select_sign<Tag>::unaryminusplus(cos30_);
    if (ido == 2) {
      for (std::size_t k = 0; k < l1; k++) {
        real_type tr2 = cc(0,1,k) + cc(0,2,k);
        real_type ti2 = cc(1,1,k) + cc(1,2,k);
        real_type cr2 = cc(0,0,k) - one_half_ * tr2;
        real_type ci2 = cc(1,0,k) - one_half_ * ti2;
        real_type cr3 = taui * (cc(0,1,k) - cc(0,2,k));
        real_type ci3 = taui * (cc(1,1,k) - cc(1,2,k));
        ch(0,k,0) = cc(0,0,k) + tr2;
        ch(1,k,0) = cc(1,0,k) + ti2;
        ch(0,k,1) = cr2 - ci3;
        ch(1,k,1) = ci2 + cr3;
        ch(0,k,2) = cr2 + ci3;
        ch(1,k,2) = ci2 - cr3;
      }
    }
    else {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 0; i0 < ido; i0 += 2) {
          std::size_t i1 = i0 + 1;
          real_type tr2 = cc(i0,1,k) + cc(i0,2,k);
          real_type ti2 = cc(i1,1,k) + cc(i1,2,k);
          real_type cr2 = cc(i0,0,k) - one_half_ * tr2;
          real_type ci2 = cc(i1,0,k) - one_half_ * ti2;
          real_type cr3 = taui * (cc(i0,1,k) - cc(i0,2,k));
          real_type ci3 = taui * (cc(i1,1,k) - cc(i1,2,k));
          real_type dr2 = cr2 - ci3;
          real_type di2 = ci2 + cr3;
          real_type dr3 = cr2 + ci3;
          real_type di3 = ci2 - cr3;
          ch(i0,k,0) = cc(i0,0,k) + tr2;
          ch(i1,k,0) = cc(i1,0,k) + ti2;
          ch(i0,k,1) = select_sign<Tag>::plusminus(wa1[i0]*dr2,wa1[i1]*di2);
          ch(i1,k,1) = select_sign<Tag>::minusplus(wa1[i0]*di2,wa1[i1]*dr2);
          ch(i0,k,2) = select_sign<Tag>::plusminus(wa2[i0]*dr3,wa2[i1]*di3);
          ch(i1,k,2) = select_sign<Tag>::minusplus(wa2[i0]*di3,wa2[i1]*dr3);
        }
      }
    }
  }

      template <class Tag>
      void pass4(select_sign<Tag>,
                 std::size_t ido,
                 std::size_t l1,
                 real_type* cc_begin,
                 real_type* ch_begin,
                 const real_type* wa1,
                 const real_type* wa2,
                 const real_type* wa3)
  // FUTURE: move out of class body
  {
    dim3 cc(cc_begin, ido, 4, l1);
    dim3 ch(ch_begin, ido, l1, 4);
    if (ido == 2) {
      for (std::size_t k = 0; k < l1; k++) {
        real_type tr1 = cc(0,0,k)-cc(0,2,k);
        real_type ti1 = cc(1,0,k)-cc(1,2,k);
        real_type tr2 = cc(0,0,k)+cc(0,2,k);
        real_type ti2 = cc(1,0,k)+cc(1,2,k);
        real_type tr3 = cc(0,1,k)+cc(0,3,k);
        real_type ti3 = cc(1,1,k)+cc(1,3,k);
        real_type
        tr4 = select_sign<Tag>::unaryplusminus(cc(1,1,k)-cc(1,3,k));
        real_type
        ti4 = select_sign<Tag>::unaryplusminus(cc(0,3,k)-cc(0,1,k));
        ch(0,k,0) = tr2+tr3;
        ch(1,k,0) = ti2+ti3;
        ch(0,k,1) = tr1+tr4;
        ch(1,k,1) = ti1+ti4;
        ch(0,k,2) = tr2-tr3;
        ch(1,k,2) = ti2-ti3;
        ch(0,k,3) = tr1-tr4;
        ch(1,k,3) = ti1-ti4;
      }
    }
    else {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 0; i0 < ido; i0 += 2) {
          std::size_t i1 = i0 + 1;
          real_type tr1 = cc(i0,0,k)-cc(i0,2,k);
          real_type ti1 = cc(i1,0,k)-cc(i1,2,k);
          real_type tr2 = cc(i0,0,k)+cc(i0,2,k);
          real_type ti2 = cc(i1,0,k)+cc(i1,2,k);
          real_type tr3 = cc(i0,1,k)+cc(i0,3,k);
          real_type ti3 = cc(i1,1,k)+cc(i1,3,k);
          real_type
          tr4 = select_sign<Tag>::unaryplusminus(cc(i1,1,k)-cc(i1,3,k));
          real_type
          ti4 = select_sign<Tag>::unaryplusminus(cc(i0,3,k)-cc(i0,1,k));
          real_type cr3 = tr2-tr3;
          real_type ci3 = ti2-ti3;
          real_type cr2 = tr1+tr4;
          real_type cr4 = tr1-tr4;
          real_type ci2 = ti1+ti4;
          real_type ci4 = ti1-ti4;
          ch(i0,k,0) = tr2+tr3;
          ch(i1,k,0) = ti2+ti3;
          ch(i0,k,1) = select_sign<Tag>::plusminus(wa1[i0]*cr2,wa1[i1]*ci2);
          ch(i1,k,1) = select_sign<Tag>::minusplus(wa1[i0]*ci2,wa1[i1]*cr2);
          ch(i0,k,2) = select_sign<Tag>::plusminus(wa2[i0]*cr3,wa2[i1]*ci3);
          ch(i1,k,2) = select_sign<Tag>::minusplus(wa2[i0]*ci3,wa2[i1]*cr3);
          ch(i0,k,3) = select_sign<Tag>::plusminus(wa3[i0]*cr4,wa3[i1]*ci4);
          ch(i1,k,3) = select_sign<Tag>::minusplus(wa3[i0]*ci4,wa3[i1]*cr4);
        }
      }
    }
  }

      template <class Tag>
      void pass5(select_sign<Tag>,
                 std::size_t ido,
                 std::size_t l1,
                 real_type* cc_begin,
                 real_type* ch_begin,
                 const real_type* wa1,
                 const real_type* wa2,
                 const real_type* wa3,
                 const real_type* wa4)
  // FUTURE: move out of class body
  {
    dim3 cc(cc_begin, ido, 5, l1);
    dim3 ch(ch_begin, ido, l1, 5);
    real_type ti11 = select_sign<Tag>::unaryminusplus(cos18_);
    real_type ti12 = select_sign<Tag>::unaryminusplus(sin36_);
    if (ido == 2) {
      for (std::size_t k = 0; k < l1; k++) {
        real_type tr2 = cc(0,1,k)+cc(0,4,k);
        real_type ti2 = cc(1,1,k)+cc(1,4,k);
        real_type tr3 = cc(0,2,k)+cc(0,3,k);
        real_type ti3 = cc(1,2,k)+cc(1,3,k);
        real_type tr4 = cc(0,2,k)-cc(0,3,k);
        real_type ti4 = cc(1,2,k)-cc(1,3,k);
        real_type tr5 = cc(0,1,k)-cc(0,4,k);
        real_type ti5 = cc(1,1,k)-cc(1,4,k);
        real_type cr2 = cc(0,0,k)+sin18_*tr2-cos36_*tr3;
        real_type ci2 = cc(1,0,k)+sin18_*ti2-cos36_*ti3;
        real_type cr3 = cc(0,0,k)-cos36_*tr2+sin18_*tr3;
        real_type ci3 = cc(1,0,k)-cos36_*ti2+sin18_*ti3;
        real_type cr5 = ti11*tr5+ti12*tr4;
        real_type ci5 = ti11*ti5+ti12*ti4;
        real_type cr4 = ti12*tr5-ti11*tr4;
        real_type ci4 = ti12*ti5-ti11*ti4;
        ch(0,k,0) = cc(0,0,k)+tr2+tr3;
        ch(1,k,0) = cc(1,0,k)+ti2+ti3;
        ch(0,k,1) = cr2-ci5;
        ch(1,k,1) = ci2+cr5;
        ch(0,k,2) = cr3-ci4;
        ch(1,k,2) = ci3+cr4;
        ch(0,k,3) = cr3+ci4;
        ch(1,k,3) = ci3-cr4;
        ch(0,k,4) = cr2+ci5;
        ch(1,k,4) = ci2-cr5;
      }
    }
    else {
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i0 = 0; i0 < ido; i0 += 2) {
          std::size_t i1 = i0 + 1;
          real_type tr2 = cc(i0,1,k)+cc(i0,4,k);
          real_type ti2 = cc(i1,1,k)+cc(i1,4,k);
          real_type tr3 = cc(i0,2,k)+cc(i0,3,k);
          real_type ti3 = cc(i1,2,k)+cc(i1,3,k);
          real_type tr4 = cc(i0,2,k)-cc(i0,3,k);
          real_type ti4 = cc(i1,2,k)-cc(i1,3,k);
          real_type tr5 = cc(i0,1,k)-cc(i0,4,k);
          real_type ti5 = cc(i1,1,k)-cc(i1,4,k);
          real_type cr2 = cc(i0,0,k)+sin18_*tr2-cos36_*tr3;
          real_type ci2 = cc(i1,0,k)+sin18_*ti2-cos36_*ti3;
          real_type cr3 = cc(i0,0,k)-cos36_*tr2+sin18_*tr3;
          real_type ci3 = cc(i1,0,k)-cos36_*ti2+sin18_*ti3;
          real_type cr4 = ti12*tr5-ti11*tr4;
          real_type ci4 = ti12*ti5-ti11*ti4;
          real_type cr5 = ti11*tr5+ti12*tr4;
          real_type ci5 = ti11*ti5+ti12*ti4;
          real_type dr3 = cr3-ci4;
          real_type dr4 = cr3+ci4;
          real_type di3 = ci3+cr4;
          real_type di4 = ci3-cr4;
          real_type dr5 = cr2+ci5;
          real_type dr2 = cr2-ci5;
          real_type di5 = ci2-cr5;
          real_type di2 = ci2+cr5;
          ch(i0,k,0) = cc(i0,0,k)+tr2+tr3;
          ch(i1,k,0) = cc(i1,0,k)+ti2+ti3;
          ch(i0,k,1) = select_sign<Tag>::plusminus(wa1[i0]*dr2,wa1[i1]*di2);
          ch(i1,k,1) = select_sign<Tag>::minusplus(wa1[i0]*di2,wa1[i1]*dr2);
          ch(i0,k,2) = select_sign<Tag>::plusminus(wa2[i0]*dr3,wa2[i1]*di3);
          ch(i1,k,2) = select_sign<Tag>::minusplus(wa2[i0]*di3,wa2[i1]*dr3);
          ch(i0,k,3) = select_sign<Tag>::plusminus(wa3[i0]*dr4,wa3[i1]*di4);
          ch(i1,k,3) = select_sign<Tag>::minusplus(wa3[i0]*di4,wa3[i1]*dr4);
          ch(i0,k,4) = select_sign<Tag>::plusminus(wa4[i0]*dr5,wa4[i1]*di5);
          ch(i1,k,4) = select_sign<Tag>::minusplus(wa4[i0]*di5,wa4[i1]*dr5);
        }
      }
    }
  }

      template <class Tag>
      void passg(select_sign<Tag>,
                 bool& nac,
                 std::size_t ido,
                 std::size_t ip,
                 std::size_t l1,
                 std::size_t idl1,
                 std::size_t iw,
                 real_type* cc_begin,
                 real_type* ch_begin,
                 const real_type* wa)
  // FUTURE: move out of class body
  {
    dim3 cc(cc_begin, ido, ip, l1);
    dim3 c1(cc_begin, ido, l1, ip);
    dim2 c2(cc_begin, idl1, ip);
    dim3 ch(ch_begin, ido, l1, ip);
    dim2 ch2(ch_begin, idl1, ip);
    std::size_t idot = ido/2;
    std::size_t ipph = (ip+1)/2;
    std::size_t idp = ip*ido;
    if (ido >= l1) {
      for (std::size_t j = 1; j < ipph; j++) {
        std::size_t jc = ip-j;
        for (std::size_t k = 0; k < l1; k++) {
          for (std::size_t i = 0; i < ido; i++) {
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k);
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k);
          }
        }
      }
      for (std::size_t k = 0; k < l1; k++) {
        for (std::size_t i = 0; i < ido; i++) {
          ch(i,k,0) = cc(i,0,k);
        }
      }
    }
    else {
      for (std::size_t j = 1; j < ipph; j++) {
        std::size_t jc = ip-j;
        for (std::size_t i = 0; i < ido; i++) {
          for (std::size_t k = 0; k < l1; k++) {
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k);
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k);
          }
        }
      }
      for (std::size_t i = 0; i < ido; i++) {
        for (std::size_t k = 0; k < l1; k++) {
          ch(i,k,0) = cc(i,0,k);
        }
      }
    }
    std::size_t idl = 0;
    for (std::size_t l = 1; l < ipph; l++) {
      std::size_t lc = ip-l;
      for (std::size_t ik = 0; ik < idl1; ik++) {
        c2(ik,l) = ch2(ik,0)+wa[iw+idl]*ch2(ik,1);
        c2(ik,lc) = select_sign<Tag>::unaryminusplus(
                                         wa[iw+idl+1]*ch2(ik,ip-1));
      }
      std::size_t idlj = idl;
      idl += ido;
      for (std::size_t j = 2; j < ipph; j++) {
        std::size_t jc = ip-j;
        idlj += idl;
        if (idlj >= idp) idlj = idlj-idp;
        real_type war = wa[iw+idlj];
        real_type wai = wa[iw+idlj+1];
        for (std::size_t ik = 0; ik < idl1; ik++) {
          c2(ik,l) = c2(ik,l)+war*ch2(ik,j);
          c2(ik,lc) = select_sign<Tag>::minusplus(c2(ik,lc),wai*ch2(ik,jc));
        }
      }
    }
    std::size_t j;
    for (j = 1; j < ipph; j++) {
      for (std::size_t ik = 0; ik < idl1; ik++) {
         ch2(ik,0) += ch2(ik,j);
      }
    }
    for (j = 1; j < ipph; j++) {
      std::size_t jc = ip-j;
      for (std::size_t ik0 = 0; ik0 < idl1; ik0 += 2) {
        std::size_t ik1 = ik0 + 1;
        ch2(ik0,j) = c2(ik0,j)-c2(ik1,jc);
        ch2(ik0,jc) = c2(ik0,j)+c2(ik1,jc);
        ch2(ik1,j) = c2(ik1,j)+c2(ik0,jc);
        ch2(ik1,jc) = c2(ik1,j)-c2(ik0,jc);
      }
    }
    if (ido == 2) {
      nac = true;
      return;
    }
    std::copy(ch_begin, ch_begin + idl1, cc_begin);
    for (j = 1; j < ip; j++) {
      for (std::size_t k = 0; k < l1; k++) {
        c1(0,k,j) = ch(0,k,j);
        c1(1,k,j) = ch(1,k,j);
      }
    }
    if (idot <= l1) {
      std::size_t idij = 2;
      for (std::size_t j = 1; j < ip; j++) {
        for (std::size_t i0 = 2; i0 < ido; i0 += 2) {
          std::size_t i1 = i0 + 1;
          for (std::size_t k = 0; k < l1; k++) {
            c1(i0,k,j) = select_sign<Tag>::plusminus(
              wa[iw+idij]*ch(i0,k,j),wa[iw+idij+1]*ch(i1,k,j));
            c1(i1,k,j) = select_sign<Tag>::minusplus(
              wa[iw+idij]*ch(i1,k,j),wa[iw+idij+1]*ch(i0,k,j));
          }
          idij += 2;
        }
        idij += 2;
      }
    }
    else {
      std::size_t idj = 2;
      for (std::size_t j = 1; j < ip; j++) {
        for (std::size_t k = 0; k < l1; k++) {
          std::size_t idij = idj;
          for (std::size_t i0 = 2; i0 < ido; i0 += 2) {
            std::size_t i1 = i0 + 1;
            c1(i0,k,j) = select_sign<Tag>::plusminus(
              wa[iw+idij]*ch(i0,k,j),wa[iw+idij+1]*ch(i1,k,j));
            c1(i1,k,j) = select_sign<Tag>::minusplus(
              wa[iw+idij]*ch(i1,k,j),wa[iw+idij+1]*ch(i0,k,j));
            idij += 2;
          }
        }
        idj += ido;
      }
    }
    nac = false;
    return;
  }
  };

  template <typename RealType,
            typename ComplexType>
  complex_to_complex<RealType, ComplexType>::complex_to_complex(std::size_t n)
    : factorization(n, false), wa_(2 * n)
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
    real_type* wa = &(*(wa_.begin()));
    real_type argh = two_pi_ / real_type(n_);
    std::size_t i = 0;
    std::size_t l1 = 1;
    for (std::size_t k1 = 0; k1 < factors_.size(); k1++) {
      std::size_t ip = factors_[k1];
      std::size_t ld = 0;
      std::size_t l2 = l1 * ip;
      std::size_t idot = 2 * (n_ / l2) + 2;
      for (std::size_t j = 0; j < ip - 1; j++) {
        std::size_t i1 = i;
        wa[i  ] = real_type(1);
        wa[i+1] = real_type(0);
        ld += l1;
        std::size_t fi = 0;
        real_type argld = real_type(ld) * argh;
        for (std::size_t ii = 4; ii <= idot; ii += 2) {
          i += 2;
          fi++;
          real_type arg = real_type(fi) * argld;
          wa[i  ] = std::cos(arg);
          wa[i+1] = std::sin(arg);
        }
        if (ip > 5) {
          wa[i1  ] = wa[i  ];
          wa[i1+1] = wa[i+1];
        }
      }
      l1 = l2;
    }
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_H
