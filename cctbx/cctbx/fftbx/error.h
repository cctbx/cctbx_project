// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 15: Created, based on cctbx/error.h (rwgk)
 */

#ifndef CCTBX_FFTBX_ERROR_H
#define CCTBX_FFTBX_ERROR_H

#include <exception>
#include <string>

namespace cctbx {

//! Fast Fourier Transform toolbox.
/*! The Fast Fourier Transform toolbox provides
    1-dimensional and 3-dimensional in-place
    complex-to-complex and real-to-complex transforms.
    There are no restrictions on the lengths of the
    sequences to be transformed.  However, the transforms
    are most efficient when the sequence length is a
    product of 2, 3 and 5.
    <p>
    The fftbx is a pure C++ adaption of the
    complex-to-complex and real-to-complex transforms of
    FFTPACK Version 4.1 (Nov 1988) by Paul Swarztrauber
    (see http://www.scd.ucar.edu/softlib/FFTPACK.html
    and the references below).
    <p>
    The core of the fftbx is implemented exclusively in C++
    header files. Therefore it is not necessary to link
    with a support library. Only the optional Python
    bindings that are provided for convenience need to be
    compiled.
    <p>
    Currently, the fftbx is templated on the vector type
    that contains the sequence to be transformed. At the
    point of this writing this has only been tested with
    std::vector<double>, but should work with other
    vector types that implement a similar interface.
    <p>
    The vector type is expected to contain instances
    of a floating-point type, such as float or double.
    Currently, vectors of std::complex<> are not
    directly supported. This will be resolved as
    the fftbx matures.
    <p>
    A related problem is that the in-place real-to-complex
    forward transform expects as input an array of
    real numbers and transforms this to an array of
    complex numbers (or vice versa for the backward
    transform). It is not yet clear how this could
    be handled in a type-safe way.
    <p>
    The fftbx is reasonably fast on the machines where it
    was tested (Tru64 Unix, Linux, Windows with Visual C++
    6), but is not as fast as FFTW (http://www.fftw.org/).
    Depending on the machine, compiler settings, transform
    lengths and the dimensionality (1D or 3D), the fftbx is
    about 30%-300% slower than FFTW for transform lengths
    that are a multiple of 2, 3 and 5. This is no surprise,
    since FFTW is a very sophisticated library that adjusts
    itself to a particular hardware, and was developed
    specifically for modern, chached-based architectures.
    In particular, the performance on Athlon Linux machines is
    very impressive. Unfortunately, the FFTW license does
    not permit the inclusion in an Open Source project
    (the GNU General Public License is too restrictive).
    Therefore the public domain FFTPACK was chosen as
    the basis for the fftbx.
    <p>
    It is hoped that the fftbx facilitates the development
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
namespace fftbx {

  //! Exeption type used by fftbx.
  class error : public std::exception {
    public:
      //! Constructor.
      error(const std::string& msg) throw() {
        message = std::string("fftbx Error: ") + msg;
      }
      //! Virtual destructor.
      virtual ~error() throw() {}
      //! Access to the error messages.
      virtual const char* what() const throw() {
        return message.c_str();
      }
    protected:
      //! The held message for the exception.
      std::string message;
  };

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_ERROR_H
