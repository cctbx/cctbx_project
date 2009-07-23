#ifndef SCITBX_MATH_IMAGINARY_UNIT_H
#define SCITBX_MATH_IMAGINARY_UNIT_H

#include <complex>
#include <boost/operators.hpp>

namespace scitbx { namespace math {

  /// The complex number i
  struct imaginary_unit_t {};

  template <typename T>
  inline
  std::complex<T> operator*(imaginary_unit_t, std::complex<T> const &z) {
    return std::complex<T>(-z.imag(), z.real());
  }

  template <typename T>
  inline
  std::complex<T> operator*(std::complex<T> const &z, imaginary_unit_t) {
    return std::complex<T>(-z.imag(), z.real());
  }

  /// An imaginary number
  /** This allows natural notations to be used resulting in more performant code
      too. For example,

      @code

      imaginary_unit_t i;
      std::complex<int> z = 1 + 2*i;
      z *= 3*i; // z == (-6, 3)
                // and only 2 multiplications were performed (and no addition)
      @endcode

      This is the scheme that W. Kahan has advocated for years, which has
      eventually been adopted in C99, and which under discussion for C++
      (c.f. http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2004/n1612.pdf ).

      All meaningful operations involving imaginary_unit_t, imaginary and
      std::complex are implemented, except -imaginary_unit_t.
  */
  template <typename T>
  class imaginary
    : boost::partially_ordered< imaginary<T>
    , boost::addable< imaginary<T>
    , boost::subtractable< imaginary<T>
    , boost::dividable2< imaginary<T>, T
    , boost::multipliable2< imaginary<T>, T
    > > > > >
  {
    T b_;
    typedef std::complex<T> c_t;
    typedef T re_t;
    typedef imaginary im_t;

    public:

    /// Classic construction
    imaginary(re_t b=0) : b_(b) {}

    /// Multiplier of i
    re_t b() {
      return b_;
    }

    /** @name Group2
        group for +
    */
    //@{
    im_t operator+=(im_t z) { b_ += z.b_; return *this; }
    im_t operator-=(im_t z) { b_ -= z.b_; return *this; }
    im_t operator-() { return im_t(-b_); }
    //@}

    /** @name Group3
        vector space over real numbers
    */
    //@{
    im_t operator*=(re_t s) { b_ *= s; return *this; }
    im_t operator/=(re_t s) { b_ /= s; return *this; }
    //@}

    /** @name Group4
        partial order
    */
    //@{
    bool operator<(im_t z) { return b_ < z.b_; }
    bool operator==(im_t z) { return b_ == z.b_; }
    //@}

  };

  /** @defgroup Group1 Construction through multiplication
                   of real and imaginary_unit_t
      @{
  */
  template <typename T>
  imaginary<T> operator*(imaginary_unit_t, T b) {
    return imaginary<T>(b);
  }

  template <typename T>
  imaginary<T> operator*(T b, imaginary_unit_t) {
    return imaginary<T>(b);
  }
  /** @} */

  /// real / imaginary
  template <typename T>
  imaginary<T> operator/(T a, imaginary<T> z) {
    return imaginary<T>(-a/z.b());
  }

  /** @defgroup Group2 Construction of std::complex with natural notation
      @{
  */
  template <typename T>
  std::complex<T> operator+(T a, imaginary<T> b) {
    return std::complex<T>(a, b.b());
  }

  template <typename T>
  std::complex<T> operator+(imaginary<T> b, T a) {
    return std::complex<T>(a, b.b());
  }
  /** @} */

  /** @defgroup Group3 Mixed operations with std::complex
      @{
  */
  template <typename T>
  std::complex<T> operator+(std::complex<T> const &z, imaginary<T> y) {
    return std::complex<T>(z.real(), z.imag() + y.b());
  }

  template <typename T>
  std::complex<T> operator+(imaginary<T> y, std::complex<T> const &z) {
    return std::complex<T>(z.real(), z.imag() + y.b());
  }

  template <typename T>
  std::complex<T> operator-(std::complex<T> const &z, imaginary<T> y) {
    return std::complex<T>(z.real(), z.imag() - y.b());
  }

  template <typename T>
  std::complex<T> operator-(imaginary<T> y, std::complex<T> const &z) {
    return std::complex<T>(-z.real(), y.b() - z.imag());
  }

  template <typename T>
  std::complex<T> operator*(std::complex<T> const &z, imaginary<T> y) {
    return std::complex<T>(-y.b()*z.imag(), y.b()*z.real());
  }

  template <typename T>
  std::complex<T> operator*(imaginary<T> y, std::complex<T> const &z) {
    return std::complex<T>(-y.b()*z.imag(), y.b()*z.real());
  }

  template <typename T>
  std::complex<T> operator/(std::complex<T> const &z, imaginary<T> y) {
    return std::complex<T>(z.imag()/y.b(), -z.real()/y.b());
  }

  template <typename T>
  std::complex<T> operator/(imaginary<T> y, std::complex<T> const &z) {
    T d2 = std::norm(z);
    return std::complex<T>(y.b()*z.imag()/d2, y.b()*z.real()/d2);
  }
  //@}

}}

#endif // GUARD
