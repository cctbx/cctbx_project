#ifndef CCTBX_SGTBX_ROT_MX_H
#define CCTBX_SGTBX_ROT_MX_H

#include <cctbx/sgtbx/tr_vec.h>
#include <boost/rational.hpp>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_rot_mx(const char* file, long line);

  class rot_mx_info; // forward declaration

  //! 3x3 rotation matrix.
  /*! The elements of the matrix are stored as integers and a common
      denominator. The actual value of an element is obtained by
      dividing the integer number by the denominator.
   */
  class rot_mx
  {
    public:
      //! Initialization of a diagonal matrix with the given denominator.
      /*! The diagonal elements are defined as diagonal * denominator.
       */
      explicit
      rot_mx(int denominator = 1, int diagonal = 1)
      : num_(diagonal * denominator), den_(denominator)
      {}

      //! Initialization with numerator m and the given denominator.
      explicit
      rot_mx(sg_mat3 const& m, int denominator = 1)
      : num_(m), den_(denominator)
      {}

      //! Initialization with the given elements and denominator.
      rot_mx(int m00, int m01, int m02,
             int m10, int m11, int m12,
             int m20, int m21, int m22,
             int denominator = 1)
        : num_(m00,m01,m02, m10,m11,m12, m20,m21, m22),
          den_(denominator)
      {}

      //! Numerator of the rotation matrix.
      sg_mat3 const& num() const { return num_; }
      //! Numerator of the rotation matrix.
      sg_mat3&       num()       { return num_; }

      //! i'th element of the numerator of the rotation matrix.
      int const& operator[](std::size_t i) const { return num_[i]; }
      //! i'th element of the numerator of the rotation matrix.
      int&       operator[](std::size_t i)       { return num_[i]; }

      //! (r*3+c)'th element of the numerator of the rotation matrix.
      int const& operator()(int r, int c) const { return num_(r, c); }
      //! (r*3+c)'th element of the numerator of the rotation matrix.
      int&       operator()(int r, int c)       { return num_(r, c); }

      //! Denominator of the rotation matrix.
      int const& den() const { return den_; }
      //! Denominator of the rotation matrix.
      int&       den()       { return den_; }

      //! True only if both the numerators and the denominators are equal.
      bool operator==(rot_mx const& rhs) const
      {
        if (den_ != rhs.den_) return false;
        return num_.const_ref().all_eq(rhs.num_.const_ref());
      }

      //! False only if both the numerators and the denominators are equal.
      bool operator!=(rot_mx const& rhs) const
      {
        return !((*this) == rhs);
      }

      //! True only if den() != 0.
      bool is_valid() const { return den_ != 0; }

      /*! \brief True only if this is a diagonal matrix and all diagonal
          elements are equal to the denominator.
       */
      bool is_unit_mx() const { return num_ == sg_mat3(den_); }

      //! This minus the unit matrix.
      rot_mx minus_unit_mx() const
      {
        rot_mx result(*this);
        for (std::size_t i=0;i<9;i+=4) result[i] -= den_;
        return result;
      }

      //! New rotation matrix with denominator new_den.
      /*! An exception is thrown if the old rotation matrix
          cannot be represented using the new denominator.
       */
      rot_mx new_denominator(int new_den) const;

      //! New rotation matrix with num()*factor and den()*factor.
      rot_mx scale(int factor) const
      {
        if (factor == 1) return *this;
        return rot_mx(num_ * factor, den_ * factor);
      }

      //! Determinant as rational number.
      /*! An exception is thrown if den() <= 0.
       */
      boost::rational<int>
      determinant() const
      {
        CCTBX_ASSERT(den_ > 0);
        return boost::rational<int>(num_.determinant(), den_*den_*den_);
      }

      //! Inverse of this matrix.
      /*! An exception is thrown if the result
          cannot be represented using the denominator den().
       */
      rot_mx inverse() const;

      //! New rotation matrix with -num(), den().
      rot_mx operator-() const { return rot_mx(-num_, den_); }

      //! Addition of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      friend rot_mx operator+(rot_mx const& lhs, rot_mx const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return rot_mx(lhs.num_ + rhs.num_, lhs.den_);
      }

      //! Subtraction of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      friend rot_mx operator-(rot_mx const& lhs, rot_mx const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return rot_mx(lhs.num_ - rhs.num_, lhs.den_);
      }

      //! In-place addition of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      rot_mx& operator+=(rot_mx const& rhs)
      {
        CCTBX_ASSERT(den_ == rhs.den_);
        num_ += rhs.num_;
        return *this;
      }

      //! Matrix multiplication.
      /*! The denominator of the result is the product lhs.den() * rhs.den().
       */
      friend rot_mx operator*(rot_mx const& lhs, rot_mx const& rhs)
      {
        return rot_mx(lhs.num_ * rhs.num_, lhs.den_ * rhs.den_);
      }

      //! Matrix*vector multiplication.
      /*! The denominator of the result is the product lhs.den() * rhs.den().
       */
      friend tr_vec operator*(rot_mx const& lhs, tr_vec const& rhs)
      {
        return tr_vec(lhs.num_ * rhs.num(), lhs.den_ * rhs.den());
      }

      //! Vector*matrix multiplication.
      /*! The denominator of the result is the product lhs.den() * rhs.den().
       */
      friend tr_vec operator*(tr_vec const& lhs, rot_mx const& rhs)
      {
        return tr_vec(lhs.num() * rhs.num_, lhs.den() * rhs.den_);
      }

      //! Matrix*vector multiplication, numerator only.
      friend sg_vec3 operator*(rot_mx const& lhs, sg_vec3 const& rhs)
      {
        return sg_vec3(lhs.num_ * rhs);
      }

      //! New rotation matrix with num()*rhs, den().
      friend rot_mx operator*(rot_mx const& lhs, int rhs)
      {
        return rot_mx(lhs.num_ * rhs, lhs.den_);
      }

      //! New rotation matrix with rhs.num()*lhs, rhs.den().
      friend rot_mx operator*(int lhs, rot_mx const& rhs)
      {
        return rhs * lhs;
      }

      //! In-place multiplication of numerator with rhs.
      rot_mx& operator*=(int rhs)
      {
        num_ *= rhs;
        return *this;
      }

      //! New rotation matrix with num() / rhs, den().
      /*! An exception is thrown if the result cannot be represented
          using den().
       */
      friend rot_mx operator/(rot_mx const& lhs, int rhs);

      //! Cancellation of factors.
      /*! The denominator of the new rotation matrix is made
          as small as possible.
       */
      rot_mx cancel() const;

      //! Matrix inversion with cancellation of factors.
      rot_mx inverse_cancel() const;

      //! Matrix multiplication with cancellation of factors.
      rot_mx multiply(rot_mx const& rhs) const
      {
        return ((*this) * rhs).cancel();
      }

      //! Matrix*vector multiplication with cancellation of factors.
      tr_vec multiply(tr_vec const& rhs) const
      {
        return ((*this) * rhs).cancel();
      }

      //! num()/rhs,den() with cancellation of factors.
      rot_mx divide(int rhs) const;

      //! Rotation-part type (1, 2, 3, 4, 6, -1, -2=m, -3, -4, -6)
      /*! See also: cctbx::sgtbx::rot_mx_info
       */
      int type() const;

      //! Rotational order.
      /*! Number of times the matrix must be multiplied with itself
          in order to obtain the unit matrix.

          See also: cctbx::sgtbx::rot_mx_info
       */
      int order(int type=0) const;

      //! Sum of repeated products of this matrix with itself.
      /*! this + this*this + this*this*this + ... + this**order()

          Restriction: the numerator must be one.
       */
      rot_mx accumulate(int type=0) const;

      //! Convenience method for constructing a rot_mx_info instance.
      rot_mx_info info() const;

      //! Conversion to a floating-point array.
      template <typename FloatType>
      scitbx::mat3<FloatType>
      as_floating_point(scitbx::type_holder<FloatType>) const
      {
        return scitbx::mat3<FloatType>(num_) / FloatType(den_);
      }

      //! Conversion to an array with element type double.
      /*! Python: float()
       */
      scitbx::mat3<double>
      as_double() const
      {
        return as_floating_point(scitbx::type_holder<double>());
      }

    private:
      sg_mat3 num_;
      int den_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_ROT_MX_H
