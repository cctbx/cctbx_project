#ifndef CCTBX_SGTBX_TR_VEC_H
#define CCTBX_SGTBX_TR_VEC_H

#include <cctbx/sgtbx/basic.h>
#include <scitbx/math/modulo.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/type_holder.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_tr_vec(const char* file, long line);

  //! Translation vector.
  /*! The elements of the vector are stored as integers and a common
      denominator. The actual value of an element is obtained by
      dividing the integer number by the denominator.
   */
  class tr_vec
  {
    public:
      /*! \brief Initialization of a 0,0,0-translation with the denominator
          tr_den.
       */
      explicit
      tr_vec(int tr_den=sg_t_den)
      :
        num_(0,0,0), den_(tr_den)
      {}

      /*! \brief Initialization a translation with the integer numerator v
          and the denominator tr_den.
       */
      explicit
      tr_vec(sg_vec3 const& v, int tr_den=sg_t_den)
      :
        num_(v), den_(tr_den)
      {}

      /*! \brief Initialization of a v0,v1,v2-translation with the
          denominator tr_den.
       */
      tr_vec(int v0, int v1, int v2, int tr_den=sg_t_den)
      :
        num_(v0,v1,v2), den_(tr_den)
      {}

      //! Numerator of the translation vector.
      sg_vec3 const&
      num() const { return num_; }
      //! Numerator of the translation vector.
      sg_vec3&
      num()       { return num_; }

      //! i'th element of the numerator of the translation vector.
      int const&
      operator[](std::size_t i) const { return num_[i]; }
      //! i'th element of the numerator of the translation vector.
      int&
      operator[](std::size_t i)       { return num_[i]; }

      //! Denominator of the translation vector.
      int const&
      den() const { return den_; }
      //! Denominator of the translation vector.
      int&
      den()       { return den_; }

      //! True only if both the numerators and the denominators are equal.
      bool
      operator==(tr_vec const& rhs) const
      {
        if (den_ != rhs.den_) return false;
        return num_.const_ref().all_eq(rhs.num_.const_ref());
      }

      //! False only if both the numerators and the denominators are equal.
      bool
      operator!=(tr_vec const& rhs) const
      {
        return !((*this) == rhs);
      }

      //! True only if den() != 0.
      bool
      is_valid() const { return den_ != 0; }

      //! True only if all elements of the numertor are equal to zero.
      bool
      is_zero() const { return num_.is_zero(); }

      //! New translation vector with denominator new_den.
      /*! An exception is thrown if the old translation vector
          cannot be represented using the new denominator.
       */
      tr_vec
      new_denominator(int new_den) const;

      //! New translation vector with num()*factor and den()*factor.
      tr_vec
      scale(int factor) const
      {
        if (factor == 1) return *this;
        return tr_vec(num_ * factor, den_ * factor);
      }

      //! New translation vector with 0 <= num()[i] < den().
      tr_vec
      mod_positive() const
      {
        tr_vec result(num_, den_);
        for(std::size_t i=0;i<3;i++) {
          result.num_[i] = scitbx::math::mod_positive(result.num_[i], den_);
        }
        return result;
      }

      //! New translation vector with -den()/2 < num()[i] <= den()/2.
      tr_vec
      mod_short() const
      {
        tr_vec result(num_, den_);
        for(std::size_t i=0;i<3;i++) {
          result.num_[i] = scitbx::math::mod_short(result.num_[i], den_);
        }
        return result;
      }

      //! New translation vector with -num(), den().
      tr_vec
      operator-() const { return tr_vec(-num_, den_); }

      //! Addition of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      friend tr_vec
      operator+(tr_vec const& lhs, tr_vec const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return tr_vec(lhs.num_ + rhs.num_, lhs.den_);
      }

      //! Subtraction of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      friend tr_vec
      operator-(tr_vec const& lhs, tr_vec const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return tr_vec(lhs.num_ - rhs.num_, lhs.den_);
      }

      //! In-place addition of numerators.
      /*! An exception is thrown if the denominators are not equal.
       */
      tr_vec&
      operator+=(tr_vec const& rhs)
      {
        CCTBX_ASSERT(den_ == rhs.den_);
        num_ += rhs.num_;
        return *this;
      }

      //! New translation vector with num() * rhs, den().
      friend tr_vec
      operator*(tr_vec const& lhs, int rhs)
      {
        return tr_vec(lhs.num_ * rhs, lhs.den_);
      }

      //! New translation vector with num() / rhs, den().
      /*! An exception is thrown if the result cannot be represented
          using den().
       */
      friend tr_vec
      operator/(tr_vec const& lhs, int rhs);

      //! New translation vector with rhs.num() * lhs, rhs.den().
      friend tr_vec
      operator*(int const& lhs, tr_vec const& rhs)
      {
        return rhs * lhs;
      }

      //! Cancellation of factors.
      /*! The denominator of the new translation vector is made
          as small as possible.
       */
      tr_vec cancel() const;

      //! Addition with cancellation of factors.
      tr_vec
      plus(tr_vec const& rhs) const;

      //! Subtraction with cancellation of factors.
      tr_vec
      minus(tr_vec const& rhs) const;

      //! Conversion to a floating-point array.
      template <typename FloatType>
      scitbx::vec3<FloatType>
      as_floating_point(scitbx::type_holder<FloatType>) const
      {
        return scitbx::vec3<FloatType>(num_) / FloatType(den_);
      }

      //! Conversion to an array with element type double.
      scitbx::vec3<double>
      as_double() const
      {
        return as_floating_point(scitbx::type_holder<double>());
      }

      //! Conversion to a string with rational or decimal numbers.
      /*! Python: str()
       */
      std::string
      as_string(bool decimal=false, const char* separator=",") const;

    private:
      sg_vec3 num_;
      int den_;
  };

  //! Constructor for initialization of constants.
  struct tr_vec_12 : tr_vec
  {
    tr_vec_12(int v0, int v1, int v2)
    : tr_vec(tr_vec(v0,v1,v2) * (sg_t_den / 12))
    {}
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_TR_VEC_H
