#ifndef SCITBX_VEC2_H
#define SCITBX_VEC2_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/operator_traits_builtin.h>

namespace scitbx {

  //! Two-dimensional vector.
  /*! This class can be used to represent points, vectors, normals.
      The usual vector operations are available.

      A 2D adaptation of the scitbx::vec3 class
      Future: Unify this with vec3 class using dimensionality
      as a template argument.
   */
  template <typename NumType>
  class vec2 : public af::tiny_plain<NumType, 2>
  {
    public:
      typedef typename af::tiny_plain<NumType, 2> base_type;

      //! Default constructor. Elements are not initialized.
      vec2() {}
      //! All elements are initialized with e.
      vec2(NumType const& e)
        : base_type(e, e)
      {}
      //! Constructor.
      vec2(NumType const& e0, NumType const& e1)
        : base_type(e0, e1)
      {}
      //! Constructor.
      vec2(base_type const& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      vec2(const NumType* a)
      {
        for(std::size_t i=0;i<2;i++) this->elems[i] = a[i];
      }

      //! Test if all elements are 0.
      bool is_zero() const
      {
        return !(this->elems[0] || this->elems[1]);
      }

      //! Minimum of vector elements.
      NumType
      min() const
      {
        NumType result = this->elems[0];
        if (result > this->elems[1]) result = this->elems[1];
        return result;
      }

      //! Maximum of vector elements.
      NumType
      max() const
      {
        NumType result = this->elems[0];
        if (result < this->elems[1]) result = this->elems[1];
        return result;
      }

      //! In-place minimum of this and other for each element.
      void
      each_update_min(vec2 const& other)
      {
        for(std::size_t i=0;i<2;i++) {
          if (this->elems[i] > other[i]) this->elems[i] = other[i];
        }
      }

      //! In-place maximum of this and other for each element.
      void
      each_update_max(vec2 const& other)
      {
        for(std::size_t i=0;i<2;i++) {
          if (this->elems[i] < other[i]) this->elems[i] = other[i];
        }
      }

      //! Sum of vector elements.
      NumType
      sum() const
      {
        return this->elems[0] + this->elems[1];
      }

      //! Product of vector elements.
      NumType
      product() const
      {
        return this->elems[0] * this->elems[1];
      }

      //! Return the length squared of the vector.
      NumType length_sq() const
      {
        return (*this) * (*this);
      }

      //! Return the length of the vector.
      NumType length() const
      {
        return NumType(std::sqrt(length_sq()));
      }

      //! Return normalized vector.
      vec2 normalize() const
      {
        return (*this) / length();
      }

      //! Return angle (in radians) between this and other.
      NumType angle(vec2 const& other) const
      {
        return NumType(
          std::acos(((*this) * other) / (length() * other.length())));
      }

      //! Return the reflection vector.
      /*! @param n is the surface normal which has to be of unit length.
       */
      vec2 reflect(vec2 const& n) const
      {
        return (*this) - NumType(2) * ((*this) * n) * n;
      }

      //! Return the transmitted vector.
      /*! @param n is the surface normal which has to be of unit length.
          @param eta is the relative index of refraction. If the returned
          vector is zero then there is no transmitted light because
          of total internal reflection.
       */
      vec2 refract(vec2 const& n, NumType const& eta) const
      {
        NumType dot = (*this) * n;
        NumType k = 1. - eta*eta*(1. - dot*dot);
        if (k < NumType(0)) return vec2(0,0);
        return eta * (*this) - (eta * dot + std::sqrt(k)) * n;
      }

      //! Returns an orthogonal vector.
      /*! Returns a vector that is orthogonal to this.
       */
      vec2 ortho() const
      {
        NumType x = abs_(this->elems[0]);
        NumType y = abs_(this->elems[1]);
        return vec2(-this->elems[1], this->elems[0]);
      }

      //! Copies the vector to a tiny instance.
      af::tiny<NumType, 2>
      as_tiny() const { return af::tiny<NumType, 2>(*this); }

    private:
      static NumType abs_(NumType const& x)
      {
        if (x < NumType(0)) return -x;
        return x;
      }
  };

  //! Test equality.
  template <typename NumType>
  inline
  bool
  operator==(
    vec2<NumType> const& lhs,
    vec2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType>
  inline
  bool
  operator==(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      if (lhs[i] != rhs   ) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of rhs == lhs.
  template <typename NumType>
  inline
  bool
  operator==(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      if (lhs    != rhs[i]) return false;
    }
    return true;
  }

  //! Test inequality.
  template <typename NumType>
  inline
  bool
  operator!=(
    vec2<NumType> const& lhs,
    vec2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of rhs != lhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec2<NumType>
  operator+(
    vec2<NumType> const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec2<NumType>
  operator+(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec2<NumType>
  operator+(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec2<NumType>
  operator-(
    vec2<NumType> const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec2<NumType>
  operator-(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec2<NumType>
  operator-(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Dot product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  typename af::binary_operator_traits<NumTypeLhs, NumTypeRhs>::arithmetic
  operator*(
    vec2<NumTypeLhs> const& lhs,
    vec2<NumTypeRhs> const& rhs)
  {
    typename af::binary_operator_traits<NumTypeLhs, NumTypeRhs>::arithmetic
    result(0);
    for(std::size_t i=0;i<2;i++) {
      result += lhs[i] * rhs[i];
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  vec2<NumType>
  operator*(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  vec2<NumType>
  operator*(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec2<NumType>
  operator/(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec2<NumType>
  operator/(
    vec2<NumType> const& lhs,
    std::size_t const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec2<NumType>
  operator/(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  vec2<NumType>
  operator%(
    vec2<NumType> const& lhs,
    NumType const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  vec2<NumType>
  operator%(
    NumType const& lhs,
    vec2<NumType> const& rhs)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator+=(
    vec2<NumType>& lhs,
    vec2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator+=(
    vec2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator-=(
    vec2<NumType>& lhs,
    vec2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator-=(
    vec2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place multiplication.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator*=(
    vec2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place division.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator/=(
    vec2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place modulus operation.
  template <typename NumType>
  inline
  vec2<NumType>&
  operator%=(
    vec2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<2;i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Element-wise unary minus.
  template <typename NumType>
  inline
  vec2<NumType>
  operator-(
    vec2<NumType> const& v)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  vec2<NumType>
  operator+(
    vec2<NumType> const& v)
  {
    vec2<NumType> result;
    for(std::size_t i=0;i<2;i++) {
      result[i] = +v[i];
    }
    return result;
  }

  //! Return the length of the vector.
  template <typename NumType>
  inline
  NumType
  abs(
    vec2<NumType> const& v)
  {
    return v.length();
  }

} // namespace scitbx

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::vec2<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_VEC2_H
