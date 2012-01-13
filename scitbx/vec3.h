#ifndef SCITBX_VEC3_H
#define SCITBX_VEC3_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/operator_traits_builtin.h>
#include <boost/optional.hpp>

namespace scitbx {

  //! Three-dimensional vector.
  /*! This class can be used to represent points, vectors, normals
      or even colors. The usual vector operations are available.

      C++ version of the python vec3 class by
      Matthias Baas (baas@ira.uka.de). See
      http://cgkit.sourceforge.net/
      for more information.
   */
  template <typename NumType>
  class vec3 : public af::tiny_plain<NumType, 3>
  {
    public:
      typedef typename af::tiny_plain<NumType, 3> base_type;

      //! Default constructor. Elements are not initialized.
      vec3() {}
      //! All elements are initialized with e.
      vec3(NumType const& e)
        : base_type(e, e, e)
      {}
      //! Constructor.
      vec3(NumType const& e0, NumType const& e1, NumType const& e2)
        : base_type(e0, e1, e2)
      {}
      //! Constructor.
      vec3(base_type const& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      vec3(const NumType* a)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = a[i];
      }

      //! Test if all elements are 0.
      bool is_zero() const
      {
        return !(this->elems[0] || this->elems[1] || this->elems[2]);
      }

      //! Minimum of vector elements.
      NumType
      min() const
      {
        NumType result = this->elems[0];
        if (result > this->elems[1]) result = this->elems[1];
        if (result > this->elems[2]) result = this->elems[2];
        return result;
      }

      //! Maximum of vector elements.
      NumType
      max() const
      {
        NumType result = this->elems[0];
        if (result < this->elems[1]) result = this->elems[1];
        if (result < this->elems[2]) result = this->elems[2];
        return result;
      }

      //! In-place minimum of this and other for each element.
      void
      each_update_min(vec3 const& other)
      {
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[i] > other[i]) this->elems[i] = other[i];
        }
      }

      //! In-place maximum of this and other for each element.
      void
      each_update_max(vec3 const& other)
      {
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[i] < other[i]) this->elems[i] = other[i];
        }
      }

      //! Sum of vector elements.
      NumType
      sum() const
      {
        return this->elems[0] + this->elems[1] + this->elems[2];
      }

      //! Product of vector elements.
      NumType
      product() const
      {
        return this->elems[0] * this->elems[1] * this->elems[2];
      }

      //! Cross product.
      vec3 cross(vec3 const& other) const
      {
        vec3 const& a = *this;
        vec3 const& b = other;
        return vec3(
          a[1] * b[2] - b[1] * a[2],
          a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]);
      }

      //! Element-wise multiplication.
      template <typename NumTypeRhs>
      vec3<
        typename af::binary_operator_traits<NumType, NumTypeRhs>::arithmetic>
      each_mul(
        vec3<NumTypeRhs> const& rhs) const
      {
        vec3<
          typename af::binary_operator_traits<NumType, NumTypeRhs>::arithmetic>
            result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = (*this)[i] * rhs[i];
        }
        return result;
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
      /*! length() == 0 will lead to a division-by-zero exception.
          Explicit checks are avoided to maximize performance.
          A maximally robust (but slower) calculation of the
          ||.||_2 norm is implemented in
          scitbx::math::accumulator::norm_accumulator .
       */
      vec3 normalize() const
      {
        return (*this) / length();
      }

      //! Return angle (in radians) between this and other.
      /*! This implementation is unsafe as it does not check
          if both vector lengths are different from zero.
       */
      NumType
      angle(
        vec3 const& other) const
      {
        return NumType(
          std::acos(((*this) * other) / (length() * other.length())));
      }

      //! Return angle (in radians) between this and other.
      /*! Safe implementation guarding against division-by-zero
          and rounding errors.
       */
      boost::optional<NumType>
      angle_rad(
        vec3 const& other) const
      {
        NumType den = length() * other.length();
        if (den == 0) return boost::optional<NumType>();
        NumType c = ((*this) * other) / den;
        if      (c < -1) c = -1;
        else if (c >  1) c =  1;
        return boost::optional<NumType>(NumType(std::acos(c)));
      }

      //! Return the reflection vector.
      /*! @param n is the surface normal which has to be of unit length.
       */
      vec3 reflect(vec3 const& n) const
      {
        return (*this) - NumType(2) * ((*this) * n) * n;
      }

      //! Return the transmitted vector.
      /*! @param n is the surface normal which has to be of unit length.
          @param eta is the relative index of refraction. If the returned
          vector is zero then there is no transmitted light because
          of total internal reflection.
       */
      vec3 refract(vec3 const& n, NumType const& eta) const
      {
        NumType dot = (*this) * n;
        NumType k = 1. - eta*eta*(1. - dot*dot);
        if (k < NumType(0)) return vec3(0,0,0);
        return eta * (*this) - (eta * dot + std::sqrt(k)) * n;
      }

      //! Returns a vector that is orthogonal to this.
      /*! If normalize is true, given the null vector (within
          floating-point precision) the result is (0,0,1).
       */
      vec3 ortho(bool normalize=false) const;

      //! Copies the vector to a tiny instance.
      af::tiny<NumType, 3>
      as_tiny() const { return af::tiny<NumType, 3>(*this); }

      //! Rotate the vector.
      /*! Use the rotation formula to rotate the vector around the direction
          vector through an angle given in radians.  The direction vector is
          normalized so need not be given as a vector of unit length.
       */
      vec3
      rotate_around_origin(vec3 const& direction, NumType const& angle) const {
        vec3 unit = direction.normalize();
        return unit_rotate_around_origin(unit,angle);
      }

      //! Rotate the vector.
      /*! Use the rotation formula to rotate the vector around the unit
          direction through an angle given in radians.  The caller
          guarantees that the direction vector is of unit length.
       */
      vec3
      unit_rotate_around_origin(vec3 const& unit, NumType const& angle) const {
        NumType cosang = std::cos(angle);
        return (*this)*cosang +
               unit*(unit*(*this))*(1.0-cosang)+
               (unit.cross(*this))*std::sin(angle);
      }

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
    vec3<NumType> const& lhs,
    vec3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType>
  inline
  bool
  operator==(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
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
    vec3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      if (lhs    != rhs[i]) return false;
    }
    return true;
  }

  //! Test inequality.
  template <typename NumType>
  inline
  bool
  operator!=(
    vec3<NumType> const& lhs,
    vec3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    vec3<NumType> const& lhs,
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
    vec3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec3<NumType>
  operator+(
    vec3<NumType> const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec3<NumType>
  operator+(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  vec3<NumType>
  operator+(
    NumType const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec3<NumType>
  operator-(
    vec3<NumType> const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec3<NumType>
  operator-(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  vec3<NumType>
  operator-(
    NumType const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Dot product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  typename af::binary_operator_traits<NumTypeLhs, NumTypeRhs>::arithmetic
  operator*(
    vec3<NumTypeLhs> const& lhs,
    vec3<NumTypeRhs> const& rhs)
  {
    typename af::binary_operator_traits<NumTypeLhs, NumTypeRhs>::arithmetic
    result(0);
    for(std::size_t i=0;i<3;i++) {
      result += lhs[i] * rhs[i];
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  vec3<NumType>
  operator*(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  vec3<NumType>
  operator*(
    NumType const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec3<NumType>
  operator/(
    vec3<NumType> const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] / rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec3<NumType>
  operator/(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec3<NumType>
  operator/(
    vec3<NumType> const& lhs,
    std::size_t const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  vec3<NumType>
  operator/(
    NumType const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  vec3<NumType>
  operator%(
    vec3<NumType> const& lhs,
    NumType const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  vec3<NumType>
  operator%(
    NumType const& lhs,
    vec3<NumType> const& rhs)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator+=(
    vec3<NumType>& lhs,
    vec3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator+=(
    vec3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator-=(
    vec3<NumType>& lhs,
    vec3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator-=(
    vec3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place multiplication.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator*=(
    vec3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place division.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator/=(
    vec3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place modulus operation.
  template <typename NumType>
  inline
  vec3<NumType>&
  operator%=(
    vec3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Element-wise unary minus.
  template <typename NumType>
  inline
  vec3<NumType>
  operator-(
    vec3<NumType> const& v)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  vec3<NumType>
  operator+(
    vec3<NumType> const& v)
  {
    vec3<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = +v[i];
    }
    return result;
  }

  //! Return the length of the vector.
  template <typename NumType>
  inline
  NumType
  abs(
    vec3<NumType> const& v)
  {
    return v.length();
  }

  template <typename NumType>
  vec3<NumType>
  vec3<NumType>::ortho(bool normalize) const
  {
    const NumType* e = this->elems;
    NumType x = abs_(e[0]);
    NumType y = abs_(e[1]);
    NumType z = abs_(e[2]);
    NumType u, v, w;
    // Is z the smallest element? Then use x and y.
    if (z <= x && z <= y) {
      u = e[1]; v = -e[0]; w = 0;
      if (normalize) {
        NumType d = std::sqrt(x*x + y*y);
        if (d != 0) { u /= d; v /= d; }
        else { u = 0; v = 0; w = 1; }
      }
    }
    // Is y smallest element? Then use x and z.
    else if (y <= x && y <= z) {
      u = -e[2]; v = 0; w = e[0];
      if (normalize) {
        NumType d = std::sqrt(x*x + z*z);
        if (d != 0) { u /= d; w /= d; }
        else { u = 0; v = 0; w = 1; }
      }
    }
    // x is smallest.
    else {
      u = 0; v = e[2]; w = -e[1];
      if (normalize) {
        NumType d = std::sqrt(y*y + z*z);
        if (d != 0) { v /= d; w /= d; }
        else { u = 0; v = 0; w = 1; }
      }
    }
    return vec3<NumType>(u,v,w);
  }

} // namespace scitbx

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::vec3<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_VEC3_H
