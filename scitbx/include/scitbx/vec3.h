/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy of cctbx/vec3.h (Ralf W. Grosse-Kunstleve)
     2002 May: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_VEC3_H
#define SCITBX_VEC3_H

#include <scitbx/array_family/tiny_plain.h>
#include <scitbx/array_family/operator_traits_builtin.h>

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
      vec3 normalize() const
      {
        return (*this) / length();
      }

      //! Return angle (in radians) between this and other.
      NumType angle(vec3 const& other) const
      {
        return NumType(
          std::acos(((*this) * other) / (length() * other.length())));
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

      //! Returns an orthogonal vector.
      /*! Returns a vector that is orthogonal to this.
       */
      vec3 ortho() const
      {
        NumType x = abs_(this->elems[0]);
        NumType y = abs_(this->elems[1]);
        NumType z = abs_(this->elems[2]);
        // Is z the smallest element? Then use x and y.
        if (z <= x && z <= y) return vec3(-this->elems[1], this->elems[0], 0);
        // Is y smallest element? Then use x and z.
        if (y <= x && y <= z) return vec3(-this->elems[2], 0, this->elems[0]);
        // x is smallest.
        return vec3(0, -this->elems[2], this->elems[1]);
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

  //! Element-wise inplace addition.
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

  //! Element-wise inplace addition.
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

  //! Element-wise inplace difference.
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

  //! Element-wise inplace difference.
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

  //! Element-wise inplace multiplication.
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

  //! Element-wise inplace division.
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

  //! Element-wise inplace modulus operation.
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

} // namespace scitbx

#endif // SCITBX_VEC3_H
