// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Created 2002 May (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAT3_H
#define CCTBX_MAT3_H

#include <utility>
#include <cctbx/error.h>
#include <cctbx/vec3.h>
#include <cctbx/array_family/tiny_reductions.h>

namespace cctbx {

  //! Matrix class (3x3).
  /*! This class represents a 3x3 matrix that can be used to store
      linear transformations.

      Enhanced version of the python mat3 class by
      Matthias Baas (baas@ira.uka.de). See
      http://cgkit.sourceforge.net/
      for more information.
   */
  template <typename NumType>
  class mat3 : public af::tiny_plain<NumType, 9>
  {
    public:
      typedef typename af::tiny_plain<NumType, 9> base_type;

      //! Default constructor. Elements are not initialized.
      mat3() {}
      //! Constructor.
      mat3(const NumType& e00, const NumType& e01, const NumType& e02,
           const NumType& e10, const NumType& e11, const NumType& e12,
           const NumType& e20, const NumType& e21, const NumType& e22)
        : base_type(e00, e01, e02, e10, e11, e12, e20, e21, e22)
      {}
      //! Constructor.
      mat3(const base_type& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      mat3(const NumType* a)
      {
        for(std::size_t i=0;i<9;i++) this->elems[i] = a[i];
      }
      //! Constructor for diagonal matrix.
      explicit
      mat3(const NumType& diag)
        : base_type(diag,0,0,0,diag,0,0,0,diag)
      {}
      //! Constructor for diagonal matrix.
      explicit
      mat3(const af::tiny_plain<NumType,3>& diag)
        : base_type(diag[0],0,0,0,diag[1],0,0,0,diag[2])
      {}
      //! Constructor for a rotation matrix.
      /*! @param angle is in radians.
       */
      mat3(const NumType& angle, const af::tiny_plain<NumType,3>& axis);

      //! Access elements with 2-dimensional indices.
      const NumType&
      operator()(std::size_t r, std::size_t c) const
      {
        return this->elems[r * 3 + c];
      }
      //! Access elements with 2-dimensional indices.
            NumType&
      operator()(std::size_t r, std::size_t c)
      {
        return this->elems[r * 3 + c];
      }

      //! Return a row.
      vec3<NumType>
      get_row(std::size_t i) const
      {
        return vec3<NumType>(&this->elems[i * 3]);
      }

      //! Set a row.
      void
      set_row(std::size_t i, const af::tiny_plain<NumType,3>& v)
      {
        std::copy(v.begin(), v.end(), &this->elems[i * 3]);
      }

      //! Swap two rows in place.
      void
      swap_rows(std::size_t i1, std::size_t i2) {
        std::swap_ranges(&(*this)(i1,0), &(*this)(i1+1,0), &(*this)(i2,0));
      }

      //! Return a column.
      vec3<NumType>
      get_column(std::size_t i) const
      {
        vec3<NumType> result;
        for(std::size_t j=0;j<3;j++) result[j] = this->elems[j * 3 + i];
        return result;
      }

      //! Set a column.
      void
      set_column(std::size_t i, const af::tiny_plain<NumType,3>& v)
      {
        for(std::size_t j=0;j<3;j++) this->elems[j * 3 + i] = v[j];
      }

      //! Swap two columns in place.
      void
      swap_columns(std::size_t i1, std::size_t i2) {
        for(std::size_t i=0;i<9;i+=3) {
          std::swap(this->elems[i + i1], this->elems[i + i2]);
        }
      }

      //! Return diagonal elements.
      vec3<NumType>
      diagonal() const
      {
        return MatrixLite::DiagonalElements(*this);
      }

      //! Return the transposed matrix.
      mat3
      transpose() const
      {
        mat3 result;
        MatrixLite::transpose(this->begin(), 3, 3, result.begin());
        return result;
      }

      //! Return trace (sum of diagonal elements).
      NumType
      trace() const
      {
        return af::sum(diagonal());
      }

      //! Return determinant.
      NumType
      determinant() const
      {
        return MatrixLite::Determinant(*this);
      }

      //! Return the transposed of the co-factor matrix.
      /*! The inverse matrix is obtained by dividing the result
          by the determinant().
       */
      mat3
      co_factor_matrix_transposed() const
      {
        return MatrixLite::CoFactorMxTp(*this);
      }

      //! Return the inverse matrix.
      /*! An exception is thrown if the matrix is not invertible,
          i.e. if the determinant() is zero.
       */
      mat3
      inverse() const
      {
        NumType d = determinant();
        if (d == NumType(0)) throw error("Matrix is not invertible.");
        return mat3(MatrixLite::CoFactorMxTp(*this)) / d;
      }

      //! Scale matrix in place.
      /*! Each row of this is multiplied element-wise with v.
       */
      mat3&
      scale(const af::tiny_plain<NumType,3>& v)
      {
        for(std::size_t i=0;i<9;) {
          for(std::size_t j=0;j<3;j++,i++) {
            this->elems[i] *= v[j];
          }
        }
        return *this;
      }

      //! Return a matrix with orthogonal base vectors.
      mat3 ortho() const;

      //! Decomposes the matrix into a rotation and scaling part.
      std::pair<mat3, vec3<NumType> >
      decompose() const;
  };

  // non-inline constructor
  template <typename NumType>
  mat3<NumType>::mat3(
    const NumType& angle, const af::tiny_plain<NumType,3>& axis)
  {
    NumType sqr_a = axis[0] * axis[0];
    NumType sqr_b = axis[1] * axis[1];
    NumType sqr_c = axis[2] * axis[2];
    NumType len2 = sqr_a + sqr_b + sqr_c;
    NumType k2 = std::cos(angle);
    NumType k1 = (1. - k2) / len2;
    NumType k3 = std::sin(angle) / std::sqrt(len2);
    NumType k1ab = k1 * axis[0] * axis[1];
    NumType k1ac = k1 * axis[0] * axis[2];
    NumType k1bc = k1 * axis[1] * axis[2];
    NumType k3a = k3 * axis[0];
    NumType k3b = k3 * axis[1];
    NumType k3c = k3 * axis[2];
    base_type result(
      k1*sqr_a+k2, k1ab-k3c, k1ac+k3b,
      k1ab+k3c, k1*sqr_b+k2, k1bc-k3a,
      k1ac-k3b, k1bc+k3a, k1*sqr_c+k2);
    std::copy(result.begin(), result.end(), this->begin());
  }

  // non-inline member function
  template <typename NumType>
  mat3<NumType>
  mat3<NumType>::ortho() const
  {
    vec3<NumType> x = get_column(0);
    vec3<NumType> y = get_column(1);
    vec3<NumType> z = get_column(2);
    NumType xl = x.length2();
    y = y - ((x * y) / xl) * x;
    z = z - ((x * z) / xl) * x;
    NumType yl = y.length2();
    z = z - ((y * z) / yl) * y;
    return mat3(x[0], y[0], z[0],
                x[1], y[1], z[1],
                x[2], y[2], z[2]);
  }

  // non-inline member function
  template <typename NumType>
  std::pair<mat3<NumType>, vec3<NumType> >
  mat3<NumType>::decompose() const
  {
    mat3 ro = ortho();
    vec3<NumType> x = ro.get_column(0);
    vec3<NumType> y = ro.get_column(1);
    vec3<NumType> z = ro.get_column(2);
    NumType xl = x.length();
    NumType yl = y.length();
    NumType zl = z.length();
    vec3<NumType> sc(xl, yl, zl);
    x /= xl;
    y /= yl;
    z /= zl;
    ro.set_column(0, x);
    ro.set_column(1, y);
    ro.set_column(2, z);
    if (ro.determinant() < NumType(0)) {
      ro.set_column(0, -x);
      sc[0] = -sc[0];
    }
    return std::make_pair(ro, sc);
  }

  //! Test equality.
  template <typename NumType>
  inline
  bool
  operator==(
    const mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType>
  inline
  bool
  operator==(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      if (lhs[i] != rhs   ) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of rhs == lhs.
  template <typename NumType>
  inline
  bool
  operator==(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      if (lhs    != rhs[i]) return false;
    }
    return true;
  }

  //! Test inequality.
  template <typename NumType>
  inline
  bool
  operator!=(
    const mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of rhs != lhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat3<NumType>
  operator+(
    const mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat3<NumType>
  operator+(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat3<NumType>
  operator+(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat3<NumType>
  operator-(
    const mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat3<NumType>
  operator-(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat3<NumType>
  operator-(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Matrix * matrix product.
  template <typename NumType>
  inline
  mat3<NumType>
  operator*(
    const mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    MatrixLite::multiply(lhs.begin(), rhs.begin(), 3, 3, 3, result.begin());
    return result;
  }

  //! Matrix * vector product.
  template <typename NumTypeMatrix,
            typename NumTypeVector>
  inline
  vec3<NumTypeVector>
  operator*(
    const mat3<NumTypeMatrix>& lhs,
    const af::tiny_plain<NumTypeVector,3>& rhs)
  {
    vec3<NumTypeVector> result;
    MatrixLite::multiply(lhs.begin(), rhs.begin(), 3, 3, 1, result.begin());
    return result;
  }

  //! Vector * matrix product.
  template <typename NumTypeVector,
            typename NumTypeMatrix>
  inline
  vec3<NumTypeVector>
  operator*(
    const af::tiny_plain<NumTypeVector,3>& lhs,
    const mat3<NumTypeMatrix>& rhs)
  {
    vec3<NumTypeVector> result;
    MatrixLite::multiply(lhs.begin(), rhs.begin(), 1, 3, 3, result.begin());
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  mat3<NumType>
  operator*(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  mat3<NumType>
  operator*(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  mat3<NumType>
  operator/(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  mat3<NumType>
  operator/(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  mat3<NumType>
  operator%(
    const mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  mat3<NumType>
  operator%(
    const      NumType & lhs,
    const mat3<NumType>& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise inplace addition.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator+=(
          mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise inplace addition.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator+=(
          mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace difference.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator-=(
          mat3<NumType>& lhs,
    const mat3<NumType>& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise inplace difference.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator-=(
          mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace multiplication.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator*=(
          mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace division.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator/=(
          mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace modulus operation.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator%=(
          mat3<NumType>& lhs,
    const      NumType & rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Element-wise unary minus.
  template <typename NumType>
  inline
  mat3<NumType>
  operator-(
    const mat3<NumType>& v)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  mat3<NumType>
  operator+(
    const mat3<NumType>& v)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = +v[i];
    }
    return result;
  }

} // namespace cctbx

#endif // CCTBX_MAT3_H
