/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy of cctbx/sym_mat3.h (Ralf W. Grosse-Kunstleve)
     2002 May: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_SYM_MAT3_H
#define SCITBX_SYM_MAT3_H

#include <scitbx/mat3.h>

namespace scitbx {

  static const std::size_t sym_mat3_storage_mapping[] = {0,3,4,
                                                         3,1,5,
                                                         4,5,2};
  //! Symmetric 3x3 matrix class.
  template <typename NumType>
  class sym_mat3 : public af::tiny_plain<NumType, 6>
  {
    public:
      typedef typename af::tiny_plain<NumType, 6> base_type;

      //! Default constructor. Elements are not initialized.
      sym_mat3() {}
      //! Constructor.
      sym_mat3(NumType const& e00, NumType const& e11, NumType const& e22,
               NumType const& e01, NumType const& e02, NumType const& e12)
        : base_type(e00, e11, e22, e01, e02, e12)
      {}
      //! Constructor.
      sym_mat3(base_type const& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      sym_mat3(const NumType* a)
      {
        for(std::size_t i=0;i<6;i++) this->elems[i] = a[i];
      }
      //! Constructor for diagonal matrix.
      explicit
      sym_mat3(NumType const& diag)
        : base_type(diag,diag,diag,0,0,0)
      {}
      //! Constructor for diagonal matrix.
      explicit
      sym_mat3(af::tiny_plain<NumType,3> const& diag)
        : base_type(diag[0],diag[1],diag[2],0,0,0)
      {}
      //! Construction from full 3x3 matrix.
      /*! The off-diagonal elements of the new sym_mat3 are copied from
          the upper-right triangle of m.
          <p>
          An exception is thrown if the absolute value of the
          difference between any pair of off-diagonal elements is
          different from zero.
          <p>
          See also: mat3<>::is_symmetric()
       */
      explicit
      sym_mat3(mat3<NumType> const& m)
        : base_type(m[0],m[4],m[8],m[1],m[2],m[5])
      {
        SCITBX_ASSERT(m.is_symmetric());
      }
      //! Construction from full 3x3 matrix.
      /*! The off-diagonal elements of the new sym_mat3 are determined
          as the averages of the corresponding off-diagonal elements
          of the input matrix m.
          <p>
          If relative_tolerance is greater than or equal to zero, it is
          used to check the input matrix m. An exception is thrown if
          the absolute value of the difference between any pair of
          off-diagonal elements is larger than
          max_abs*relative_tolerance, where max_abs is the maximum of
          the absolute values of the elements of m.
          <p>
          See also: mat3<>::is_symmetric()
       */
      explicit
      sym_mat3(mat3<NumType> const& m, NumType const& relative_tolerance)
        : base_type(m[0],m[4],m[8],
                    (m[1]+m[3])/2,
                    (m[2]+m[6])/2,
                    (m[5]+m[7])/2)
      {
        SCITBX_ASSERT(relative_tolerance < 0
          || m.is_symmetric(relative_tolerance));
      }

      //! Access elements with 2-dimensional indices.
      NumType const&
      operator()(std::size_t r, std::size_t c) const
      {
        return this->begin()[sym_mat3_storage_mapping[r * 3 + c]];
      }
      //! Access elements with 2-dimensional indices.
      NumType&
      operator()(std::size_t r, std::size_t c)
      {
        return this->begin()[sym_mat3_storage_mapping[r * 3 + c]];
      }

      //! Return diagonal elements.
      vec3<NumType>
      diagonal() const
      {
        return vec3<NumType>(this->begin());
      }

      //! Return the transposed matrix.
      /*! In analogy with mat3, return a new instance.
       */
      sym_mat3
      transpose()
      {
        return *this;
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
        sym_mat3 const& m = *this;
        return   m(0,0) * (m(1,1) * m(2,2) - m(1,2) * m(2,1))
               - m(0,1) * (m(1,0) * m(2,2) - m(1,2) * m(2,0))
               + m(0,2) * (m(1,0) * m(2,1) - m(1,1) * m(2,0));
      }

      //! Return the transposed of the co-factor matrix.
      /*! The inverse matrix is obtained by dividing the result
          by the determinant().
       */
      sym_mat3
      co_factor_matrix_transposed() const
      {
        sym_mat3 const& m = *this;
        sym_mat3 result;
        result(0,0) =  m(1,1) * m(2,2) - m(1,2) * m(2,1);
        result(1,1) =  m(0,0) * m(2,2) - m(0,2) * m(2,0);
        result(2,2) =  m(0,0) * m(1,1) - m(0,1) * m(1,0);
        result(0,1) = -m(0,1) * m(2,2) + m(0,2) * m(2,1);
        result(0,2) =  m(0,1) * m(1,2) - m(0,2) * m(1,1);
        result(1,2) = -m(0,0) * m(1,2) + m(0,2) * m(1,0);
        return result;
      }

      //! Return the inverse matrix.
      /*! An exception is thrown if the matrix is not invertible,
          i.e. if the determinant() is zero.
       */
      sym_mat3
      inverse() const
      {
        NumType d = determinant();
        if (d == NumType(0)) throw error("Matrix is not invertible.");
        return co_factor_matrix_transposed() / d;
      }
  };

  // Constructor for mat3.
  template <typename NumType>
  inline
  mat3<NumType>::mat3(sym_mat3<NumType> const& m)
  : base_type(m[0], m[3], m[4],
              m[3], m[1], m[5],
              m[4], m[5], m[2])
  {}

  //! Test equality.
  template <typename NumType>
  inline
  bool
  operator==(
    sym_mat3<NumType> const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType>
  inline
  bool
  operator==(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
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
    sym_mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      if (lhs    != rhs[i]) return false;
    }
    return true;
  }

  //! Test inequality.
  template <typename NumType>
  inline
  bool
  operator!=(
    sym_mat3<NumType> const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    sym_mat3<NumType> const& lhs,
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
    sym_mat3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator+(
    sym_mat3<NumType> const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator+(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator+(
    NumType const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator-(
    sym_mat3<NumType> const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator-(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator-(
    NumType const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Matrix * matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    sym_mat3<NumTypeLhs> const& lhs,
    sym_mat3<NumTypeRhs> const& rhs)
  {
    mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic> result;
    result.fill(0);
    for(std::size_t i=0;i<3;i++) {
    for(std::size_t j=0;j<3;j++) {
    for(std::size_t k=0;k<3;k++) {
      result(i,j) += lhs(i,k) * rhs(k,j);
    }}}
    return result;
  }

  //! Matrix * vector product.
  template <typename NumTypeMatrix,
            typename NumTypeVector>
  inline
  vec3<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    sym_mat3<NumTypeMatrix> const& lhs,
    af::tiny_plain<NumTypeVector,3> const& rhs)
  {
    vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic> result(0,0,0);
    for(std::size_t r=0;r<3;r++) {
    for(std::size_t c=0;c<3;c++) {
      result[r] += lhs(r,c) * rhs[c];
    }}
    return result;
  }

  //! Vector * matrix product.
  template <typename NumTypeVector,
            typename NumTypeMatrix>
  inline
  vec3<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    af::tiny_plain<NumTypeVector,3> const& lhs,
    sym_mat3<NumTypeMatrix> const& rhs)
  {
    vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic> result(0,0,0);
    for(std::size_t c=0;c<3;c++) {
    for(std::size_t r=0;r<3;r++) {
      result[c] += lhs[r] * rhs(r,c);
    }}
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator*(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator*(
    NumType const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator/(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator/(
    NumType const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator%(
    sym_mat3<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator%(
    NumType const& lhs,
    sym_mat3<NumType> const& rhs)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise inplace addition.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator+=(
          sym_mat3<NumType>& lhs,
    sym_mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise inplace addition.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator+=(
          sym_mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace difference.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator-=(
          sym_mat3<NumType>& lhs,
    sym_mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise inplace difference.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator-=(
          sym_mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace multiplication.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator*=(
          sym_mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace division.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator/=(
          sym_mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise inplace modulus operation.
  template <typename NumType>
  inline
  sym_mat3<NumType>&
  operator%=(
          sym_mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<6;i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Element-wise unary minus.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator-(
    sym_mat3<NumType> const& v)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  sym_mat3<NumType>
  operator+(
    sym_mat3<NumType> const& v)
  {
    sym_mat3<NumType> result;
    for(std::size_t i=0;i<6;i++) {
      result[i] = +v[i];
    }
    return result;
  }

} // namespace scitbx

#endif // SCITBX_SYM_MAT3_H
