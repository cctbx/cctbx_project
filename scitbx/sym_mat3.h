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

      //! Tensor transform: c * (*this) * c.transpose()
      sym_mat3
      tensor_transform(mat3<NumType> const& c) const;

      //! Antisymmetric tensor transform: c * (*this) * c.transpose()
      /*! c is the antisymmetric matrix
            {{  0,  v0,  v1},
             {-v0,   0,  v2}
             {-v1, -v2,   0}}
       */
      sym_mat3
      antisymmetric_tensor_transform(
        NumType const& v0,
        NumType const& v1,
        NumType const& v2) const;

      //! Antisymmetric tensor transform: c * (*this) * c.transpose()
      /*! c is the antisymmetric matrix
            {{  0,    v[0],  v[1]},
             {-v[0],   0,    v[2]}
             {-v[1], -v[2],   0  }}
       */
      sym_mat3
      antisymmetric_tensor_transform(vec3<NumType> const& v) const
      {
        return antisymmetric_tensor_transform(v[0], v[1], v[2]);
      }

      //! Tensor transform: c.transpose() * (*this) * c
      sym_mat3
      tensor_transpose_transform(mat3<NumType> const& c) const;

      //! Sum of 9 element-wise products.
      inline
      NumType
      dot(sym_mat3 const& other) const
      {
        sym_mat3 const& m = *this;
        return m[0] * other[0]
             + m[1] * other[1]
             + m[2] * other[2]
             + 2 * (  m[3] * other[3]
                    + m[4] * other[4]
                    + m[5] * other[5]);
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

  template <typename NumType>
  inline
  sym_mat3<NumType>
  mat3<NumType>::self_times_self_transpose() const
  {
    mat3<NumType> const& m = *this;
    return sym_mat3<NumType>(
      m[0]*m[0]+m[1]*m[1]+m[2]*m[2],
      m[3]*m[3]+m[4]*m[4]+m[5]*m[5],
      m[6]*m[6]+m[7]*m[7]+m[8]*m[8],
      m[0]*m[3]+m[1]*m[4]+m[2]*m[5],
      m[0]*m[6]+m[1]*m[7]+m[2]*m[8],
      m[3]*m[6]+m[4]*m[7]+m[5]*m[8]);
  }

  template <typename NumType>
  inline
  sym_mat3<NumType>
  mat3<NumType>::self_transpose_times_self() const
  {
    mat3<NumType> const& m = *this;
    return sym_mat3<NumType>(
      m[0]*m[0]+m[3]*m[3]+m[6]*m[6],
      m[1]*m[1]+m[4]*m[4]+m[7]*m[7],
      m[2]*m[2]+m[5]*m[5]+m[8]*m[8],
      m[0]*m[1]+m[3]*m[4]+m[6]*m[7],
      m[0]*m[2]+m[3]*m[5]+m[6]*m[8],
      m[1]*m[2]+m[4]*m[5]+m[7]*m[8]);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat3<NumType>
  sym_mat3<NumType>
  ::antisymmetric_tensor_transform(
    NumType const& v0,
    NumType const& v1,
    NumType const& v2) const
  {
    sym_mat3<NumType> const& t = *this;
    NumType v00 = v0 * v0;
    NumType v11 = v1 * v1;
    NumType v22 = v2 * v2;
    NumType v01 = v0 * v1;
    NumType v02 = v0 * v2;
    NumType v12 = v1 * v2;
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat3<NumType>(
       t[2]*v11 + 2*t[5]*v01 + t[1]*v00,
       t[2]*v22 - 2*t[4]*v02 + t[0]*v00,
       t[1]*v22 + 2*t[3]*v12 + t[0]*v11,
       t[2]*v12 + t[5]*v02 - t[4]*v01 - t[3]*v00,
      -t[5]*v12 - t[1]*v02 - t[4]*v11 - t[3]*v01,
      -t[5]*v22 - t[4]*v12 + t[3]*v02 + t[0]*v01);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat3<NumType>
  sym_mat3<NumType>
  ::tensor_transform(mat3<NumType> const& c) const
  {
    mat3<NumType> ct = c * (*this);
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat3<NumType>(
      ct[0]*c[0]+ct[1]*c[1]+ct[2]*c[2],
      ct[3]*c[3]+ct[4]*c[4]+ct[5]*c[5],
      ct[6]*c[6]+ct[7]*c[7]+ct[8]*c[8],
      ct[0]*c[3]+ct[1]*c[4]+ct[2]*c[5],
      ct[0]*c[6]+ct[1]*c[7]+ct[2]*c[8],
      ct[3]*c[6]+ct[4]*c[7]+ct[5]*c[8]);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat3<NumType>
  sym_mat3<NumType>
  ::tensor_transpose_transform(mat3<NumType> const& c) const
  {
    sym_mat3<NumType> const& t = *this;
    mat3<NumType> ctt( // c.transpose() * (*this)
      c[0]*t[0]+c[3]*t[3]+c[6]*t[4],
      c[3]*t[1]+c[0]*t[3]+c[6]*t[5],
      c[6]*t[2]+c[0]*t[4]+c[3]*t[5],
      c[1]*t[0]+c[4]*t[3]+c[7]*t[4],
      c[4]*t[1]+c[1]*t[3]+c[7]*t[5],
      c[7]*t[2]+c[1]*t[4]+c[4]*t[5],
      c[2]*t[0]+c[5]*t[3]+c[8]*t[4],
      c[5]*t[1]+c[2]*t[3]+c[8]*t[5],
      c[8]*t[2]+c[2]*t[4]+c[5]*t[5]);
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat3<NumType>(
      ctt[0]*c[0]+ctt[1]*c[3]+ctt[2]*c[6],
      ctt[3]*c[1]+ctt[4]*c[4]+ctt[5]*c[7],
      ctt[6]*c[2]+ctt[7]*c[5]+ctt[8]*c[8],
      ctt[0]*c[1]+ctt[1]*c[4]+ctt[2]*c[7],
      ctt[0]*c[2]+ctt[1]*c[5]+ctt[2]*c[8],
      ctt[3]*c[2]+ctt[4]*c[5]+ctt[5]*c[8]);
  }

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

  //! Symmetric * symmetric matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    sym_mat3<NumTypeLhs> const& lhs,
    sym_mat3<NumTypeRhs> const& rhs)
  {
    typedef
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic
          result_element_type;
    result_element_type lhs_3__rhs_3_ = lhs[3]*rhs[3];
    result_element_type lhs_4__rhs_4_ = lhs[4]*rhs[4];
    result_element_type lhs_5__rhs_5_ = lhs[5]*rhs[5];
    return mat3<result_element_type>(
      lhs[0]*rhs[0]+lhs_3__rhs_3_+lhs_4__rhs_4_,
      lhs[3]*rhs[1]+lhs[0]*rhs[3]+lhs[4]*rhs[5],
      lhs[4]*rhs[2]+lhs[0]*rhs[4]+lhs[3]*rhs[5],
      lhs[3]*rhs[0]+lhs[1]*rhs[3]+lhs[5]*rhs[4],
      lhs[1]*rhs[1]+lhs_3__rhs_3_+lhs_5__rhs_5_,
      lhs[5]*rhs[2]+lhs[3]*rhs[4]+lhs[1]*rhs[5],
      lhs[4]*rhs[0]+lhs[5]*rhs[3]+lhs[2]*rhs[4],
      lhs[5]*rhs[1]+lhs[4]*rhs[3]+lhs[2]*rhs[5],
      lhs[2]*rhs[2]+lhs_4__rhs_4_+lhs_5__rhs_5_);
  }

  //! Square * symmetric matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    mat3<NumTypeLhs> const& lhs,
    sym_mat3<NumTypeRhs> const& rhs)
  {
    return mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[4],
          lhs[1]*rhs[1]+lhs[0]*rhs[3]+lhs[2]*rhs[5],
          lhs[2]*rhs[2]+lhs[0]*rhs[4]+lhs[1]*rhs[5],
          lhs[3]*rhs[0]+lhs[4]*rhs[3]+lhs[5]*rhs[4],
          lhs[4]*rhs[1]+lhs[3]*rhs[3]+lhs[5]*rhs[5],
          lhs[5]*rhs[2]+lhs[3]*rhs[4]+lhs[4]*rhs[5],
          lhs[6]*rhs[0]+lhs[7]*rhs[3]+lhs[8]*rhs[4],
          lhs[7]*rhs[1]+lhs[6]*rhs[3]+lhs[8]*rhs[5],
          lhs[8]*rhs[2]+lhs[6]*rhs[4]+lhs[7]*rhs[5]);
  }

  //! Symmetric * square matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    sym_mat3<NumTypeLhs> const& lhs,
    mat3<NumTypeRhs> const& rhs)
  {
    return mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          rhs[0]*lhs[0]+rhs[3]*lhs[3]+rhs[6]*lhs[4],
          rhs[1]*lhs[0]+rhs[4]*lhs[3]+rhs[7]*lhs[4],
          rhs[2]*lhs[0]+rhs[5]*lhs[3]+rhs[8]*lhs[4],
          rhs[3]*lhs[1]+rhs[0]*lhs[3]+rhs[6]*lhs[5],
          rhs[4]*lhs[1]+rhs[1]*lhs[3]+rhs[7]*lhs[5],
          rhs[5]*lhs[1]+rhs[2]*lhs[3]+rhs[8]*lhs[5],
          rhs[6]*lhs[2]+rhs[0]*lhs[4]+rhs[3]*lhs[5],
          rhs[7]*lhs[2]+rhs[1]*lhs[4]+rhs[4]*lhs[5],
          rhs[8]*lhs[2]+rhs[2]*lhs[4]+rhs[5]*lhs[5]);
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
    return vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[3]*rhs[1]+lhs[4]*rhs[2],
          lhs[3]*rhs[0]+lhs[1]*rhs[1]+lhs[5]*rhs[2],
          lhs[4]*rhs[0]+lhs[5]*rhs[1]+lhs[2]*rhs[2]);
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
    return vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[4],
          lhs[1]*rhs[1]+lhs[0]*rhs[3]+lhs[2]*rhs[5],
          lhs[2]*rhs[2]+lhs[0]*rhs[4]+lhs[1]*rhs[5]);
  }

  /// Linear form * (matrix as 6-vector)
  template <typename NumTypeVector,
            typename NumTypeMatrix>
  inline
  typename af::binary_operator_traits<NumTypeVector,
                                      NumTypeMatrix>::arithmetic
  operator*(af::tiny_plain<NumTypeVector, 6> const &form,
            sym_mat3<NumTypeMatrix> const &matrix_as_vec6)
  {
    typename af::binary_operator_traits<NumTypeVector,
                                        NumTypeMatrix>::arithmetic
    result = 0;
    for (int i=0; i<6; ++i) result += form[i]*matrix_as_vec6[i];
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

  //! Element-wise in-place addition.
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

  //! Element-wise in-place addition.
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

  //! Element-wise in-place difference.
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

  //! Element-wise in-place difference.
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

  //! Element-wise in-place multiplication.
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

  //! Element-wise in-place division.
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

  //! Element-wise in-place modulus operation.
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

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::sym_mat3<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_SYM_MAT3_H
