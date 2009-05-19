#ifndef SCITBX_SYM_MAT2_H
#define SCITBX_SYM_MAT2_H

#include <scitbx/mat2.h>

namespace scitbx {

  static const std::size_t sym_mat2_storage_mapping[] = {0,2,
                                                         2,1};
  //! Symmetric 2x2 matrix class.
  template <typename NumType>
  class sym_mat2 : public af::tiny_plain<NumType, 3>
  {
    public:
      typedef typename af::tiny_plain<NumType, 3> base_type;

      //! Default constructor. Elements are not initialized.
      sym_mat2() {}
      //! Constructor.
      sym_mat2(NumType const& e00, NumType const& e11, NumType const& e01)
        : base_type(e00, e11, e01)
      {}
      //! Constructor.
      sym_mat2(base_type const& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      sym_mat2(const NumType* a)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = a[i];
      }
      //! Constructor for diagonal matrix.
      explicit
      sym_mat2(NumType const& diag)
        : base_type(diag,diag,0)
      {}
      //! Constructor for diagonal matrix.
      explicit
      sym_mat2(af::tiny_plain<NumType,2> const& diag)
        : base_type(diag[0],diag[1],0)
      {}
      //! Construction from full 2x2 matrix.
      /*! The off-diagonal elements of the new sym_mat2 are copied from
          the upper-right triangle of m.
          <p>
          An exception is thrown if the absolute value of the
          difference between any pair of off-diagonal elements is
          different from zero.
          <p>
          See also: mat2<>::is_symmetric()
       */
      explicit
      sym_mat2(mat2<NumType> const& m)
        : base_type(m[0],m[3],m[1])
      {
        SCITBX_ASSERT(m.is_symmetric());
      }
      //! Construction from full 2x2 matrix.
      /*! The off-diagonal elements of the new sym_mat2 are determined
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
          See also: mat2<>::is_symmetric()
       */
      explicit
      sym_mat2(mat2<NumType> const& m, NumType const& relative_tolerance)
        : base_type(m[0],m[3],
                    (m[1]+m[2])/2)
      {
        SCITBX_ASSERT(relative_tolerance < 0
          || m.is_symmetric(relative_tolerance));
      }

      //! Access elements with 2-dimensional indices.
      NumType const&
      operator()(std::size_t r, std::size_t c) const
      {
        return this->begin()[sym_mat2_storage_mapping[r * 2 + c]];
      }
      //! Access elements with 2-dimensional indices.
      NumType&
      operator()(std::size_t r, std::size_t c)
      {
        return this->begin()[sym_mat2_storage_mapping[r * 2 + c]];
      }

      //! Return diagonal elements.
      vec2<NumType>
      diagonal() const
      {
        return vec2<NumType>(this->begin());
      }

      //! Return the transposed matrix.
      /*! In analogy with mat2, return a new instance.
       */
      sym_mat2
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
        sym_mat2 const& m = *this;
        return   m(0,0) * m(1,1) - m(0,1) * m(1,0);
      }

      //! Return the transposed of the co-factor matrix.
      /*! The inverse matrix is obtained by dividing the result
          by the determinant().
       */
      sym_mat2
      co_factor_matrix_transposed() const
      {
        sym_mat2 const& m = *this;
        sym_mat2 result;
        result(0,0) =  m(1,1);
        result(1,1) =  m(0,0);
        result(0,1) = -m(0,1);
        return result;
      }

      //! Return the inverse matrix.
      /*! An exception is thrown if the matrix is not invertible,
          i.e. if the determinant() is zero.
       */
      sym_mat2
      inverse() const
      {
        NumType d = determinant();
        if (d == NumType(0)) throw error("Matrix is not invertible.");
        return co_factor_matrix_transposed() / d;
      }

      //! Tensor transform: c * (*this) * c.transpose()
      sym_mat2
      tensor_transform(mat2<NumType> const& c) const;

      //! Antisymmetric tensor transform: c * (*this) * c.transpose()
      /*! c is the antisymmetric matrix
            {{  0,  v0,  v1},
             {-v0,   0,  v2}
             {-v1, -v2,   0}}
       */
      sym_mat2
      antisymmetric_tensor_transform(
        NumType const& v0) const;

      //! Antisymmetric tensor transform: c * (*this) * c.transpose()
      /*! c is the antisymmetric matrix
            {{  0,    v[0],  v[1]},
             {-v[0],   0,    v[2]}
             {-v[1], -v[2],   0  }}
       */
      sym_mat2
      antisymmetric_tensor_transform(vec2<NumType> const& v) const
      {
        return antisymmetric_tensor_transform(v[0], v[1], v[2]);
      }

      //! Tensor transform: c.transpose() * (*this) * c
      sym_mat2
      tensor_transpose_transform(mat2<NumType> const& c) const;

      //! Sum of 4 element-wise products.
      inline
      NumType
      dot(sym_mat2 const& other) const
      {
        sym_mat2 const& m = *this;
        return m[0] * other[0]
             + m[1] * other[1]
             + 2 * (  m[2] * other[2] );
      }
  };

  // Constructor for mat2.
  template <typename NumType>
  inline
  mat2<NumType>::mat2(sym_mat2<NumType> const& m)
  : base_type(m[0], m[2],
              m[2], m[1])
  {}

  template <typename NumType>
  inline
  sym_mat2<NumType>
  mat2<NumType>::self_times_self_transpose() const
  {
    mat2<NumType> const& m = *this;
    return sym_mat2<NumType>(
      m[0]*m[0]+m[1]*m[1],
      m[2]*m[2]+m[3]*m[3],
      m[0]*m[2]+m[1]*m[3]);
  }

  template <typename NumType>
  inline
  sym_mat2<NumType>
  mat2<NumType>::self_transpose_times_self() const
  {
    mat2<NumType> const& m = *this;
    return sym_mat2<NumType>(
      m[0]*m[0]+m[2]*m[2],
      m[1]*m[1]+m[3]*m[3],
      m[0]*m[1]+m[2]*m[3]);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat2<NumType>
  sym_mat2<NumType>
  ::antisymmetric_tensor_transform(
    NumType const& v0) const
  {
    sym_mat2<NumType> const& t = *this;
    NumType v00 = v0 * v0;
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat2<NumType>(
       t[1]*v00,
       t[0]*v00,
      -t[2]*v00);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat2<NumType>
  sym_mat2<NumType>
  ::tensor_transform(mat2<NumType> const& c) const
  {
    mat2<NumType> ct = c * (*this);
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat2<NumType>(
      ct[0]*c[0]+ct[1]*c[1],
      ct[2]*c[2]+ct[3]*c[3],
      ct[0]*c[2]+ct[1]*c[3]);
  }

  // non-inline member function
  template <typename NumType>
  sym_mat2<NumType>
  sym_mat2<NumType>
  ::tensor_transpose_transform(mat2<NumType> const& c) const
  {
    sym_mat2<NumType> const& t = *this;
    mat2<NumType> ctt( // c.transpose() * (*this)
      c[0]*t[0]+c[2]*t[2],
      c[2]*t[1]+c[0]*t[2],
      c[1]*t[0]+c[3]*t[2],
      c[3]*t[1]+c[1]*t[2]);
    // The result is guaranteed to be a symmetric matrix.
    return sym_mat2<NumType>(
      ctt[0]*c[0]+ctt[1]*c[2],
      ctt[2]*c[1]+ctt[3]*c[3],
      ctt[0]*c[1]+ctt[1]*c[3]);
  }

  //! Test equality.
  template <typename NumType>
  inline
  bool
  operator==(
    sym_mat2<NumType> const& lhs,
    sym_mat2<NumType> const& rhs)
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
    sym_mat2<NumType> const& lhs,
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
    sym_mat2<NumType> const& rhs)
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
    sym_mat2<NumType> const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    sym_mat2<NumType> const& lhs,
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
    sym_mat2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator+(
    sym_mat2<NumType> const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator+(
    sym_mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator+(
    NumType const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator-(
    sym_mat2<NumType> const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator-(
    sym_mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator-(
    NumType const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Symmetric * symmetric matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat2<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    sym_mat2<NumTypeLhs> const& lhs,
    sym_mat2<NumTypeRhs> const& rhs)
  {
    typedef
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic
          result_element_type;
    result_element_type lhs_2__rhs_2_ = lhs[2]*rhs[2];
    return mat2<result_element_type>(
      lhs[0]*rhs[0]+lhs_2__rhs_2_,
      lhs[2]*rhs[1]+lhs[0]*rhs[2],
      lhs[2]*rhs[0]+lhs[1]*rhs[2],
      lhs[1]*rhs[1]+lhs_2__rhs_2_);
  }

  //! Square * symmetric matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat2<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    mat2<NumTypeLhs> const& lhs,
    sym_mat2<NumTypeRhs> const& rhs)
  {
    return mat2<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[2],
          lhs[1]*rhs[1]+lhs[0]*rhs[2],
          lhs[2]*rhs[0]+lhs[3]*rhs[2],
          lhs[2]*rhs[2]+lhs[3]*rhs[1]);
  }

  //! Symmetric * square matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat2<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    sym_mat2<NumTypeLhs> const& lhs,
    mat2<NumTypeRhs> const& rhs)
  {
    return mat2<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          rhs[0]*lhs[0]+rhs[2]*lhs[2],
          rhs[1]*lhs[0]+rhs[3]*lhs[2],
          rhs[2]*lhs[1]+rhs[0]*lhs[2],
          rhs[1]*lhs[2]+rhs[3]*lhs[1]);
  }

  //! Matrix * vector product.
  template <typename NumTypeMatrix,
            typename NumTypeVector>
  inline
  vec2<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    sym_mat2<NumTypeMatrix> const& lhs,
    af::tiny_plain<NumTypeVector,2> const& rhs)
  {
    return vec2<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[2]*rhs[1],
          lhs[2]*rhs[0]+lhs[1]*rhs[1]);
  }

  //! Vector * matrix product.
  template <typename NumTypeVector,
            typename NumTypeMatrix>
  inline
  vec2<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    af::tiny_plain<NumTypeVector,2> const& lhs,
    sym_mat2<NumTypeMatrix> const& rhs)
  {
    return vec2<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[2],
          lhs[1]*rhs[1]+lhs[0]*rhs[2]);
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator*(
    sym_mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator*(
    NumType const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator/(
    sym_mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator/(
    NumType const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator%(
    sym_mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator%(
    NumType const& lhs,
    sym_mat2<NumType> const& rhs)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  sym_mat2<NumType>&
  operator+=(
          sym_mat2<NumType>& lhs,
    sym_mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  sym_mat2<NumType>&
  operator+=(
          sym_mat2<NumType>& lhs,
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
  sym_mat2<NumType>&
  operator-=(
          sym_mat2<NumType>& lhs,
    sym_mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<3;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  sym_mat2<NumType>&
  operator-=(
          sym_mat2<NumType>& lhs,
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
  sym_mat2<NumType>&
  operator*=(
          sym_mat2<NumType>& lhs,
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
  sym_mat2<NumType>&
  operator/=(
          sym_mat2<NumType>& lhs,
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
  sym_mat2<NumType>&
  operator%=(
          sym_mat2<NumType>& lhs,
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
  sym_mat2<NumType>
  operator-(
    sym_mat2<NumType> const& v)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  sym_mat2<NumType>
  operator+(
    sym_mat2<NumType> const& v)
  {
    sym_mat2<NumType> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = +v[i];
    }
    return result;
  }

} // namespace scitbx

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::sym_mat2<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_SYM_MAT2_H
