#ifndef SCITBX_MAT2_H
#define SCITBX_MAT2_H

#include <utility>
#include <scitbx/error.h>
#include <scitbx/vec2.h>
#include <scitbx/array_family/tiny_reductions.h>

namespace scitbx {

  // forward declaration
  template <typename NumType>
  class sym_mat2;

  //! Matrix class (2x2).
  /*! This class represents a 2x2 matrix that can be used to store
      linear transformations.

      A 2D adaptation of the scitbx::mat3 class
   */
  template <typename NumType>
  class mat2 : public af::tiny_plain<NumType, 4>
  {
    public:
      typedef typename af::tiny_plain<NumType, 4> base_type;

      //! Default constructor. Elements are not initialized.
      mat2() {}
      //! Constructor.
      mat2(NumType const& e00, NumType const& e01,
           NumType const& e10, NumType const& e11)
        : base_type(e00, e01, e10, e11)
      {}
      //! Constructor.
      mat2(base_type const& a)
        : base_type(a)
      {}
      //! Constructor.
      explicit
      mat2(const NumType* a)
      {
        for(std::size_t i=0;i<4;i++) this->elems[i] = a[i];
      }
      //! Constructor for diagonal matrix.
      explicit
      mat2(NumType const& diag)
        : base_type(diag,0,0,diag)
      {}
      //! Constructor for diagonal matrix.
      explicit
      mat2(NumType const& diag0, NumType const& diag1)
        : base_type(diag0,0,0,diag1)
      {}
      //! Constructor for diagonal matrix.
      explicit
      mat2(af::tiny_plain<NumType,2> const& diag)
        : base_type(diag[0],0,0,diag[1])
      {}
      //! Construction from symmetric matrix.
      explicit
      inline
      mat2(sym_mat2<NumType> const& m);

      //! Access elements with 2-dimensional indices.
      NumType const&
      operator()(std::size_t r, std::size_t c) const
      {
        return this->elems[r * 2 + c];
      }
      //! Access elements with 2-dimensional indices.
      NumType&
      operator()(std::size_t r, std::size_t c)
      {
        return this->elems[r * 2 + c];
      }

      //! Return a row.
      vec2<NumType>
      get_row(std::size_t i) const
      {
        return vec2<NumType>(&this->elems[i * 2]);
      }

      //! Set a row.
      void
      set_row(std::size_t i, af::tiny_plain<NumType,2> const& v)
      {
        std::copy(v.begin(), v.end(), &this->elems[i * 2]);
      }

      //! Swap two rows in place.
      void
      swap_rows(std::size_t i1, std::size_t i2)
      {
        std::swap_ranges(&(*this)(i1,0), &(*this)(i1+1,0), &(*this)(i2,0));
      }

      //! Return a column.
      vec2<NumType>
      get_column(std::size_t i) const
      {
        vec2<NumType> result;
        for(std::size_t j=0;j<2;j++) result[j] = this->elems[j * 2 + i];
        return result;
      }

      //! Set a column.
      void
      set_column(std::size_t i, af::tiny_plain<NumType,2> const& v)
      {
        for(std::size_t j=0;j<2;j++) this->elems[j * 2 + i] = v[j];
      }

      //! Swap two columns in place.
      void
      swap_columns(std::size_t i1, std::size_t i2)
      {
        for(std::size_t i=0;i<4;i+=2) {
          std::swap(this->elems[i + i1], this->elems[i + i2]);
        }
      }

      //! Return diagonal elements.
      vec2<NumType>
      diagonal() const
      {
        mat2 const& m = *this;
        return vec2<NumType>(m[0], m[3]);
      }

      //! Return the transposed matrix.
      mat2
      transpose() const
      {
        mat2 const& m = *this;
        return mat2(m[0], m[2],
                    m[1], m[3]);
      }

      //! Return trace (sum of diagonal elements).
      NumType
      trace() const
      {
        mat2 const& m = *this;
        return m[0] + m[3];
      }

      //! Return determinant.
      NumType
      determinant() const
      {
        mat2 const& m = *this;
        return   m[0] * m[3] - m[1] * m[2];
      }

      //! Maximum of the absolute values of the elements of this matrix.
      NumType
      max_abs() const
      {
        return af::max_absolute(this->const_ref());
      }

      //! Test for symmetric matrix.
      /*! Returns false iff the absolute value of the difference between
          any pair of off-diagonal elements is different from zero.
       */
      bool
      is_symmetric() const
      {
        mat2 const& m = *this;
        return    m[1] == m[2];
      }

      //! Test for symmetric matrix.
      /*! Returns false iff the absolute value of the difference between
          any pair of off-diagonal elements is larger than
          max_abs()*relative_tolerance.
       */
      bool
      is_symmetric(NumType const& relative_tolerance) const
      {
        mat2 const& m = *this;
        NumType tolerance = max_abs() * relative_tolerance;
        return    fn::approx_equal(m[1], m[2], tolerance);
      }

      bool
      is_diagonal() const
      {
        mat2 const& m = *this;
        return (              m[1]==0
                && m[2]==0);
      }

      //! Return the transposed of the co-factor matrix.
      /*! The inverse matrix is obtained by dividing the result
          by the determinant().
       */
      mat2
      co_factor_matrix_transposed() const
      {
        mat2 const& m = *this;
        return mat2(m[3],-m[1],-m[2],m[0]);
      }

      //! Return the inverse matrix.
      /*! An exception is thrown if the matrix is not invertible,
          i.e. if the determinant() is zero.
       */
      mat2
      inverse() const
      {
        NumType d = determinant();
        if (d == NumType(0)) throw error("Matrix is not invertible.");
        return co_factor_matrix_transposed() / d;
      }

      //! Returns the inverse matrix, after minimizing the error numerically.
      /*! Here's the theory:
          M*M^-1 = I-E, where E is the error
          M*M^-1*(I+E) = (I-E)*(I+E)
          M*(M^-1*(I+E)) = I^2-E^2
          M*(M^-1*(I+E)) = I-E^2
                      let M^-1*(I+E)  = M1
                      let E^2         = E2
          M*M1*(I+E2) = (I-E2)*(I+E2)
          M*M2 = I-E4
          M*Mi = I-E2^i
          Supposedly this will drive the error pretty low after
          only a few repetitions. The error rate should be ~E^(2^iterations),
          which I think is pretty good. This assumes that E is "<< 1",
          whatever that means. Attributed to Judah I. Rosenblatt.

          2*I - (I-E) ==> 2*I - I + E = I + E
       */
      mat2
      error_minimizing_inverse ( std::size_t iterations ) const
      {
        mat2 inverse = this->inverse();
        if ( 0 == iterations )
                      return inverse;
              mat2 two_diagonal(2);
        for ( std::size_t i=0; i<iterations; ++i )
                inverse = inverse * (two_diagonal - this*inverse);
        return inverse;
      }

      //! Scale matrix in place.
      /*! Each row of this is multiplied element-wise with v.
       */
      mat2&
      scale(af::tiny_plain<NumType,2> const& v)
      {
        for(std::size_t i=0;i<4;) {
          for(std::size_t j=0;j<2;j++,i++) {
            this->elems[i] *= v[j];
          }
        }
        return *this;
      }

      //! Return a matrix with orthogonal base vectors.
      mat2 ortho() const;

      //! Decomposes the matrix into a rotation and scaling part.
      std::pair<mat2, vec2<NumType> >
      decompose() const;

      //! (*this) * this->transpose().
      inline
      sym_mat2<NumType>
      self_times_self_transpose() const;

      //! this->transpose() * (*this).
      inline
      sym_mat2<NumType>
      self_transpose_times_self() const;

      //! Sum of element-wise products.
      inline
      NumType
      dot(mat2 const& other) const
      {
        mat2 const& m = *this;
        return m[0] * other[0]
             + m[1] * other[1]
             + m[2] * other[2]
             + m[3] * other[3];
      }
  };

  // non-inline member function
  template <typename NumType>
  mat2<NumType>
  mat2<NumType>::ortho() const
  {
    vec2<NumType> x = get_column(0);
    vec2<NumType> y = get_column(1);
    NumType xl = x.length_sq();
    y = y - ((x * y) / xl) * x;
    return mat2(x[0], y[0],
                x[1], y[1]);
  }

  // non-inline member function
  template <typename NumType>
  std::pair<mat2<NumType>, vec2<NumType> >
  mat2<NumType>::decompose() const
  {
    mat2 ro = ortho();
    vec2<NumType> x = ro.get_column(0);
    vec2<NumType> y = ro.get_column(1);
    NumType xl = x.length();
    NumType yl = y.length();
    vec2<NumType> sc(xl, yl);
    x /= xl;
    y /= yl;
    ro.set_column(0, x);
    ro.set_column(1, y);
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
    mat2<NumType> const& lhs,
    mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType>
  inline
  bool
  operator==(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
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
    mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      if (lhs    != rhs[i]) return false;
    }
    return true;
  }

  //! Test inequality.
  template <typename NumType>
  inline
  bool
  operator!=(
    mat2<NumType> const& lhs,
    mat2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    mat2<NumType> const& lhs,
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
    mat2<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat2<NumType>
  operator+(
    mat2<NumType> const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] + rhs[i];
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat2<NumType>
  operator+(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] + rhs   ;
    }
    return result;
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat2<NumType>
  operator+(
    NumType const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs    + rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat2<NumType>
  operator-(
    mat2<NumType> const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] - rhs[i];
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat2<NumType>
  operator-(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] - rhs   ;
    }
    return result;
  }

  //! Element-wise difference.
  template <typename NumType>
  inline
  mat2<NumType>
  operator-(
    NumType const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs    - rhs[i];
    }
    return result;
  }

  //! Matrix * matrix product.
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat2<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  operator*(
    mat2<NumTypeLhs> const& lhs,
    mat2<NumTypeRhs> const& rhs)
  {
    return mat2<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[2],
          lhs[0]*rhs[1]+lhs[1]*rhs[3],
          lhs[2]*rhs[0]+lhs[3]*rhs[2],
          lhs[2]*rhs[1]+lhs[3]*rhs[3]);
  }

  //! Matrix * vector product.
  template <typename NumTypeMatrix,
            typename NumTypeVector>
  inline
  vec2<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    mat2<NumTypeMatrix> const& lhs,
    af::tiny_plain<NumTypeVector,2> const& rhs)
  {
    return vec2<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[1],
          lhs[2]*rhs[0]+lhs[3]*rhs[1]);
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
    mat2<NumTypeMatrix> const& rhs)
  {
    return vec2<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[2],
          lhs[0]*rhs[1]+lhs[1]*rhs[3]);
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  mat2<NumType>
  operator*(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] * rhs   ;
    }
    return result;
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  mat2<NumType>
  operator*(
    NumType const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs    * rhs[i];
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  mat2<NumType>
  operator/(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] / rhs   ;
    }
    return result;
  }

  //! Element-wise division.
  template <typename NumType>
  inline
  mat2<NumType>
  operator/(
    NumType const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs    / rhs[i];
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  mat2<NumType>
  operator%(
    mat2<NumType> const& lhs,
    NumType const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs[i] % rhs   ;
    }
    return result;
  }

  //! Element-wise modulus operation.
  template <typename NumType>
  inline
  mat2<NumType>
  operator%(
    NumType const& lhs,
    mat2<NumType> const& rhs)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator+=(
    mat2<NumType>& lhs,
    mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator+=(
    mat2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator-=(
    mat2<NumType>& lhs,
    mat2<NumType> const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator-=(
    mat2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place multiplication.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator*=(
    mat2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place division.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator/=(
    mat2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place modulus operation.
  template <typename NumType>
  inline
  mat2<NumType>&
  operator%=(
    mat2<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<4;i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Element-wise unary minus.
  template <typename NumType>
  inline
  mat2<NumType>
  operator-(
    mat2<NumType> const& v)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = -v[i];
    }
    return result;
  }

  //! Element-wise unary plus.
  template <typename NumType>
  inline
  mat2<NumType>
  operator+(
    mat2<NumType> const& v)
  {
    mat2<NumType> result;
    for(std::size_t i=0;i<4;i++) {
      result[i] = +v[i];
    }
    return result;
  }

} // namespace scitbx

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::mat2<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_MAT2_H
