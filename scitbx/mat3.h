#ifndef SCITBX_MAT3_H
#define SCITBX_MAT3_H

#include <utility>
#include <scitbx/error.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_reductions.h>

namespace scitbx {

  // forward declaration
  template <typename NumType>
  class sym_mat3;

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
      mat3(NumType const& e00, NumType const& e01, NumType const& e02,
           NumType const& e10, NumType const& e11, NumType const& e12,
           NumType const& e20, NumType const& e21, NumType const& e22)
        : base_type(e00, e01, e02, e10, e11, e12, e20, e21, e22)
      {}
      //! Constructor.
      mat3(base_type const& a)
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
      mat3(NumType const& diag)
        : base_type(diag,0,0,0,diag,0,0,0,diag)
      {}
      //! Constructor for diagonal matrix.
      explicit
      mat3(NumType const& diag0, NumType const& diag1, NumType const& diag2)
        : base_type(diag0,0,0,0,diag1,0,0,0,diag2)
      {}
      //! Constructor for diagonal matrix.
      explicit
      mat3(af::tiny_plain<NumType,3> const& diag)
        : base_type(diag[0],0,0,0,diag[1],0,0,0,diag[2])
      {}
      //! Construction from symmetric matrix.
      explicit
      inline
      mat3(sym_mat3<NumType> const& m);

      //! Access elements with 2-dimensional indices.
      NumType const&
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
      set_row(std::size_t i, af::tiny_plain<NumType,3> const& v)
      {
        std::copy(v.begin(), v.end(), &this->elems[i * 3]);
      }

      //! Swap two rows in place.
      void
      swap_rows(std::size_t i1, std::size_t i2)
      {
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
      set_column(std::size_t i, af::tiny_plain<NumType,3> const& v)
      {
        for(std::size_t j=0;j<3;j++) this->elems[j * 3 + i] = v[j];
      }

      //! Swap two columns in place.
      void
      swap_columns(std::size_t i1, std::size_t i2)
      {
        for(std::size_t i=0;i<9;i+=3) {
          std::swap(this->elems[i + i1], this->elems[i + i2]);
        }
      }

      //! Return diagonal elements.
      vec3<NumType>
      diagonal() const
      {
        mat3 const& m = *this;
        return vec3<NumType>(m[0], m[4], m[8]);
      }

      //! Return the transposed matrix.
      mat3
      transpose() const
      {
        mat3 const& m = *this;
        return mat3(m[0], m[3], m[6],
                    m[1], m[4], m[7],
                    m[2], m[5], m[8]);
      }

      //! Return trace (sum of diagonal elements).
      NumType
      trace() const
      {
        mat3 const& m = *this;
        return m[0] + m[4] + m[8];
      }

      //! Return determinant.
      NumType
      determinant() const
      {
        mat3 const& m = *this;
        return   m[0] * (m[4] * m[8] - m[5] * m[7])
               - m[1] * (m[3] * m[8] - m[5] * m[6])
               + m[2] * (m[3] * m[7] - m[4] * m[6]);
      }

      //! Maximum of the absolute values of the elements of this matrix.
      NumType
      max_abs() const
      {
        return af::max_absolute(this->const_ref());
      }

      //! Test for symmetric matrix.
      /*! Returns false if the absolute value of the difference between
          any pair of off-diagonal elements is different from zero.
       */
      bool
      is_symmetric() const
      {
        mat3 const& m = *this;
        return    m[1] == m[3]
               && m[2] == m[6]
               && m[5] == m[7];
      }

      //! Test for symmetric matrix.
      /*! Returns false if the absolute value of the difference between
          any pair of off-diagonal elements is larger than
          max_abs()*relative_tolerance.
       */
      bool
      is_symmetric(NumType const& relative_tolerance) const
      {
        mat3 const& m = *this;
        NumType tolerance = max_abs() * relative_tolerance;
        return    fn::approx_equal(m[1], m[3], tolerance)
               && fn::approx_equal(m[2], m[6], tolerance)
               && fn::approx_equal(m[5], m[7], tolerance);
      }

      bool
      is_diagonal() const
      {
        mat3 const& m = *this;
        return (              m[1]==0 && m[2]==0
                && m[3]==0            && m[5]==0
                && m[6]==0 && m[7]==0);
      }

      //! Return the transposed of the co-factor matrix.
      /*! The inverse matrix is obtained by dividing the result
          by the determinant().
       */
      mat3
      co_factor_matrix_transposed() const
      {
        mat3 const& m = *this;
        return mat3(
           m[4] * m[8] - m[5] * m[7],
          -m[1] * m[8] + m[2] * m[7],
           m[1] * m[5] - m[2] * m[4],
          -m[3] * m[8] + m[5] * m[6],
           m[0] * m[8] - m[2] * m[6],
          -m[0] * m[5] + m[2] * m[3],
           m[3] * m[7] - m[4] * m[6],
          -m[0] * m[7] + m[1] * m[6],
           m[0] * m[4] - m[1] * m[3]);
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
      mat3
      error_minimizing_inverse ( std::size_t iterations ) const
      {
        mat3 inverse = this->inverse();
        if ( 0 == iterations )
                      return inverse;
              mat3 two_diagonal(2);
        for ( std::size_t i=0; i<iterations; ++i )
                inverse = inverse * (two_diagonal - this*inverse);
        return inverse;
      }

      //! Scale matrix in place.
      /*! Each row of this is multiplied element-wise with v.
       */
      mat3&
      scale(af::tiny_plain<NumType,3> const& v)
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

      //! (*this) * this->transpose().
      inline
      sym_mat3<NumType>
      self_times_self_transpose() const;

      //! this->transpose() * (*this).
      inline
      sym_mat3<NumType>
      self_transpose_times_self() const;

      //! Sum of element-wise products.
      inline
      NumType
      dot(mat3 const& other) const
      {
        mat3 const& m = *this;
        return m[0] * other[0]
             + m[1] * other[1]
             + m[2] * other[2]
             + m[3] * other[3]
             + m[4] * other[4]
             + m[5] * other[5]
             + m[6] * other[6]
             + m[7] * other[7]
             + m[8] * other[8];
      }

      //! Matrix associated with vector cross product.
      /*! a.cross(b) is equivalent to cross_product_matrix(a) * b.
          Useful for simplification of equations. Used frequently in
          robotics and classical mechanics literature.
       */
      static
      mat3
      cross_product_matrix(
        vec3<NumType> const& v)
      {
        return mat3(
              0, -v[2],  v[1],
           v[2],     0, -v[0],
          -v[1],  v[0],     0);
      }

    /// Matrix associated with symmetric rank-2 tensor transform
    /** The change of basis of such a tensor U reads
        \f[ U' = R U R^T \f]
        where \f$U, U'\f$ are represented as symmetric matrices.
        If, on the other hand, \f$U\f$ is represented as a vector
        \f[ u = (u_{11}, u_{22}, u_{33}, u_{12}, u_{13}, u_{23}) \f]
        and similarly \f$U'\f$, then the change of basis may be written
        as a linear transformation
        \f[ u' = P u \]
        where \f$P\f$ reads
        \f[ P = \begin{bmatrix}
         r_{11}^2 & r_{12}^2 & r_{13}^2 &
          2 r_{11} r_{12} & 2 r_{11} r_{13} & 2 r_{12} r_{13} \\
         r_{21}^2 & r_{22}^2 & r_{23}^2 &
          2 r_{21} r_{22} & 2 r_{21} r_{23} & 2 r_{22} r_{23} \\
         r_{31}^2 & r_{32}^2 & r_{33}^2 &
          2 r_{31} r_{32} & 2 r_{31} r_{33} & 2 r_{32} r_{33} \\
         r_{11} r_{21} & r_{12} r_{22} & r_{13} r_{23} &
          r_{12} r_{21} + r_{11} r_{22} &
          r_{13} r_{21} + r_{11} r_{23} &
          r_{13} r_{22} + r_{12} r_{23} \\
         r_{11} r_{31} & r_{12} r_{32} & r_{13} r_{33} &
          r_{12} r_{31} + r_{11} r_{32} &
          r_{13} r_{31} + r_{11} r_{33} &
          r_{13} r_{32} + r_{12} r_{33} \\
         r_{21} r_{31} & r_{22} r_{32} & r_{23} r_{33} &
          r_{22} r_{31} + r_{21} r_{32} &
          r_{23} r_{31} + r_{21} r_{33} &
          r_{23} r_{32} + r_{22} r_{33}
        \end{bmatrix}
        C.f. scitbx/matrix/tensor_transform_as_linear_map.nb
        for the Mathematica code used to obtain this result.

        The result is a matrix stored by row.
     */
    af::tiny<NumType, 6*6> tensor_transform_matrix() const {
      af::tiny<NumType, 6*6> result;
      NumType *p = result.begin();

      mat3 const &m = *this;
      NumType r11 = m[0], r12 = m[1], r13 = m[2],
              r21 = m[3], r22 = m[4], r23 = m[5],
              r31 = m[6], r32 = m[7], r33 = m[8];

      // 1st row
      *p++ = r11*r11;
      *p++ = r12*r12;
      *p++ = r13*r13;
      *p++ = 2*r11*r12;
      *p++ = 2*r11*r13;
      *p++ = 2*r12*r13;

      // 2nd row
      *p++ = r21*r21;
      *p++ = r22*r22;
      *p++ = r23*r23;
      *p++ = 2*r21*r22;
      *p++ = 2*r21*r23;
      *p++ = 2*r22*r23;

      // 3rd row
      *p++ = r31*r31;
      *p++ = r32*r32;
      *p++ = r33*r33;
      *p++ = 2*r31*r32;
      *p++ = 2*r31*r33;
      *p++ = 2*r32*r33;

      // 4th row
      *p++ = r11*r21;
      *p++ = r12*r22;
      *p++ = r13*r23;
      *p++ = r12*r21 + r11*r22;
      *p++ = r13*r21 + r11*r23;
      *p++ = r13*r22 + r12*r23;

      // 5th row
      *p++ = r11*r31;
      *p++ = r12*r32;
      *p++ = r13*r33;
      *p++ = r12*r31 + r11*r32;
      *p++ = r13*r31 + r11*r33;
      *p++ = r13*r32 + r12*r33;

      // 6th row
      *p++ = r21*r31;
      *p++ = r22*r32;
      *p++ = r23*r33;
      *p++ = r22*r31 + r21*r32;
      *p++ = r23*r31 + r21*r33;
      *p++ = r23*r32 + r22*r33;

      return result;
    }
  };

  // non-inline member function
  template <typename NumType>
  mat3<NumType>
  mat3<NumType>::ortho() const
  {
    vec3<NumType> x = get_column(0);
    vec3<NumType> y = get_column(1);
    vec3<NumType> z = get_column(2);
    NumType xl = x.length_sq();
    y = y - ((x * y) / xl) * x;
    z = z - ((x * z) / xl) * x;
    NumType yl = y.length_sq();
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
    mat3<NumType> const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    mat3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType>
  inline
  bool
  operator!=(
    mat3<NumType> const& lhs,
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
    mat3<NumType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise addition.
  template <typename NumType>
  inline
  mat3<NumType>
  operator+(
    mat3<NumType> const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
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
    mat3<NumTypeLhs> const& lhs,
    mat3<NumTypeRhs> const& rhs)
  {
    return mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
          lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],
          lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8],
          lhs[3]*rhs[0]+lhs[4]*rhs[3]+lhs[5]*rhs[6],
          lhs[3]*rhs[1]+lhs[4]*rhs[4]+lhs[5]*rhs[7],
          lhs[3]*rhs[2]+lhs[4]*rhs[5]+lhs[5]*rhs[8],
          lhs[6]*rhs[0]+lhs[7]*rhs[3]+lhs[8]*rhs[6],
          lhs[6]*rhs[1]+lhs[7]*rhs[4]+lhs[8]*rhs[7],
          lhs[6]*rhs[2]+lhs[7]*rhs[5]+lhs[8]*rhs[8]);
  }

  //! lhs.transpose() * rhs
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  transpose_mul(
    mat3<NumTypeLhs> const& lhs,
    mat3<NumTypeRhs> const& rhs)
  {
    return mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[3]*rhs[3]+lhs[6]*rhs[6],
          lhs[0]*rhs[1]+lhs[3]*rhs[4]+lhs[6]*rhs[7],
          lhs[0]*rhs[2]+lhs[3]*rhs[5]+lhs[6]*rhs[8],
          lhs[1]*rhs[0]+lhs[4]*rhs[3]+lhs[7]*rhs[6],
          lhs[1]*rhs[1]+lhs[4]*rhs[4]+lhs[7]*rhs[7],
          lhs[1]*rhs[2]+lhs[4]*rhs[5]+lhs[7]*rhs[8],
          lhs[2]*rhs[0]+lhs[5]*rhs[3]+lhs[8]*rhs[6],
          lhs[2]*rhs[1]+lhs[5]*rhs[4]+lhs[8]*rhs[7],
          lhs[2]*rhs[2]+lhs[5]*rhs[5]+lhs[8]*rhs[8]);
  }

  //! lhs * rhs.transpose()
  template <typename NumTypeLhs, typename NumTypeRhs>
  inline
  mat3<
    typename af::binary_operator_traits<
      NumTypeLhs, NumTypeRhs>::arithmetic>
  mul_transpose(
    mat3<NumTypeLhs> const& lhs,
    mat3<NumTypeRhs> const& rhs)
  {
    return mat3<
      typename af::binary_operator_traits<
        NumTypeLhs, NumTypeRhs>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2],
          lhs[0]*rhs[3]+lhs[1]*rhs[4]+lhs[2]*rhs[5],
          lhs[0]*rhs[6]+lhs[1]*rhs[7]+lhs[2]*rhs[8],
          lhs[3]*rhs[0]+lhs[4]*rhs[1]+lhs[5]*rhs[2],
          lhs[3]*rhs[3]+lhs[4]*rhs[4]+lhs[5]*rhs[5],
          lhs[3]*rhs[6]+lhs[4]*rhs[7]+lhs[5]*rhs[8],
          lhs[6]*rhs[0]+lhs[7]*rhs[1]+lhs[8]*rhs[2],
          lhs[6]*rhs[3]+lhs[7]*rhs[4]+lhs[8]*rhs[5],
          lhs[6]*rhs[6]+lhs[7]*rhs[7]+lhs[8]*rhs[8]);
  }

  //! Matrix * vector product.
  template <typename NumTypeMatrix,
            typename NumTypeVector>
  inline
  vec3<
    typename af::binary_operator_traits<
      NumTypeMatrix, NumTypeVector>::arithmetic>
  operator*(
    mat3<NumTypeMatrix> const& lhs,
    af::tiny_plain<NumTypeVector,3> const& rhs)
  {
    return vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2],
          lhs[3]*rhs[0]+lhs[4]*rhs[1]+lhs[5]*rhs[2],
          lhs[6]*rhs[0]+lhs[7]*rhs[1]+lhs[8]*rhs[2]);
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
    mat3<NumTypeMatrix> const& rhs)
  {
    return vec3<
      typename af::binary_operator_traits<
        NumTypeMatrix, NumTypeVector>::arithmetic>(
          lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
          lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],
          lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8]);
  }

  //! Element-wise multiplication.
  template <typename NumType>
  inline
  mat3<NumType>
  operator*(
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
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
    mat3<NumType> const& lhs,
    NumType const& rhs)
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
    NumType const& lhs,
    mat3<NumType> const& rhs)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = lhs    % rhs[i];
    }
    return result;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator+=(
    mat3<NumType>& lhs,
    mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator+=(
    mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator-=(
    mat3<NumType>& lhs,
    mat3<NumType> const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator-=(
    mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place multiplication.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator*=(
    mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place division.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator/=(
    mat3<NumType>& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<9;i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place modulus operation.
  template <typename NumType>
  inline
  mat3<NumType>&
  operator%=(
    mat3<NumType>& lhs,
    NumType const& rhs)
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
    mat3<NumType> const& v)
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
    mat3<NumType> const& v)
  {
    mat3<NumType> result;
    for(std::size_t i=0;i<9;i++) {
      result[i] = +v[i];
    }
    return result;
  }

} // namespace scitbx

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename NumType>
  struct has_trivial_destructor<scitbx::mat3<NumType> > {
    static const bool value = ::boost::has_trivial_destructor<NumType>::value;
  };
}

#endif

#endif // SCITBX_MAT3_H
