#ifndef SCITBX_MATRIX_TENSOR_RANK_2_H
#define SCITBX_MATRIX_TENSOR_RANK_2_H

#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace scitbx { namespace matrix {

//! Computations involving tensors of rank 2.
namespace tensor_rank_2 {

  /*! \brief Transformation of gradients (first derivatives)
      w.r.t. the elements of a symmetric tensor.
   */
  /*! @param a is a linear transformation.
      @param g are the gradients.
      <p>
      Mathematica script used to determine the transformation law:<pre>
        g={{g0,g3,g4},{g3,g1,g5},{g4,g5,g2}}
        a={{a0,a1,a2},{a3,a4,a5},{a6,a7,a8}}
        agat = a.g.Transpose[a]
        D[agat, g0]
        D[agat, g1]
        D[agat, g2]
        D[agat, g3]
        D[agat, g4]
        D[agat, g5]</pre>
      Output:<pre>
            2                           2                           2
        {{a0 , a0 a3, a0 a6}, {a0 a3, a3 , a3 a6}, {a0 a6, a3 a6, a6 }}

            2                           2                           2
        {{a1 , a1 a4, a1 a7}, {a1 a4, a4 , a4 a7}, {a1 a7, a4 a7, a7 }}

            2                           2                           2
        {{a2 , a2 a5, a2 a8}, {a2 a5, a5 , a5 a8}, {a2 a8, a5 a8, a8 }}

        {{2 a0 a1, a1 a3 + a0 a4, a1 a6 + a0 a7},
         {a1 a3 + a0 a4, 2 a3 a4, a4 a6 + a3 a7},
         {a1 a6 + a0 a7, a4 a6 + a3 a7, 2 a6 a7}}

        {{2 a0 a2, a2 a3 + a0 a5, a2 a6 + a0 a8},
         {a2 a3 + a0 a5, 2 a3 a5, a5 a6 + a3 a8},
         {a2 a6 + a0 a8, a5 a6 + a3 a8, 2 a6 a8}}

        {{2 a1 a2, a2 a4 + a1 a5, a2 a7 + a1 a8},
         {a2 a4 + a1 a5, 2 a4 a5, a5 a7 + a4 a8},
         {a2 a7 + a1 a8, a5 a7 + a4 a8, 2 a7 a8}}
      </pre>
      The transformed gradients are obtained by the dot products of
      the upper triangles of the derivative matrices with the
      original gradients g.
   */
  template <typename NumTypeA, typename NumTypeG>
  inline sym_mat3<NumTypeG>
  gradient_transform(
    mat3<NumTypeA> const& a,
    sym_mat3<NumTypeG> const& g)
  {
    return sym_mat3<NumTypeG>(
       a[0]*a[0]*g[0]
      +a[3]*a[3]*g[1]
      +a[6]*a[6]*g[2]
      +a[0]*a[3]*g[3]
      +a[0]*a[6]*g[4]
      +a[3]*a[6]*g[5],
       a[1]*a[1]*g[0]
      +a[4]*a[4]*g[1]
      +a[7]*a[7]*g[2]
      +a[1]*a[4]*g[3]
      +a[1]*a[7]*g[4]
      +a[4]*a[7]*g[5],
       a[2]*a[2]*g[0]
      +a[5]*a[5]*g[1]
      +a[8]*a[8]*g[2]
      +a[2]*a[5]*g[3]
      +a[2]*a[8]*g[4]
      +a[5]*a[8]*g[5],
       a[0]*a[1]*2*g[0]
      +a[3]*a[4]*2*g[1]
      +a[6]*a[7]*2*g[2]
      +(a[1]*a[3]+a[0]*a[4])*g[3]
      +(a[1]*a[6]+a[0]*a[7])*g[4]
      +(a[4]*a[6]+a[3]*a[7])*g[5],
       a[0]*a[2]*2*g[0]
      +a[3]*a[5]*2*g[1]
      +a[6]*a[8]*2*g[2]
      +(a[2]*a[3]+a[0]*a[5])*g[3]
      +(a[2]*a[6]+a[0]*a[8])*g[4]
      +(a[5]*a[6]+a[3]*a[8])*g[5],
       a[1]*a[2]*2*g[0]
      +a[4]*a[5]*2*g[1]
      +a[7]*a[8]*2*g[2]
      +(a[2]*a[4]+a[1]*a[5])*g[3]
      +(a[2]*a[7]+a[1]*a[8])*g[4]
      +(a[5]*a[7]+a[4]*a[8])*g[5]);
  }

  //! See overload.
  template <typename NumTypeA>
  void
  gradient_transform_matrix(
    NumTypeA* result,
    NumTypeA const* a)
  {
    *result++ = a[0]*a[0];
    *result++ = a[3]*a[3];
    *result++ = a[6]*a[6];
    *result++ = a[0]*a[3];
    *result++ = a[0]*a[6];
    *result++ = a[3]*a[6];
    *result++ = a[1]*a[1];
    *result++ = a[4]*a[4];
    *result++ = a[7]*a[7];
    *result++ = a[1]*a[4];
    *result++ = a[1]*a[7];
    *result++ = a[4]*a[7];
    *result++ = a[2]*a[2];
    *result++ = a[5]*a[5];
    *result++ = a[8]*a[8];
    *result++ = a[2]*a[5];
    *result++ = a[2]*a[8];
    *result++ = a[5]*a[8];
    *result++ = a[0]*a[1]*2;
    *result++ = a[3]*a[4]*2;
    *result++ = a[6]*a[7]*2;
    *result++ = a[1]*a[3]+a[0]*a[4];
    *result++ = a[1]*a[6]+a[0]*a[7];
    *result++ = a[4]*a[6]+a[3]*a[7];
    *result++ = a[0]*a[2]*2;
    *result++ = a[3]*a[5]*2;
    *result++ = a[6]*a[8]*2;
    *result++ = a[2]*a[3]+a[0]*a[5];
    *result++ = a[2]*a[6]+a[0]*a[8];
    *result++ = a[5]*a[6]+a[3]*a[8];
    *result++ = a[1]*a[2]*2;
    *result++ = a[4]*a[5]*2;
    *result++ = a[7]*a[8]*2;
    *result++ = a[2]*a[4]+a[1]*a[5];
    *result++ = a[2]*a[7]+a[1]*a[8];
    *result   = a[5]*a[7]+a[4]*a[8];
  }

  /*! \brief Transformation matrix for gradients (first derivatives)
      w.r.t. the elements of a symmetric tensor.
   */
  /*! The result is a 6x6 matrix of the coefficients as used in
      gradient_transform(). The transformed gradients are obtained
      by multiplying this matrix with the vector of original
      gradients.
   */
  template <typename NumTypeA>
  af::versa<NumTypeA, af::c_grid<2> >
  gradient_transform_matrix(
    mat3<NumTypeA> const& a)
  {
    af::versa<NumTypeA, af::c_grid<2> > result(
      af::c_grid<2>(6,6), af::init_functor_null<NumTypeA>());
    gradient_transform_matrix(result.begin(), a.begin());
    return result;
  }

}}} // namespace scitbx::matrix::tensor_rank_2

#endif // SCITBX_MATRIX_TENSOR_RANK_2_H
