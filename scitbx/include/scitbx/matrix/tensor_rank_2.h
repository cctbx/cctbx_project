#ifndef SCITBX_MATRIX_TENSOR_RANK_2_H
#define SCITBX_MATRIX_TENSOR_RANK_2_H

#include <scitbx/sym_mat3.h>

namespace scitbx { namespace matrix {

//! Computations involving tensors of rank 2.
namespace tensor_rank_2 {

  /*! \brief Transformation of gradients (first derivatives)
      w.r.t. the elements of a symmetric tensor.
   */
  /*! @param a is the transformation matrix.
      @param g are the gradients.
      <p>
      Mathematica script used to determine the transformation law:<pre>
        SetOptions["stdout", PageWidth -> 50]
        us={{s00,s01,s02},{s01,s11,s12},{s02,s12,s22}}
        hs={hs0,hs1,hs2}
        fs=Exp[cb hs.us.hs]
        gs={{D[fs,s00],D[fs,s01],D[fs,s02]},
            {D[fs,s01],D[fs,s11],D[fs,s12]},
            {D[fs,s02],D[fs,s12],D[fs,s22]}}/fs
        uc={{c00,c01,c02},{c01,c11,c12},{c02,c12,c22}}
        hc={hc0,hc1,hc2}
        fc=Exp[cb hc.uc.hc]
        gc={{D[fc,c00],D[fc,c01],D[fc,c02]},
            {D[fc,c01],D[fc,c11],D[fc,c12]},
            {D[fc,c02],D[fc,c12],D[fc,c22]}}/fc
        a={{a00,a01,a02},
           {a10,a11,a12},
           {a20,a21,a22}}
        hc=Transpose[a].hs
        hc0=hc[[1]]
        hc1=hc[[2]]
        hc2=hc[[3]]
        FullSimplify[gs]
        gcs=Expand[Expand[FullSimplify[gc]]
          /. cb hs0^2 -> g00
          /. cb hs1^2 -> g11
          /. cb hs2^2 -> g22
          /. cb hs0 hs1 -> g01/2
          /. cb hs0 hs2 -> g02/2
          /. cb hs1 hs2 -> g12/2]
        InputForm[gcs[[1,1]]]
        InputForm[gcs[[2,2]]]
        InputForm[gcs[[3,3]]]
        InputForm[gcs[[1,2]]]
        InputForm[gcs[[1,3]]]
        InputForm[gcs[[2,3]]]</pre>
   */
  template <typename NumTypeA, typename NumTypeG>
  inline sym_mat3<NumTypeG>
  gradient_transform(
    mat3<NumTypeA> const& a,
    sym_mat3<NumTypeG> const& g)
  {
    return sym_mat3<NumTypeG>(
        a[0]*a[0]*g[0] + a[0]*a[3]*g[3] + a[0]*a[6]*g[4]
      + a[3]*a[3]*g[1] + a[3]*a[6]*g[5] + a[6]*a[6]*g[2],
        a[1]*a[1]*g[0] + a[1]*a[4]*g[3] + a[1]*a[7]*g[4]
      + a[4]*a[4]*g[1] + a[4]*a[7]*g[5] + a[7]*a[7]*g[2],
        a[2]*a[2]*g[0] + a[2]*a[5]*g[3] + a[2]*a[8]*g[4]
      + a[5]*a[5]*g[1] + a[5]*a[8]*g[5] + a[8]*a[8]*g[2],
        2*a[0]*a[1]*g[0] + a[1]*a[3]*g[3] +   a[0]*a[4]*g[3]
      +   a[1]*a[6]*g[4] + a[0]*a[7]*g[4] + 2*a[3]*a[4]*g[1]
      +   a[4]*a[6]*g[5] + a[3]*a[7]*g[5] + 2*a[6]*a[7]*g[2],
        2*a[0]*a[2]*g[0] + a[2]*a[3]*g[3] +   a[0]*a[5]*g[3]
      +   a[2]*a[6]*g[4] + a[0]*a[8]*g[4] + 2*a[3]*a[5]*g[1]
      +   a[5]*a[6]*g[5] + a[3]*a[8]*g[5] + 2*a[6]*a[8]*g[2],
        2*a[1]*a[2]*g[0] + a[2]*a[4]*g[3] +   a[1]*a[5]*g[3]
      +   a[2]*a[7]*g[4] + a[1]*a[8]*g[4] + 2*a[4]*a[5]*g[1]
      +   a[5]*a[7]*g[5] + a[4]*a[8]*g[5] + 2*a[7]*a[8]*g[2]);
  }

  //! Refactored version of gradient_transform() above.
  template <typename NumTypeA>
  struct gradient_average_cache
  {
    //! Array of 45 accumulated coefficients.
    /*! Not available in Python.
     */
    af::tiny<NumTypeA, 45> c;

    //! Initialization of all coefficients with zero.
    gradient_average_cache() { c.fill(0); }

    //! Accumulation of coefficients.
    void
    accumulate(mat3<NumTypeA> const& a)
    {
      c[ 0]+=a[0]*a[0]; c[ 1]+=a[0]*a[3]; c[ 2]+=a[0]*a[6];
      c[ 3]+=a[3]*a[3]; c[ 4]+=a[3]*a[6]; c[ 5]+=a[6]*a[6];
      c[ 6]+=a[1]*a[1]; c[ 7]+=a[1]*a[4]; c[ 8]+=a[1]*a[7];
      c[ 9]+=a[4]*a[4]; c[10]+=a[4]*a[7]; c[11]+=a[7]*a[7];
      c[12]+=a[2]*a[2]; c[13]+=a[2]*a[5]; c[14]+=a[2]*a[8];
      c[15]+=a[5]*a[5]; c[16]+=a[5]*a[8]; c[17]+=a[8]*a[8];
      c[18]+=2*a[0]*a[1]; c[19]+=a[1]*a[3]; c[20]+=  a[0]*a[4];
      c[21]+=  a[1]*a[6]; c[22]+=a[0]*a[7]; c[23]+=2*a[3]*a[4];
      c[24]+=  a[4]*a[6]; c[25]+=a[3]*a[7]; c[26]+=2*a[6]*a[7];
      c[27]+=2*a[0]*a[2]; c[28]+=a[2]*a[3]; c[29]+=  a[0]*a[5];
      c[30]+=  a[2]*a[6]; c[31]+=a[0]*a[8]; c[32]+=2*a[3]*a[5];
      c[33]+=  a[5]*a[6]; c[34]+=a[3]*a[8]; c[35]+=2*a[6]*a[8];
      c[36]+=2*a[1]*a[2]; c[37]+=a[2]*a[4]; c[38]+=  a[1]*a[5];
      c[39]+=  a[2]*a[7]; c[40]+=a[1]*a[8]; c[41]+=2*a[4]*a[5];
      c[42]+=  a[5]*a[7]; c[43]+=a[4]*a[8]; c[44]+=2*a[7]*a[8];
    }

    //! Application of accumulated coefficients.
    template <typename NumTypeG>
    sym_mat3<NumTypeG>
    average(sym_mat3<NumTypeG> const& g, NumTypeG const& denominator) const
    {
      return sym_mat3<NumTypeG>(
        (  c[ 0]*g[0] + c[ 1]*g[3] + c[ 2]*g[4]
         + c[ 3]*g[1] + c[ 4]*g[5] + c[ 5]*g[2]) / denominator,
        (  c[ 6]*g[0] + c[ 7]*g[3] + c[ 8]*g[4]
         + c[ 9]*g[1] + c[10]*g[5] + c[11]*g[2]) / denominator,
        (  c[12]*g[0] + c[13]*g[3] + c[14]*g[4]
         + c[15]*g[1] + c[16]*g[5] + c[17]*g[2]) / denominator,
        (  c[18]*g[0] + c[19]*g[3] + c[20]*g[3]
         + c[21]*g[4] + c[22]*g[4] + c[23]*g[1]
         + c[24]*g[5] + c[25]*g[5] + c[26]*g[2]) / denominator,
        (  c[27]*g[0] + c[28]*g[3] + c[29]*g[3]
         + c[30]*g[4] + c[31]*g[4] + c[32]*g[1]
         + c[33]*g[5] + c[34]*g[5] + c[35]*g[2]) / denominator,
        (  c[36]*g[0] + c[37]*g[3] + c[38]*g[3]
         + c[39]*g[4] + c[40]*g[4] + c[41]*g[1]
         + c[42]*g[5] + c[43]*g[5] + c[44]*g[2]) / denominator);
    }
  };

}}} // namespace scitbx::matrix::tensor_rank_2

#endif // SCITBX_MATRIX_TENSOR_RANK_2_H
