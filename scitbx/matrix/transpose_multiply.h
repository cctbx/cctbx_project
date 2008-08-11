#ifndef SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H
#define SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H

#include <scitbx/mat3.h>
#include <scitbx/mat2.h>

namespace scitbx { namespace matrix {

  template <typename NumType>
  mat3<NumType>
  transpose_multiply(
    af::const_ref<vec3<NumType> > const& lhs,
    af::const_ref<vec3<NumType> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    mat3<NumType> result(static_cast<NumType>(0));
    for(std::size_t i=0;i<lhs.size();i++) {
      NumType* r_jk = result.begin();
      for(unsigned j=0;j<3;j++) {
        NumType lhs_ij = lhs[i][j];
        const NumType* rhs_ik = rhs[i].begin();
        *(r_jk++) += lhs_ij * *(rhs_ik++);
        *(r_jk++) += lhs_ij * *(rhs_ik++);
        *(r_jk++) += lhs_ij * *(rhs_ik  );
      }
    }
    return result;
  }

  template <typename NumType>
  mat2<NumType>
  transpose_multiply(
    af::const_ref<vec2<NumType> > const& lhs,
    af::const_ref<vec2<NumType> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    mat2<NumType> result(static_cast<NumType>(0));
    for(std::size_t i=0;i<lhs.size();i++) {
      NumType* r_jk = result.begin();
      for(unsigned j=0;j<2;j++) {
        NumType lhs_ij = lhs[i][j];
        const NumType* rhs_ik = rhs[i].begin();
        *(r_jk++) += lhs_ij * *(rhs_ik++);
        *(r_jk++) += lhs_ij * *(rhs_ik  );
      }
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H
