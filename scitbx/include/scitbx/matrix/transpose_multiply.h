#ifndef SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H
#define SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H

#include <scitbx/mat3.h>

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
      for(unsigned j=0;j<3;j++) {
        for(unsigned k=0;k<3;k++) {
          result[j*3+k] += lhs[i][j] * rhs[i][k];
        }
      }
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_TRANSPOSE_MULTIPLY_H
