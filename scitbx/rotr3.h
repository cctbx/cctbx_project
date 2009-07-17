#ifndef SCITBX_ROTR3_H
#define SCITBX_ROTR3_H

#include <scitbx/mat3.h>

namespace scitbx {

  template <typename FloatType>
  struct rotr3
  {
    mat3<FloatType> r;
    vec3<FloatType> t;

    rotr3() {}

    rotr3(
      mat3<FloatType> const& r_,
      vec3<FloatType> const& t_)
    :
      r(r_),
      t(t_)
    {}

    static
    rotr3
    identity()
    {
      return rotr3(mat3<FloatType>(1,1,1), vec3<FloatType>(0,0,0));
    }
  };

} // namespace scitbx

#endif // GUARD
