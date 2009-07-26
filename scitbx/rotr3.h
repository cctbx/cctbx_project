#ifndef SCITBX_ROTR3_H
#define SCITBX_ROTR3_H

#include <scitbx/mat3.h>

namespace scitbx {

  template <typename FloatType>
  struct rotr3
  {
    typedef FloatType ft;

    mat3<ft> r;
    vec3<ft> t;

    rotr3() {}

    rotr3(
      mat3<ft> const& r_,
      vec3<ft> const& t_)
    :
      r(r_),
      t(t_)
    {}

    static
    rotr3
    identity()
    {
      return rotr3(mat3<ft>(1,1,1), vec3<ft>(0,0,0));
    }

    vec3<ft>
    operator*(
      vec3<ft> const& rhs) const
    {
      return r * rhs + t;
    }

    rotr3
    operator*(
      rotr3 const& rhs) const
    {
      return rotr3(r * rhs.r, r * rhs.t + t);
    }
  };

} // namespace scitbx

#endif // GUARD
