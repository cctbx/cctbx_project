#ifndef SCITBX_MATH_LEAST_SQUARE_PLANE_H
#define SCITBX_MATH_LEAST_SQUARE_PLANE_H

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>

namespace scitbx { namespace math {

template<typename FloatType=double>
class least_square_plane
{
  public:
    least_square_plane(af::const_ref<vec3<FloatType> > const &points,
                       vec3<FloatType> const &origin)
    {
      FloatType xx=0, yy=0, zz=0, xy=0, xz=0, yz=0;
      vec3<FloatType> b(0,0,0);
      for(int i=0; i < points.size(); ++i) {
        vec3<FloatType> p = points[i] - origin;
        b += p;
        FloatType x=p[0], y=p[1], z=p[2];
        xx += x*x;
        yy += y*y;
        zz += z*z;
        xy += x*y;
        xz += x*z;
        yz += y*z;
      }
      sym_mat3<FloatType> m(xx, yy, zz, xy, xz, yz);
      vec3<FloatType> u = m.inverse()*b;
      d = 1/u.length();
      n = d*u;
    }

    vec3<FloatType> const &normal() { return n; }

    FloatType distance_to_origin() { return d; }

  private:
    vec3<FloatType> n;
    FloatType d;
};

}}

#endif // GUARD
