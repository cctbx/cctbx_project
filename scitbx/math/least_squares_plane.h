#ifndef SCITBX_MATH_LEAST_SQUARE_PLANE_H
#define SCITBX_MATH_LEAST_SQUARE_PLANE_H

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>

namespace scitbx { namespace math {

template<typename ScalarType=double>
class least_squares_plane
{
  public:
    typedef ScalarType scalar_type;
    typedef vec3<scalar_type> vector_type;
    typedef mat3<scalar_type> matrix_type;

    least_squares_plane(af::const_ref<vector_type> const &points)
    {
      float_type xx=0, yy=0, zz=0, xy=0, xz=0, yz=0;
      vector_type b(0,0,0);
      for(int i=0; i < points.size(); ++i) {
        vector_type const &p = points[i];
        b += p;
        float_type x=p[0], y=p[1], z=p[2];
        xx += x*x;
        yy += y*y;
        zz += z*z;
        xy += x*y;
        xz += x*z;
        yz += y*z;
      }
      sym_mat3<float_type> m(xx, yy, zz, xy, xz, yz);
      vector_type u = m.inverse()*b;
      d = 1/u.length();
      n = d*u;
    }

    vector_type const &normal() { return n; }

    scalar_type distance_to_origin() { return d; }

  private:
    vector_type n;
    scalar_type d;
};

}}

#endif // GUARD
