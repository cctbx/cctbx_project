#ifndef SCITBX_VEC3_ORTHONORMALISE_H
#define SCITBX_VEC3_ORTHONORMALISE_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/error.h>
#include <scitbx/math/modulo.h>

namespace scitbx { namespace math {

  /**
   @defgroup orthonormal_basis Orthonormal basis construction in dimension 3.
   i.e. three unit 3-vector e0, e1, e2 orthogonal to each others.
   */
  /*@{*/

  /** Construct a basis by orthonormalizing the given pair of vectors.
      The construction is such that:
        - e0 is parallel to and in the same direction as v0;
        - e1 is in the plane spanned by v0 and v1;
        - e1 and v1 are on the same side of v0.
      Implementation note: although this is the famous Gram-Schmidt procedure,
      it is not implemented in the classic manner with projections. We take
      advantage of dimension 3 with cross-products.
   */
  template <typename T>
  af::tiny<vec3<T>, 3>
  orthonormal_basis(vec3<T> const &v0, vec3<T> const &v1,
                    bool right_handed=true)
  {
    af::tiny<vec3<T>, 3> e;
    e[0] = v0.normalize();
    e[2] = e[0].cross(v1);
    T l2 = e[2].length();
    SCITBX_ASSERT(l2 > 0)(l2);
    e[2] /= l2;
    e[1] = e[2].cross(e[0]);
    if(!right_handed) e[2] = -e[2];
    return e;
  }

  /// Construct a basis whose i0-th and i1-th vectors are built from v0 and v1.
  /** The procedure is the same as the other overload, with e0 and e1
      replaced by respectively the i0-th and i1-th basis vector.
   */
  template <typename T>
  af::tiny<vec3<T>, 3>
  orthonormal_basis(vec3<T> const &v0, int i0,
                    vec3<T> const &v1, int i1,
                    bool right_handed=true)
  {
    SCITBX_ASSERT(   i0 != i1
                  && 0 <= i0 && i0 < 3
                  && 0 <= i1 && i1 < 3)(i0)(i1);
    af::tiny<vec3<T>, 3> f = orthonormal_basis(v0, v1, right_handed);
    int i2 = 3 - i0 - i1;
    af::tiny<vec3<T>, 3> e;
    e[i0] = f[0];
    e[i1] = f[1];
    e[i2] = f[2];
    /** if i0 --> i1 is a move by -1, the orientation being given by
     1 -> 2 -> 3
     ^         |
     |_________|
     */
    if (mod_short(i1 - i0, 3) == -1) e[i2] = -e[i2];
    return e;
  }

  /*@}*/

}}

#endif // GUARD
