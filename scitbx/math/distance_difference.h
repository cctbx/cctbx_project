#ifndef SCITBX_MATH_DISTANCE_DIFFERENCE_H
#define SCITBX_MATH_DISTANCE_DIFFERENCE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/vec3.h>
#include <scitbx/error.h>

namespace scitbx { namespace math {

  template <typename FloatType>
  af::versa< FloatType, af::c_grid<2> > distance_difference_matrix (
    af::shared< vec3<FloatType> > sites1,
    af::shared< vec3<FloatType> > sites2)
  {
    SCITBX_ASSERT(sites1.size() == sites2.size());
    af::c_grid<2> ddm_grid(sites1.size(), sites1.size());
    af::versa< FloatType, af::c_grid<2> > ddm(ddm_grid, 0.0);
    for (unsigned i = 0; i < sites1.size(); i++) {
      for (unsigned j = 0; j < sites1.size(); j++) {
        FloatType d1 = (sites1[i] - sites1[j]).length();
        FloatType d2 = (sites2[i] - sites2[j]).length();
        ddm(i,j) = (d2 - d1);
      }
    }
    return ddm;
  }

}} // namespace scitbx::math

#endif
