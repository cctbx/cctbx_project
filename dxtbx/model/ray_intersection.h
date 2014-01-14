/*
 * ray_intersection.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_RAY_INTERSECTION_H
#define DXTBX_MODEL_RAY_INTERSECTION_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;

  /** Get the coordinate of a ray intersecting with the detector */
  inline
  vec2<double> plane_ray_intersection(mat3<double> D, vec3<double> s1) {
    vec3 <double> v = D * s1;
    DXTBX_ASSERT(v[2] > 0);
    return vec2<double>(v[0] / v[2], v[1] / v[2]);
  }

  /** Get the coordinate of a ray intersecting with the detector */
  inline
  vec2<double> bidirectional_plane_ray_intersection(
      mat3<double> D, vec3<double> s1) {
    vec3 <double> v = D * s1;
    DXTBX_ASSERT(v[2] != 0);
    return vec2<double>(v[0] / v[2], v[1] / v[2]);
  }

  /** Get world coordinate of plane xy */
  inline
  vec3<double> plane_world_coordinate(mat3<double> d, vec2<double> xy) {
    return d * vec3<double>(xy[0], xy[1], 1.0);
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_RAY_INTERSECTION_H

