/*
 * parallax_correction.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_PARALLAX_CORRECTION_H
#define DXTBX_MODEL_PARALLAX_CORRECTION_H

#include <scitbx/vec2.h>
#include <cmath>

namespace dxtbx { namespace model {

  using std::sqrt;
  using scitbx::vec2;

  /**
   * Function to perform a parallax correction on a given coordinate. Find the
   * xy coordinate by projecting the vector at the given coordinate la (the
   * attenuation length) into the detector.
   * @param d The distance from the detector to the source along the normal
   * @param la The attenuation length
   * @param xy0 The detector (mm) coordinate of the origin at the normal
   * @param xy The coordinate to correct
   */
  inline
  vec2<double> parallax_correction(double d, double la, vec2<double> xy0,
      vec2<double> xy) {
    vec2<double> xyp = xy - xy0;
    return xy + la * xyp / sqrt(d * d + xyp.length_sq());
  }

  /**
   * Function to perform the inverse parallax correction.
   *
   * N.B This function is an approximation. The real solution is to the
   * equation:
   *
   *    xy' = xy + xy * l / sqrt(h^2 + xy^2)
   *
   * Solving this equation for xy in terms of xy' gives an extremely complicated
   * solution. My maths wan't good enough to get a simpler solution. Might have
   * another look if the approximation becomes an issue but with reasonable
   * geometry, the approximation is correction to more than 3 decimal places.
   *
   * The actual equation used is xy = xy' - xy' * l / sqrt(h^2 + xy'^2)
   *
   * @param d The distance from the detector to the source along the normal
   * @param la The attenuation length
   * @param xy0 The detector (mm) coordinate of the origin at the normal
   * @param xy The coordinate to correct
   */
  inline
  vec2<double> parallax_correction_inv(double d, double la, vec2<double> xy0,
      vec2<double> xy) {
    vec2<double> xyp = xy - xy0;
    return xy - la * xyp / sqrt(d * d + xyp.length_sq());
  }


}} // namespace dxtbx::model

#endif /* DXTBX_MODEL_PARALLAX_CORRECTION_H */
