/*
* model_helpers.h
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_MODEL_HELPERS_H
#define DXTBX_MODEL_MODEL_HELPERS_H

#include <cmath>
#include <scitbx/vec3.h>

namespace dxtbx { namespace model {

  using scitbx::vec3;

  inline
  double angle_safe(const vec3<double> &a, const vec3<double> &b) {
    double den = a.length() * b.length();
    if (den <= 0) return 0.0;
    double c = (a * b) / den;
    if      (c < -1) c = -1;
    else if (c >  1) c =  1;
    return std::acos(c);
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_MODEL_HELPERS_H
