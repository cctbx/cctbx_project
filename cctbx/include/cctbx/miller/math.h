// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_MATH_H
#define CCTBX_MILLER_MATH_H

#include <cctbx/sgtbx/groups.h>

namespace cctbx { namespace miller {

  /*! sum((w/epsilon)*f) / sum(w)
      friedel_flag = false: w = 1
      friedel_flag = true: w = 1 for centric reflections
                           w = 2 for acentric reflections
   */
  template <typename FloatType>
  FloatType
  statistical_mean(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    af::shared<Index> miller_indices,
    af::shared<FloatType> data)
  {
    FloatType sum_numerator = 0;
    FloatType sum_denominator = 0;
    FloatType w = 1;
    bool fixed_w = (!friedel_flag || SgOps.isCentric());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      FloatType e = SgOps.epsilon(miller_indices[i]);
      if (!fixed_w) {
        w = 1;
        if (!SgOps.isCentric(miller_indices[i])) w = 2;
        sum_denominator += w;
      }
      sum_numerator += w / e * data[i];
    }
    if (fixed_w) sum_denominator = miller_indices.size();
    if (sum_denominator == 0) return 0;
    return sum_numerator / sum_denominator;
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATH_H
