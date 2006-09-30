#ifndef CCTBX_MILLER_MATH_H
#define CCTBX_MILLER_MATH_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace miller {

  /*! sum((w/epsilon)*data) / sum(w)
      anomalous_flag = true: w = 1
      anomalous_flag = false: w = 1 for centric reflections
                              w = 2 for acentric reflections
      <p>
      See also: sgtbx::space_group::epsilon()
   */
  template <typename FloatType>
  FloatType
  statistical_mean(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& miller_indices,
    af::const_ref<FloatType> const& data)
  {
    FloatType sum_numerator = 0;
    FloatType sum_denominator = 0;
    FloatType w = 1;
    bool fixed_w = (anomalous_flag || space_group.is_centric());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      FloatType e = space_group.epsilon(miller_indices[i]);
      if (!fixed_w) {
        w = 1;
        if (!space_group.is_centric(miller_indices[i])) w = 2;
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
