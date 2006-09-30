#ifndef CCTBX_MATH_ARRAY_UTILS_H
#define CCTBX_MATH_ARRAY_UTILS_H

#include <cctbx/fixes/cmath>
#include <cctbx/array_family/ref_reductions.h>

namespace cctbx { namespace math {

  template <typename FloatType>
  class array_statistics
  {
    public:
      array_statistics() {}
      array_statistics(af::const_ref<FloatType> const& a)
      {
        m_min = af::min(a);
        m_max = af::max(a);
        m_mean = af::mean(a);
        m_mean2 = af::mean_sq(a);
        m_sigma = m_mean2 - m_mean * m_mean;
        if (m_sigma < FloatType(0)) m_sigma = 0;
        m_sigma = std::sqrt(m_sigma);
      }
      FloatType const& min() const { return m_min; }
      FloatType const& max() const { return m_max; }
      FloatType const& mean() const { return m_mean; }
      FloatType const& mean2() const { return m_mean2; }
      FloatType const& sigma() const { return m_sigma; }
    protected:
      FloatType m_min;
      FloatType m_max;
      FloatType m_mean;
      FloatType m_mean2;
      FloatType m_sigma;
  };

}} // namespace cctbx::math

#endif // CCTBX_MATH_ARRAY_UTILS_H
