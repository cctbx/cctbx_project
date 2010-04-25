#ifndef SCITBX_MATH_WEIGHTED_LINEAR_CORRELATION_H
#define SCITBX_MATH_WEIGHTED_LINEAR_CORRELATION_H

#include <scitbx/array_family/ref.h>
#include <boost/optional.hpp>

namespace scitbx { namespace math {

  /// Variance-Covariance of a 2-dimensional distribution (x,y)
  template <typename FloatType>
  class weighted_covariance
  {
  public:
    typedef FloatType float_type;

    /// No points: use accumulate
    weighted_covariance()
     : sum_w(0), mean_x_(0), mean_y_(0), m_xx(0), m_xy(0), m_yy(0)
    {}

    /// Variance-covariance of the given sequence of weighted points
    /** This uses the corrected two-pass algorithm. */
    weighted_covariance(af::const_ref<float_type> const &x,
                        af::const_ref<float_type> const &y,
                        af::const_ref<float_type> const &w)
     : sum_w(0), mean_x_(0), mean_y_(0), m_xx(0), m_xy(0), m_yy(0)
    {
      SCITBX_ASSERT(x.size() == w.size());
      SCITBX_ASSERT(y.size() == w.size());
      int n = w.size();
      for (int i=0; i<n; ++i) {
        sum_w += w[i];
        mean_x_ += w[i]*x[i];
        mean_y_ += w[i]*y[i];
      }
      SCITBX_ASSERT(sum_w);
      mean_x_ /= sum_w;
      mean_y_ /= sum_w;
      float_type sum_w_delta_x = 0, sum_w_delta_y = 0;
      for (int i=0; i<n; ++i) {
        float_type delta_x = x[i] - mean_x_,
                   delta_y = y[i] - mean_y_;
        sum_w_delta_x += w[i]*delta_x;
        sum_w_delta_y += w[i]*delta_y;
        m_xx += w[i]*delta_x*delta_x;
        m_xy += w[i]*delta_x*delta_y;
        m_yy += w[i]*delta_y*delta_y;
      }
      m_xx -= sum_w_delta_x*sum_w_delta_x/sum_w;
      m_xy -= sum_w_delta_x*sum_w_delta_y/sum_w;
      m_yy -= sum_w_delta_y*sum_w_delta_y/sum_w;
    }

    /// Recompute the statistics for the sequence augmented by the given point
    /** This uses the famous updating formula usually credited to Knuth. */
    weighted_covariance &accumulate(float_type x, float_type y, float_type w=1)
    {
      sum_w += w;
      float_type w_over_sum_w = w/sum_w;
      float_type delta_x = x - mean_x_,
                 delta_y = y - mean_y_;
      mean_x_ += w_over_sum_w*delta_x;
      mean_y_ += w_over_sum_w*delta_y;
      float_type new_delta_x = x - mean_x_,
                 new_delta_y = y - mean_y_;
      m_xx += w*delta_x*new_delta_x;
      m_xy += w*delta_x*new_delta_y;
      m_yy += w*delta_y*new_delta_y;
      return *this;
    }

    float_type sum_weights() const { return sum_w; }

    float_type mean_x() const { return mean_x_; }

    float_type mean_y() const { return mean_y_; }

    float_type variance_x() const {
      SCITBX_ASSERT(sum_w);
      return m_xx/sum_w;
    }

    float_type variance_y() const {
      SCITBX_ASSERT(sum_w);
      return m_yy/sum_w;
    }

    float_type covariance_xy() const {
      SCITBX_ASSERT(sum_w);
      return m_xy/sum_w;
    }

    /// Correlation between x and y
    /** Unitialised if it is ill-defined */
    boost::optional<float_type> correlation() const {
      boost::optional<float_type> result;
      if      (m_xx != 0 && m_yy != 0) result = m_xy/std::sqrt(m_xx * m_yy);
      else if (m_xy == 0)              result = 1;
      return result;
    }

  private:
    float_type sum_w, mean_x_, mean_y_, m_xx, m_xy, m_yy;
  };

}}

#endif
