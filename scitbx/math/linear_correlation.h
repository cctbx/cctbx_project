#ifndef SCITBX_MATH_LINEAR_CORRELATION_H
#define SCITBX_MATH_LINEAR_CORRELATION_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/misc_functions.h>

namespace scitbx { namespace math {

  template <typename FloatType = double>
  class linear_correlation
  {
    public:
      typedef FloatType float_type;

      template <typename OtherFloatType>
      linear_correlation(
        af::const_ref<OtherFloatType> const& x,
        af::const_ref<OtherFloatType> const& y,
        FloatType const& epsilon=1.e-15,
        bool const& subtract_mean=true)
      :
        is_well_defined_(false),
        n_(x.size()),
        mean_x_(0),
        mean_y_(0),
        numerator_(0),
        sum_denominator_x_(0),
        sum_denominator_y_(0),
        denominator_(0),
        coefficient_(0)
      {
        SCITBX_ASSERT(x.size() == y.size());
        if (n_ != 0) {
          if(subtract_mean) {
            for(std::size_t i=0;i<n_;i++) mean_x_ += x[i];
            for(std::size_t i=0;i<n_;i++) mean_y_ += y[i];
            mean_x_ /= n_;
            mean_y_ /= n_;
          }
          for(std::size_t i=0;i<n_;i++) {
            FloatType delta_x = x[i] - mean_x_;
            FloatType delta_y = y[i] - mean_y_;
            numerator_ += delta_x * delta_y;
            sum_denominator_x_ += delta_x * delta_x;
            sum_denominator_y_ += delta_y * delta_y;
          }
          denominator_ = std::sqrt(sum_denominator_x_ * sum_denominator_y_);
          if (numerator_ == 0 && denominator_ == 0) {
            coefficient_ = 1;
            is_well_defined_ = true;
          }
          else if (denominator_ > fn::absolute(numerator_ * epsilon)) {
            coefficient_ = numerator_ / denominator_;
            is_well_defined_ = true;
          }
        }
      }

      bool
      is_well_defined() const { return is_well_defined_; }

      std::size_t
      n() const { return n_; }

      FloatType
      mean_x() const { return mean_x_; }

      FloatType
      mean_y() const { return mean_y_; }

      FloatType
      numerator() const { return numerator_; }

      FloatType
      sum_denominator_x() const { return sum_denominator_x_; }

      FloatType
      sum_denominator_y() const { return sum_denominator_y_; }

      FloatType
      denominator() const { return denominator_; }

      FloatType
      coefficient() const { return coefficient_; }

    private:
      bool is_well_defined_;
      std::size_t n_;
      FloatType mean_x_;
      FloatType mean_y_;
      FloatType numerator_;
      FloatType sum_denominator_x_;
      FloatType sum_denominator_y_;
      FloatType denominator_;
      FloatType coefficient_;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_LINEAR_CORRELATION_H
