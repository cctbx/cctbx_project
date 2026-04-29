#ifndef CCTBX_MATH_COS_SIN_TABLE_H
#define CCTBX_MATH_COS_SIN_TABLE_H

#include <cctbx/error.h>
#include <scitbx/constants.h>
#include <boost/shared_array.hpp>
#include <complex>
#include <cstdio>

namespace cctbx { namespace math {

  template <typename FloatType=double>
  class cos_sin_exact
  {
    public:
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      static
      complex_type
      get(float_type unary_angle)
      {
        float_type x = scitbx::constants::two_pi * unary_angle;
        return complex_type(std::cos(x), std::sin(x));
      }

      complex_type operator()(float_type unary_angle) const {
        return get(unary_angle);
      }
  };

  template <typename FloatType=double>
  class cos_sin_table
  {
    public:
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      cos_sin_table() {}

      cos_sin_table(int n_points, bool interpolate=false)
      :
        n_points_(n_points),
        n_points_4_(n_points/4),
        values_memory_(new float_type[n_points_+n_points_4_]),
        values_(values_memory_.get()),
        interpolate_(interpolate)
      {
        if (interpolate_) {
          diffvalues_memory_.reset(new float_type[n_points_ + n_points_4_ + 1]);
          diffvalues_ = diffvalues_memory_.get();
        }
        CCTBX_ASSERT(n_points % 4 == 0);
        using scitbx::constants::two_pi;
        FloatType two_pi_over_n_points = two_pi/n_points_;
        for(int i=0;i<n_points_+n_points_4_;i++) {
          FloatType value = std::sin(i*two_pi_over_n_points);
          values_[i] = value;
          if(interpolate_) {
            diffvalues_[i] = std::sin((i+1)*two_pi_over_n_points) - value;
          }
        }
      }

      int
      n_points() const { return n_points_; }

      complex_type
      get(float_type unary_angle) const
      {
        if(interpolate_) {
          float_type dec = unary_angle*n_points_;
          if (dec<0.0) {
            dec = -dec;
            int i = static_cast<int>(dec);
            float_type frac = dec - i;
            i = i % n_points_;
            return complex_type(
                values_[i+n_points_4_] + frac*diffvalues_[i+n_points_4_],
                - values_[i] - frac*diffvalues_[i]);
          }
          else {
            int i = static_cast<int>(dec);
            float_type frac = dec - i;
            i = i % n_points_;
            return complex_type(
                values_[i+n_points_4_] + frac*diffvalues_[i+n_points_4_],
                values_[i] + frac*diffvalues_[i]);
          }
        }
        else {
          int i = static_cast<int>(unary_angle*n_points_) % n_points_;
          if (i < 0) i += n_points_;
          return complex_type(values_[i+n_points_4_], values_[i]);
        }
      }

      complex_type operator()(float_type unary_angle) const {
        return get(unary_angle);
      }

    private:
      int n_points_;
      int n_points_4_;
      boost::shared_array<float_type> values_memory_;
      float_type* values_;
      bool interpolate_;
      boost::shared_array<float_type> diffvalues_memory_;
      float_type* diffvalues_;
  };

  template <typename FloatType=double>
  class cos_sin_lin_interp_table
  {
    public:
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      cos_sin_lin_interp_table() {}

      cos_sin_lin_interp_table(int n_points)
      :
        n_points_(n_points),
        n_points_4_(n_points/4),
        values_memory_(new float_type[n_points_+n_points_4_+1]),
        values_(values_memory_.get()),
        diffvalues_memory_(new float_type[n_points_+n_points_4_+1]),
        diffvalues_(diffvalues_memory_.get())
      {
        CCTBX_ASSERT(n_points % 4 == 0);
        using scitbx::constants::two_pi;
        for(int i=0;i<=n_points_+n_points_4_;i++)
        {
          values_[i] = std::sin(i*two_pi/n_points_);
          diffvalues_[i] = std::sin((i+1)*two_pi/n_points_)
                         - std::sin(i*two_pi/n_points_);
        }
      }

      int
      n_points() const { return n_points_; }

      complex_type
      get(float_type unary_angle) const
      {
        float_type dec = unary_angle*n_points_;
        if (dec<0.0)
        {
          dec = -dec;
          int i = static_cast<int>(dec); // floor it to an integer
          // get fraction between the ceiling and the floor of fval
          float_type frac = dec - i;
          // return interpolated value between the i-th and the (i+1)th entry
          i = i % n_points_;

          return complex_type(
              values_[i+n_points_4_] + frac*diffvalues_[i+n_points_4_],
              - values_[i] -frac*diffvalues_[i]);
        }
        else
        {
          int i = static_cast<int>(dec); // floor it to an integer
          float_type frac = dec - i;
          // return interpolated value between the i-th and the (i+1)th entry
          i = i % n_points_;

          return complex_type(
              values_[i+n_points_4_] + frac*diffvalues_[i+n_points_4_],
              values_[i] + frac*diffvalues_[i]);
        }
      }

      complex_type operator()(float_type unary_angle) const {
        return get(unary_angle);
      }

    private:
      int n_points_;
      int n_points_4_;
      boost::shared_array<float_type> values_memory_;
      float_type* values_;
      boost::shared_array<float_type> diffvalues_memory_;
      float_type* diffvalues_;
  };


#if defined(__sgi) && !defined(__GNUC__)
  namespace work_around_edg_typename_handling {
    typedef cos_sin_exact<double> cos_sin_exact_double;
    typedef cos_sin_table<double> cos_sin_table_double;
    typedef cos_sin_lin_interp_table<double> cos_sin_lin_interp_table_double;
  }
#endif

}} // namespace cctbx::math

#endif // CCTBX_MATH_COS_SIN_TABLE_H
