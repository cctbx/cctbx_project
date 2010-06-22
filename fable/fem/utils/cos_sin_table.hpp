#ifndef FEM_UTILS_COS_SIN_TABLE_HPP
#define FEM_UTILS_COS_SIN_TABLE_HPP

#include <boost_adaptbx/error_utils.h>
#include <boost/noncopyable.hpp>
#include <fem/size_t.hpp>

namespace fem { namespace utils {

  template <typename FloatType>
  struct cos_sin_table : boost::noncopyable
  {
    typedef FloatType ft;

    protected:
      size_t n_points_;
      size_t n_points_2_;
      ft* values_;

      ft rad_as_i;
      ft sin_05_neg;
      ft sin_05_pos;

      public:

    ~cos_sin_table() { delete[] values_; }

    cos_sin_table() : values_(0) {}

    cos_sin_table(
      size_t n_points)
    :
      n_points_(n_points),
      n_points_2_(n_points / 2),
      values_(new ft[n_points_2_ + 1])
    {
      ASSERTBX(n_points % 4 == 0);
      ft two_pi = 8 * std::atan(ft(1));
      for(size_t i=0;i<=n_points_2_;i++) {
        values_[i] = std::cos(i * two_pi / n_points_);
      }
      size_t n_points_4 = n_points / 4;
      values_[n_points_4] = 0;
      rad_as_i = n_points_ / two_pi;
      sin_05_neg = 0.5f + n_points_4;
      sin_05_pos = 0.5f - n_points_4 + n_points_;
    }

    size_t
    n_points() const { return n_points_; }

    ft
    max_delta() const
    {
      ft pi = 4 * std::atan(ft(1));
      return std::sin(pi / n_points_);
    }

    ft const&
    cos(
      ft const& angle_rad) const
    {
      ssize_t i;
      if (angle_rad < 0) {
        i = static_cast<ssize_t>(0.5f - angle_rad * rad_as_i);
      }
      else {
        i = static_cast<ssize_t>(0.5f + angle_rad * rad_as_i);
      }
      i %= n_points_;
      if (i > n_points_2_) i = n_points_ - i;
      return values_[i];
    }

    ft const&
    sin(
      ft const& angle_rad) const
    {
      ssize_t i;
      if (angle_rad < 0) {
        i = static_cast<ssize_t>(sin_05_neg - angle_rad * rad_as_i);
      }
      else {
        i = static_cast<ssize_t>(sin_05_pos + angle_rad * rad_as_i);
      }
      i %= n_points_;
      if (i > n_points_2_) i = n_points_ - i;
      return values_[i];
    }
  };

}} // namespace fem::utils

#endif // GUARD
