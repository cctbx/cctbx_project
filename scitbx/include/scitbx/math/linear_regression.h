/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy of cctbx/math/linear_regression.h (rwgk)
     2002 Mar: modified version of cctbx/vector/linear_regresion.h (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_MATH_LINEAR_REGRESSION_H
#define SCITBX_MATH_LINEAR_REGRESSION_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/error.h>
#include <cmath>

namespace scitbx { namespace math {

  template <typename FloatType = double,
            typename IntegerType = std::size_t>
  class linear_regression_core
  {
    public:
      typedef FloatType float_type;
      typedef std::size_t size_type;

      void
      reset()
      {
        is_well_defined_ = false;
        y_intercept_ = slope_ = FloatType(0);
      }

      void
      set(IntegerType const& n,
          FloatType const& min_x, FloatType const& max_x,
          FloatType const& min_y, FloatType const& max_y,
          FloatType const& sum_x, FloatType const& sum_x2,
          FloatType const& sum_y, FloatType const& sum_y2,
          FloatType const& sum_xy,
          FloatType const& epsilon = 1.e-15);

      //! Flag.
      bool
      is_well_defined() const { return is_well_defined_; }

      //! Linear regression coefficients b: y = mx + b.
      FloatType
      y_intercept() const { return y_intercept_; }

      //! Linear regression coefficients m: y = mx + b.
      FloatType
      slope() const { return slope_; }

    protected:
      bool is_well_defined_;
      FloatType y_intercept_;
      FloatType slope_;
  };

  template <typename FloatType,
            typename IntegerType>
  void
  linear_regression_core<FloatType, IntegerType>
  ::set(IntegerType const& n,
        FloatType const& min_x, FloatType const& max_x,
        FloatType const& min_y, FloatType const& max_y,
        FloatType const& sum_x, FloatType const& sum_x2,
        FloatType const& sum_y, FloatType const& sum_y2,
        FloatType const& sum_xy,
        FloatType const& epsilon)
  {
    reset();
    if (n < 1) return;
    if (min_x == max_x) return;
    if (min_y == max_y) {
      y_intercept_ = min_y;
      is_well_defined_ = true;
      return;
    }
    FloatType fn(n);
    FloatType dx = std::max(fn::absolute(min_x - sum_x / fn),
                            fn::absolute(max_x - sum_x / fn));
    FloatType dy = std::max(fn::absolute(min_y - sum_y / fn),
                            fn::absolute(max_y - sum_y / fn));
    if (dx == 0) return;
    if (dy == 0) {
      y_intercept_ = sum_y / fn;
      is_well_defined_ = true;
      return;
    }
    if (dx < dy * epsilon) return;
    FloatType d = fn * sum_x2 - sum_x * sum_x;
    if (d != 0) {
      y_intercept_ = (sum_x2 * sum_y - sum_x * sum_xy) / d;
      slope_ = (fn * sum_xy - sum_x * sum_y) / d;
    }
    is_well_defined_ = true;
  }

  template <typename FloatType = double>
  class linear_regression : public linear_regression_core<FloatType>
  {
    public:
      typedef FloatType float_type;
      typedef typename linear_regression_core<FloatType>::size_type size_type;

      linear_regression() {}

      linear_regression(af::const_ref<FloatType> const& x,
                        af::const_ref<FloatType> const& y,
                        FloatType const& epsilon = 1.e-15);
  };

  template <typename FloatType>
  linear_regression<FloatType>
  ::linear_regression(af::const_ref<FloatType> const& x,
                      af::const_ref<FloatType> const& y,
                      FloatType const& epsilon)
  {
    size_type n = x.size();
    SCITBX_ASSERT(n == y.size());
    if (n == 0) {
      this->reset();
      return;
    }
    float_type min_x = x[0];
    float_type max_x = x[0];
    float_type min_y = y[0];
    float_type max_y = y[0];
    float_type sum_x = x[0];
    float_type sum_x2 = x[0] * x[0];
    float_type sum_y = y[0];
    float_type sum_y2 = y[0] * y[0];
    float_type sum_xy = x[0] * y[0];
    for(size_type i=1;i<n;i++) {
      if (min_x > x[i]) min_x = x[i];
      if (max_x < x[i]) max_x = x[i];
      if (min_y > y[i]) min_y = y[i];
      if (max_y < y[i]) max_y = y[i];
      sum_x += x[i];
      sum_x2 += x[i] * x[i];
      sum_y += y[i];
      sum_y2 += y[i] * y[i];
      sum_xy += x[i] * y[i];
    }
    this->set(n, min_x, max_x, min_y, max_y,
              sum_x, sum_x2, sum_y, sum_y2, sum_xy,
              epsilon);
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_LINEAR_REGRESSION_H
