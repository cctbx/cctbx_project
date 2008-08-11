#include <scitbx/math/linear_interpolation.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>

namespace scitbx { namespace af {

  template <typename FloatType>
  FloatType
  linear_interpolation(
    const_ref<FloatType> const& table_x,
    const_ref<FloatType> const& table_y,
    FloatType const& x,
    FloatType const& tolerance=1.e-6)
  {
    SCITBX_ASSERT(table_x.size() == table_y.size());
    SCITBX_ASSERT(table_x.size() > 0);
    SCITBX_ASSERT(tolerance >= 0);
    FloatType range_x = table_x.back() - table_x.front();
    if (x <= table_x.front()) {
      if (table_x.front() - x < range_x * tolerance) {
        if (table_x.size() == 1) return table_y[0];
        return math::linear_interpolation(
          x, table_x[0], table_x[1], table_y[0], table_y[1]);
      }
      throw error("x-value smaller than smallest table_x.");
    }
    if (x >= table_x.back()) {
      if (x - table_x.back() < range_x * tolerance) {
        if (table_x.size() == 1) return table_y[0];
        std::size_t j = table_x.size()-1;
        return math::linear_interpolation(
          x, table_x[j-1], table_x[j], table_y[j-1], table_y[j]);
      }
      throw error("x-value larger than largest table_x.");
    }
    for(std::size_t j=1;j<table_x.size();j++) {
      if (table_x[j] > x) {
        return math::linear_interpolation(
          x, table_x[j-1], table_x[j], table_y[j-1], table_y[j]);
      }
      else {
        if (table_x[j-1] >= table_x[j]) {
          throw error("table_x not strictly increasing.");
        }
      }
    }
    throw SCITBX_INTERNAL_ERROR();
  }

  template <typename FloatType>
  shared<FloatType>
  linear_interpolation(
    const_ref<FloatType> const& table_x,
    const_ref<FloatType> const& table_y,
    const_ref<FloatType> const& x,
    FloatType const& tolerance=1.e-6)
  {
    shared<FloatType> result((reserve(x.size())));
    for(std::size_t i=0;i<x.size();i++) {
      result.push_back(
        linear_interpolation(table_x, table_y, x[i], tolerance));
    }
    return result;
  }

}} // namespace scitbx::af
