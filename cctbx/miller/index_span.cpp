#include <cctbx/miller/index_span.h>
#include <scitbx/math/utils.h>
#include <scitbx/array_family/misc_functions.h>

namespace cctbx { namespace miller {

  index_span::index_span(af::const_ref<index<> > const& indices)
  {
    this->fill(min_end(0, 0));
    if (indices.size()) {
      for(std::size_t j=0;j<3;j++) {
        (*this)[j] = min_end().fill(indices[0][j]);
      }
    }
    for(std::size_t i=1;i<indices.size();i++) {
      for(std::size_t j=0;j<3;j++) {
        scitbx::math::update_min((*this)[j][0], indices[i][j]);
        scitbx::math::update_max((*this)[j][1], indices[i][j]);
      }
    }
    for(std::size_t j=0;j<3;j++) (*this)[j][1]++;
  }

  af::int3
  index_span::min() const
  {
    af::int3 result;
    for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][0];
    return result;
  }

  af::int3
  index_span::max() const
  {
    af::int3 result;
    for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][1] - 1;
    return result;
  }

  af::int3
  index_span::abs_range() const
  {
    af::int3 result;
    std::size_t j;
    for(j=0;j<3;j++) {
      result[j] = scitbx::fn::absolute((*this)[j][0]);
      scitbx::math::update_max(
        result[j], scitbx::fn::absolute((*this)[j][1]-1));
    }
    for(j=0;j<3;j++) result[j] += 1;
    return result;
  }

  af::int3
  index_span::map_grid() const
  {
    af::int3 result = abs_range();
    for(std::size_t j=0;j<3;j++) {
      result[j] = (result[j] - 1) * 2 + 1;
    }
    return result;
  }

  bool
  index_span::is_in_domain(index<> const& h) const
  {
    for(std::size_t j=0;j<3;j++) {
      if (!((*this)[j][0] <= h[j] && h[j] < (*this)[j][1])) return false;
    }
    return true;
  }

}} // namespace cctbx::miller
