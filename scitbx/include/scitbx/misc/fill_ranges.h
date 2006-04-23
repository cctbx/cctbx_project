#ifndef SCITBX_MISC_FILL_RANGES_H
#define SCITBX_MISC_FILL_RANGES_H

#include <algorithm>

namespace scitbx { namespace misc {

  template <
    typename BeginEndType,
    typename BoundaryConstIteratorType,
    typename FillElementType>
  void
  fill_ranges(
    BeginEndType begin,
    BeginEndType const& end,
    BoundaryConstIteratorType boundaries,
    BoundaryConstIteratorType const& boundaries_end,
    FillElementType* fill_begin)
  {
    FillElementType id=0;
    for(; boundaries != boundaries_end; id++) {
      BeginEndType i_transition = static_cast<BeginEndType>(*boundaries++);
      if (i_transition > begin) {
        if (i_transition < end) {
          FillElementType* fill_end = fill_begin + (i_transition-begin);
          std::fill(fill_begin, fill_end, id);
          begin = i_transition;
          fill_begin = fill_end;
        }
        else {
          break;
        }
      }
    }
    FillElementType* fill_end = fill_begin + (end-begin);
    std::fill(fill_begin, fill_end, id);
  }

}} // namespace scitbx::misc

#endif // SCITBX_MISC_FILL_RANGES_H
