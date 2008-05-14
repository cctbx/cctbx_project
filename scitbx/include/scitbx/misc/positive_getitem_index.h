#ifndef SCITBX_MISC_POSITIVE_GETITEM_INDEX_H
#define SCITBX_MISC_POSITIVE_GETITEM_INDEX_H

#include <stdexcept>

namespace scitbx {

  template <typename IndexType, typename SizeType>
  SizeType
  positive_getitem_index(
    IndexType const& i,
    SizeType const& size,
    bool allow_i_eq_size=false,
    const char* index_out_of_range = "Index out of range.")
  {
    if (i >= 0) {
      SizeType result = static_cast<SizeType>(i);
      if (   result > size
          || (result == size && !allow_i_eq_size)) {
        throw std::out_of_range(index_out_of_range);
      }
      return result;
    }
    if (static_cast<SizeType>(-i) > size) {
      throw std::out_of_range(index_out_of_range);
    }
    return size + i;
  }

} // namespace scitbx

#endif // SCITBX_MISC_POSITIVE_GETITEM_INDEX_H
