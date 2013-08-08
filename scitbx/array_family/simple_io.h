#ifndef SCITBX_ARRAY_FAMILY_SIMPLE_IO_H
#define SCITBX_ARRAY_FAMILY_SIMPLE_IO_H

#include <ostream>
#include <scitbx/array_family/ref.h>

namespace scitbx { namespace af {

  template <typename ElementType, typename AccessorType>
  std::ostream&
  operator<<(std::ostream& os, const_ref<ElementType, AccessorType> const& a) {
    os << "{";
    if (a.size() > 0) {
      for (std::size_t i = 0;;) {
        os << a[i];
        i++;
        if (i == a.size()) break;
        os << ",";
      }
    }
    os << "}";
    return os;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SIMPLE_IO_H
