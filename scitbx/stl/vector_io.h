#ifndef SCITBX_STD_VECTOR_IO_H
#define SCITBX_STD_VECTOR_IO_H

#include <ostream>
#include <vector>

namespace scitbx {

template <typename ElementType>
std::ostream& operator << (std::ostream& o, std::vector<ElementType> const& a) {
  return o << af::const_ref<ElementType>(&a[0], a.size());
}

}

#endif
