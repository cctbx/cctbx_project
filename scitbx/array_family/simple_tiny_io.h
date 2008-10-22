#ifndef SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
#define SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H

#include <scitbx/array_family/simple_io.h>

namespace scitbx { namespace af {

  template <template<typename, std::size_t> class TinyOrSmall,
            typename  ElementType, std::size_t N>
  std::ostream&
  operator<<(std::ostream& os, TinyOrSmall<ElementType, N> const& a) {
    return os << a.const_ref();
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
