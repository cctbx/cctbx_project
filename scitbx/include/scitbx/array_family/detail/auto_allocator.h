#ifndef SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
#define SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H

#include <boost/type_traits/alignment_traits.hpp>

namespace scitbx { namespace af { namespace detail {

  template <typename T, std::size_t N>
  union auto_allocator
  {
    typedef typename boost::type_with_alignment<
      (boost::alignment_of<T>::value)>::type
        align_t;
    align_t dummy_;
    char buffer[N * sizeof(T)];
  };

}}} // namespace scitbx::af::detail

#endif // SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
