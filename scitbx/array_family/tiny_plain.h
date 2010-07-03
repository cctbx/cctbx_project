#ifndef SCITBX_ARRAY_FAMILY_TINY_PLAIN_H
#define SCITBX_ARRAY_FAMILY_TINY_PLAIN_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/detail/misc.h>
#include <scitbx/array_family/detail/tiny_helpers.h>
#include <scitbx/array_family/array_adaptor.h>

namespace scitbx { namespace af {

  // Automatic allocation, fixed size.
  template <typename ElementType, std::size_t N>
  class tiny_plain
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      BOOST_STATIC_CONSTANT(std::size_t, fixed_size=N);

      ElementType elems[N];

      tiny_plain() {}

      template <typename OtherArrayType>
      tiny_plain(array_adaptor<OtherArrayType> const& a_a)
      {
        OtherArrayType const& a = *(a_a.pointee);
        if (a.size() != N) throw_range_error();
        for(std::size_t i=0;i<N;i++) elems[i] = a[i];
      }

      template <typename OtherArrayType>
      tiny_plain(array_adaptor_with_static_cast<OtherArrayType> const& a_a)
      {
        OtherArrayType const& a = *(a_a.pointee);
        if (a.size() != N) throw_range_error();
        for(std::size_t i=0;i<N;i++) elems[i] = static_cast<ElementType>(a[i]);
      }

      SCITBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny_plain)
      SCITBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny_plain)

      static size_type size() { return N; }
      static bool empty() { return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(tiny_plain, elems, N)

      SCITBX_ARRAY_FAMILY_TAKE_REF(elems, N)

      void swap(ElementType* other) {
        std::swap(*this, other);
      }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_TINY_PLAIN_H
