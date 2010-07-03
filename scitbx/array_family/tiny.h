#ifndef SCITBX_ARRAY_FAMILY_TINY_H
#define SCITBX_ARRAY_FAMILY_TINY_H

#include <scitbx/array_family/tiny_plain.h>
#include <scitbx/array_family/ref_reductions.h>

namespace scitbx { namespace af {

  // Automatic allocation, fixed size, standard operators.
  template <typename ElementType, std::size_t N>
  class tiny : public tiny_plain<ElementType, N>
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef tiny_plain<ElementType, N> base_class;

      tiny() {}

      template <typename OtherArrayType>
      tiny(array_adaptor<OtherArrayType> const& a_a)
      : base_class(a_a)
      {}

      template <typename OtherArrayType>
      tiny(array_adaptor_with_static_cast<OtherArrayType> const& a_a)
      : base_class(a_a)
      {}

      SCITBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny)
      SCITBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny)

#     include <scitbx/array_family/detail/reducing_boolean_mem_fun.h>
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_TINY_H
