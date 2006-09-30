#ifndef SCITBX_ARRAY_FAMILY_SMALL_H
#define SCITBX_ARRAY_FAMILY_SMALL_H

#include <scitbx/array_family/small_plain.h>
#include <scitbx/array_family/ref_reductions.h>

namespace scitbx { namespace af {

  // Automatic allocation, fixed size, standard operators.
  template <typename ElementType, std::size_t N>
  class small : public small_plain<ElementType, N>
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef small_plain<ElementType, N> base_class;

      small()
      {}

      explicit
      small(size_type const& sz)
        : base_class(sz)
      {}

      // non-std
      small(af::reserve const& sz)
        : base_class(sz)
      {}

      small(size_type const& sz, ElementType const& x)
        : base_class(sz, x)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      // non-std
      template <typename FunctorType>
      small(size_type const& sz, init_functor<FunctorType> const& ftor)
        : base_class(sz, ftor)
      {}
#endif

      small(const ElementType* first, const ElementType* last)
        : base_class(first, last)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      small(const OtherElementType* first, const OtherElementType* last)
        : base_class(first, last)
      {}
#endif

      template <typename OtherArrayType>
      small(array_adaptor<OtherArrayType> const& a_a)
        : base_class(a_a)
      {}

#     include <scitbx/array_family/detail/reducing_boolean_mem_fun.h>
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SMALL_H
