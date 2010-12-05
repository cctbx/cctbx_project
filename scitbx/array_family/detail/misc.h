#ifndef SCITBX_ARRAY_FAMILY_MISC_H
#define SCITBX_ARRAY_FAMILY_MISC_H

#include <scitbx/array_family/type_traits.h>
#include <memory>

namespace scitbx { namespace af {

  class reserve
  {
    public:
      reserve(std::size_t size) : size_(size) {}

      std::size_t
      operator()() const { return size_; }

    private:
      std::size_t size_;
  };

  template <typename FunctorType>
  struct init_functor
  {
    init_functor() : held(0) {}
    explicit init_functor(FunctorType const& ftor) : held(&ftor) {}
    const FunctorType* held;
  };

  template <typename FunctorType>
  inline
  init_functor<FunctorType>
  make_init_functor(FunctorType const& ftor)
  {
    return init_functor<FunctorType>(ftor);
  }

  template <typename ElementType>
  struct init_functor_null : init_functor<init_functor_null<ElementType> >
  {
    init_functor_null() { this->held = this; }

    void operator()(ElementType*, std::size_t const&) const
    {
      SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
    }
  };

  namespace detail {

    template <class ElementType>
    inline
    void destroy_array_element(ElementType* elem, false_type) {
      elem->~ElementType();
    }

    template <class ElementType>
    inline
    void destroy_array_element(ElementType* /*elem*/, true_type) {
    }

    template <class ElementType>
    inline
    void destroy_array_elements(ElementType* first, ElementType* last,
                                false_type) {
      while (first != last) {
        first->~ElementType();
        ++first;
      }
    }

    template <class ElementType>
    inline
    void destroy_array_elements(ElementType*, ElementType*,
                                true_type) {
    }

  } // namespace detail

#if !defined(__GNUC__) \
    || ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 0)))
#define SCITBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
#endif

  template <typename InputElementType,
            typename OutputElementType>
  inline
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
#ifdef SCITBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
    return std::copy(first, last, result);
#else
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
#endif
  }

  template <typename InputElementType,
            typename OutputElementType>
#ifdef SCITBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
  inline
#endif
  OutputElementType*
  uninitialized_copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
#ifdef SCITBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
    return std::uninitialized_copy(first, last, result);
#else
    OutputElementType* p = result;
    try {
      for (; first != last; p++, first++) {
        new (p) OutputElementType(*first);
      }
    }
    catch (...) {
      detail::destroy_array_elements(result, p, false_type());
      throw;
    }
    return result;
#endif
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_MISC_H
