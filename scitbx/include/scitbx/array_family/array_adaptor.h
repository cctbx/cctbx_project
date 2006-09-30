#ifndef SCITBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H
#define SCITBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H

namespace scitbx { namespace af {

  template <typename T>
  struct array_adaptor
  {
    const T* pointee;
    array_adaptor(T const& a) : pointee(&a) {}
  };

  template <typename T>
  inline
  array_adaptor<T>
  adapt(T const& a)
  {
    return array_adaptor<T>(a);
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H
