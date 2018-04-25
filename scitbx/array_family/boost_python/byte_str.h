#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_BYTE_STR_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_BYTE_STR_H

#include <boost/python/str.hpp>
#include <scitbx/array_family/shared.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace af { namespace boost_python {

  template <typename ArrayType>
  boost::python::str
  copy_to_byte_str(
    ArrayType const& a)
  {
    return boost::python::str(
      reinterpret_cast<const char*>(a.begin()),
      reinterpret_cast<const char*>(a.end()));
  }

  template <typename ArrayType>
  boost::python::str
  slice_to_byte_str(
    ArrayType const& a,
    std::size_t const& offset_begin,
    std::size_t const& offset_end)
  {
    SCITBX_ASSERT(offset_end <= a.size());
    SCITBX_ASSERT(offset_begin <= offset_end);
    return boost::python::str(
      reinterpret_cast<const char*>(a.begin() + offset_begin),
      reinterpret_cast<const char*>(a.begin() + offset_end));
  }

  template <typename ElementType>
  shared<ElementType>
  shared_from_byte_str(
    boost::python::str const& byte_str)
  {
    boost::python::ssize_t
      len_byte_str = boost::python::len(byte_str);
    boost::python::ssize_t
      shared_array_size = len_byte_str / sizeof(ElementType);
    SCITBX_ASSERT(shared_array_size * sizeof(ElementType) == len_byte_str);
#ifdef IS_PY3K
    const char* str_ptr = PyBytes_AsString(byte_str.ptr());
#else
    const char* str_ptr = PyString_AsString(byte_str.ptr());
#endif
    return shared<ElementType>(
      reinterpret_cast<const ElementType*>(str_ptr),
      reinterpret_cast<const ElementType*>(str_ptr) + shared_array_size);
  }

}}} // namespace scitbx::af::boost_python

#endif // GUARD
