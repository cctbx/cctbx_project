#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_BYTE_STR_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_BYTE_STR_H

#include <boost/python/str.hpp>
#include <boost/python/handle.hpp>
#include <scitbx/array_family/shared.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace af { namespace boost_python {

  template <typename ArrayType>
#ifdef IS_PY3K
  PyObject*
  copy_to_byte_str(
    ArrayType const& a)
  {
    return PyBytes_FromStringAndSize(
      reinterpret_cast<const char*>(a.begin()),
      (a.end()-a.begin()) * a.element_size()
      );
#else
  boost::python::str
  copy_to_byte_str(
    ArrayType const& a)
  {
    return boost::python::str(
      reinterpret_cast<const char*>(a.begin()),
      reinterpret_cast<const char*>(a.end()));
#endif
  }

  template <typename ArrayType>
#ifdef IS_PY3K
  PyObject*
  slice_to_byte_str(
    ArrayType const& a,
    std::size_t const& offset_begin,
    std::size_t const& offset_end)
  {
    SCITBX_ASSERT(offset_end <= a.size());
    SCITBX_ASSERT(offset_begin <= offset_end);
    return PyBytes_FromStringAndSize(
      reinterpret_cast<const char*>(a.begin() + offset_begin),
      (offset_end-offset_begin) * a.element_size()
      );
#else
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
#endif
  }

  template <typename ElementType>
  shared<ElementType>
  shared_from_byte_str(
    boost::python::object const& byte_str)
  {
#ifdef IS_PY3K
    PyObject* o(byte_str.ptr());
    // A unicode argument is re-encoded to UTF-8 bytes; keep that temporary
    // alive (and released) with handle<>. len_byte_str must be the size of
    // the byte buffer str_ptr points at, not len(byte_str): for non-ASCII
    // unicode the latter is the code-point count and would truncate the array.
    boost::python::handle<> utf8;
    if (PyUnicode_Check(o)) {
      utf8 = boost::python::handle<>(PyUnicode_AsUTF8String(o));
      o = utf8.get();
    }
    const char* str_ptr = PyBytes_AsString(o);
    boost::python::ssize_t len_byte_str = PyBytes_GET_SIZE(o);
#else
    const char* str_ptr = PyString_AsString(byte_str.ptr());
    boost::python::ssize_t len_byte_str = boost::python::len(byte_str);
#endif
    boost::python::ssize_t
      shared_array_size = len_byte_str / sizeof(ElementType);
    SCITBX_ASSERT(shared_array_size * sizeof(ElementType) == len_byte_str);
    return shared<ElementType>(
      reinterpret_cast<const ElementType*>(str_ptr),
      reinterpret_cast<const ElementType*>(str_ptr) + shared_array_size);
  }

}}} // namespace scitbx::af::boost_python

#endif // GUARD
