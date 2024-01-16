#ifndef GLTBX_POINTER_ARGS_BPL_H
#define GLTBX_POINTER_ARGS_BPL_H

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/str.hpp>
#include <boost/shared_array.hpp>
#include <stdexcept>
#include <sstream>
#include <vector>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace gltbx { namespace boost_python {

  typedef boost::python::ssize_t py_ssize_t;

  namespace detail {

    inline
    unsigned
    consolidate_sizes(
      const char* arg_name,
      py_ssize_t expected_size,
      py_ssize_t given_size)
    {
      if (expected_size != 0
          && given_size != 0
          && given_size != expected_size) {
        std::ostringstream o;
        o << "Argument \"" << arg_name
            << "\" has the wrong number of elements:\n"
          << "  expected size: " << expected_size << "\n"
          << "     given size: " << given_size;
        throw std::runtime_error(o.str());
      }
      if (expected_size == 0) return given_size;
      return expected_size;
    }

  } // namespace detail

  template <typename T>
  struct converter
  {
    const char* arg_name_;
    boost::python::object py_arg_;
    bool is_const_;
    PyObject* py_arg_ptr_;
    py_ssize_t len_py_arg_;
    std::vector<T> data_;

    void
    process_size(py_ssize_t expected_size, py_ssize_t len_py_arg)
    {
      len_py_arg_ = len_py_arg;
      py_ssize_t data_size = detail::consolidate_sizes(
        arg_name_, expected_size, len_py_arg);
      if (len_py_arg_ == 0) data_.resize(data_size, 0);
      else                  data_.reserve(data_size);
    }

    void
    extract_element(PyObject* py_elem_ptr)
    {
      boost::python::handle<>
        py_elem_hdl(boost::python::borrowed(py_elem_ptr));
      boost::python::object py_elem_obj(py_elem_hdl);
      boost::python::extract<T> elem_proxy(py_elem_obj);
      if (!elem_proxy.check()) {
        std::ostringstream o;
        o << "Argument \""
          << arg_name_ << "\" has one or more elements of the wrong type.";
        throw std::runtime_error(o.str());
      }
      data_.push_back(elem_proxy());
    }

    converter(
      const char* arg_name,
      boost::python::object const& py_arg,
      py_ssize_t expected_size,
      bool is_const)
    :
      arg_name_(arg_name),
      py_arg_(py_arg),
      is_const_(is_const),
      py_arg_ptr_(py_arg.ptr())
    {
      if (!is_const_ && !PyList_Check(py_arg_ptr_)) {
        throw std::runtime_error(
          std::string(arg_name_)+" must be a Python list.");
      }
      if (PyList_Check(py_arg_ptr_)) {
        process_size(expected_size, PyList_GET_SIZE(py_arg_ptr_));
        for(py_ssize_t i=0;i<len_py_arg_;i++) {
          extract_element(PyList_GET_ITEM(py_arg_ptr_, i));
        }
      }
      else if (PyTuple_Check(py_arg_ptr_)) {
        process_size(expected_size, PyTuple_GET_SIZE(py_arg_ptr_));
        for(py_ssize_t i=0;i<len_py_arg_;i++) {
          extract_element(PyTuple_GET_ITEM(py_arg_ptr_, i));
        }
      }
      else {
        throw std::runtime_error(
          std::string(arg_name_)+"must be a Python list or tuple.");
      }
    }

    T*
    get() { return (data_.size() ? &*data_.begin() : 0); }

    void
    write_back()
    {
      namespace bp = boost::python;
      py_ssize_t sz = data_.size();
      for(py_ssize_t i=0;i<sz;i++) {
        bp::object item(data_[i]);
        if (len_py_arg_ == 0) {
          if (PyList_Append(py_arg_ptr_, item.ptr()) != 0) {
            bp::throw_error_already_set();
          }
        }
        else {
          if (PyList_SetItem(py_arg_ptr_, i, bp::incref(item.ptr())) != 0) {
            bp::throw_error_already_set();
          }
        }
      }
    }
  };

  template <typename T>
  struct converter_str
  {
    const char* arg_name_;
    boost::python::object py_arg_;
    bool is_const_;
    PyObject* py_arg_ptr_;
    py_ssize_t len_py_arg_;
    py_ssize_t data_size_;
    boost::shared_array<T> data_;

    void
    throw_must_be() const
    {
      if (!is_const_) {
        throw std::runtime_error(
          std::string(arg_name_)
          + " must be a Python list with one string element.");
      }
      throw std::runtime_error(
        std::string(arg_name_)
        + " must be a Python string or list with one string element.");
    }

    converter_str(
      const char* arg_name,
      boost::python::object const& py_arg,
      py_ssize_t expected_size,
      bool is_const)
    :
      arg_name_(arg_name),
      py_arg_(py_arg),
      is_const_(is_const),
      py_arg_ptr_(py_arg.ptr()),
      len_py_arg_(0)
    {
      PyObject* item = 0;
#ifdef IS_PY3K
      bool is_unicode = false;
#endif
      if (PyList_Check(py_arg_ptr_)) {
        len_py_arg_ = PyList_GET_SIZE(py_arg_ptr_);
        if (len_py_arg_ == 1) {
          item = PyList_GET_ITEM(py_arg_ptr_, 0);
#ifdef IS_PY3K
          if (!(PyBytes_Check(item) || PyUnicode_Check(item))) throw_must_be();
          if (PyUnicode_Check(item)) is_unicode = false;
#else
          if (!PyString_Check(item)) throw_must_be();
#endif
        }
        else if (len_py_arg_ != 0 || is_const_) {
          throw_must_be();
        }
      }
#ifdef IS_PY3K
      else if (is_const_ && (PyBytes_Check(py_arg_ptr_) || PyUnicode_Check(py_arg_ptr_))) {
        item = py_arg_ptr_;
        if (PyUnicode_Check(item)) is_unicode = false;
#else
      else if (is_const_ && PyString_Check(py_arg_ptr_)) {
        item = py_arg_ptr_;
#endif
      }
      else {
        throw_must_be();
      }
      py_ssize_t item_size = 0;
      if (item != 0) {
#ifdef IS_PY3K
        if (is_unicode)
          item_size = PyUnicode_GET_LENGTH(item);
        else
          item_size = PyBytes_GET_SIZE(item);
#else
        item_size = PyString_GET_SIZE(item);
#endif
        data_size_ = detail::consolidate_sizes(
          arg_name_, expected_size, item_size);
      }
      else if (expected_size != 0) {
        data_size_ = expected_size;
      }
      data_ = boost::shared_array<T>(new T[data_size_]);
      if (item != 0) {
#ifdef IS_PY3K
        char* item_char_pointer;
        if (is_unicode)
          // If encoding fails, no error will be reported. Revisit? XXX
          item_char_pointer = PyBytes_AsString(PyUnicode_EncodeLocale(item, NULL));
        else
          item_char_pointer = PyBytes_AsString(item);
#else
        char* item_char_pointer = PyString_AsString(item);
#endif
        if (item_char_pointer == 0) {
          boost::python::throw_error_already_set();
        }
        for(py_ssize_t i=0;i<item_size;i++) {
          data_[i] = static_cast<T>(item_char_pointer[i]);
        }
      }
      for(py_ssize_t i=item_size;i<data_size_;i++) {
        data_[i] = static_cast<T>(0);
      }
    }

    T*
    get() { return data_.get(); }

    void
    write_back()
    {
      namespace bp = boost::python;
#ifdef IS_PY3K
      bp::object new_item((bp::handle<>(PyBytes_FromStringAndSize(
#else
      bp::object new_item((bp::handle<>(PyString_FromStringAndSize(
#endif
        reinterpret_cast<char*>(data_.get()),
        data_size_ * sizeof(T)))));
      if (len_py_arg_ == 0) {
        if (PyList_Append(py_arg_ptr_, new_item.ptr()) != 0) {
          bp::throw_error_already_set();
        }
      }
      else {
        if (PyList_SetItem(py_arg_ptr_, 0, bp::incref(new_item.ptr())) != 0) {
          bp::throw_error_already_set();
        }
      }
    }
  };

}} // namespace gltbx::boost_python

#endif // GUARD
