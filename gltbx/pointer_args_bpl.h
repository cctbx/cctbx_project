#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/str.hpp>
#include <boost/shared_array.hpp>
#include <vector>
#include <stdexcept>
#include <cstdio>

namespace gltbx { namespace boost_python {

  namespace detail {

    inline
    unsigned
    consolidate_sizes(
      const char* arg_name,
      unsigned expected_size,
      unsigned given_size)
    {
      if (expected_size != 0
          && given_size != 0
          && given_size != expected_size) {
        char msg[512];
        std::sprintf(msg,
          "Argument \"%s\" has the wrong number of elements:\n"
          "  expected size: %d\n"
          "     given size: %d", arg_name, expected_size, given_size);
        throw std::runtime_error(msg);
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
    unsigned len_py_arg_;
    std::vector<T> data_;

    void
    process_size(unsigned expected_size, unsigned len_py_arg)
    {
      len_py_arg_ = len_py_arg;
      unsigned data_size = detail::consolidate_sizes(
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
        char msg[512];
        std::sprintf(msg,
          "Argument \"%s\" has one or more elements of the wrong type.",
          arg_name_);
        throw std::runtime_error(msg);
      }
      data_.push_back(elem_proxy());
    }

    converter(
      const char* arg_name,
      boost::python::object const& py_arg,
      unsigned expected_size,
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
        for(unsigned i=0;i<len_py_arg_;i++) {
          extract_element(PyList_GET_ITEM(py_arg_ptr_, i));
        }
      }
      else if (PyTuple_Check(py_arg_ptr_)) {
        process_size(expected_size, PyTuple_GET_SIZE(py_arg_ptr_));
        for(unsigned i=0;i<len_py_arg_;i++) {
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
      unsigned sz = data_.size();
      for(unsigned i=0;i<sz;i++) {
        boost::python::object item(data_[i]);
        PyObject* item_ptr = boost::python::incref(item.ptr());
        if (len_py_arg_ == 0) {
          if (PyList_Append(py_arg_ptr_, item_ptr) != 0) {
            boost::python::throw_error_already_set();
          }
        }
        else {
          if (PyList_SetItem(py_arg_ptr_, i, item_ptr) != 0) {
            boost::python::throw_error_already_set();
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
    unsigned data_size_;
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
      unsigned expected_size,
      bool is_const)
    :
      arg_name_(arg_name),
      py_arg_(py_arg),
      is_const_(is_const),
      py_arg_ptr_(py_arg.ptr())
    {
      if (!is_const_ && !PyList_Check(py_arg_ptr_)) throw_must_be();
      PyObject* item;
      if (PyList_Check(py_arg_ptr_)) {
        if (PyList_GET_SIZE(py_arg_ptr_) != 1) throw_must_be();
        item = PyList_GET_ITEM(py_arg_ptr_, 0);
        if (!PyString_Check(item)) throw_must_be();
      }
      else if (PyString_Check(py_arg_ptr_)) {
        item = py_arg_ptr_;
      }
      else {
        throw_must_be();
      }
      unsigned item_size = PyString_GET_SIZE(item);
      data_size_ = detail::consolidate_sizes(
        arg_name_, expected_size, item_size);
      data_ = boost::shared_array<T>(new T[data_size_]);
      char* item_char_pointer = PyString_AsString(item);
      if (item_char_pointer == 0) {
        boost::python::throw_error_already_set();
      }
      for(unsigned i=0;i<item_size;i++) {
        data_[i] = static_cast<T>(item_char_pointer[i]);
      }
      for(unsigned i=item_size;i<data_size_;i++) {
        data_[i] = static_cast<T>(0);
      }
    }

    T*
    get() { return data_.get(); }

    void
    write_back()
    {
      boost::shared_array<char> buf(new char[data_size_]);
      for(unsigned i=0;i<data_size_;i++) {
        buf[i] = static_cast<char>(data_[i]);
      }
      PyObject* new_item = PyString_FromStringAndSize(buf.get(), data_size_);
      if (new_item == 0) {
        boost::python::throw_error_already_set();
      }
      if (PyList_SetItem(py_arg_ptr_, 0, new_item) != 0) {
        boost::python::throw_error_already_set();
      }
    }
  };

}} // namespace gltbx::boost_python
