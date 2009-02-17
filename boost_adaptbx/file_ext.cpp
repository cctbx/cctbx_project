#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/handle.hpp>
#include <boost_adaptbx/python_file_stream.h>

namespace boost_adaptbx { namespace file_conversion {

  // Boost.Python conversion dark magic
  struct python_file_to_stream_buffer
  {
    static void register_conversion() {
      using namespace boost::python;
      converter::registry::push_back(
        &convertible,
        &construct,
        type_id<python_file_buffer>());
    }

    static void *convertible(PyObject *obj_ptr) {
      using namespace boost::python;
      if (!(   PyObject_HasAttrString(obj_ptr, "read")
            && PyObject_HasAttrString(obj_ptr, "readline")
            && PyObject_HasAttrString(obj_ptr, "readlines"))
          &&
          !(   PyObject_HasAttrString(obj_ptr, "write")
            && PyObject_HasAttrString(obj_ptr, "writelines"))) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject *obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data *data)
    {
      using namespace boost::python;
      typedef converter::rvalue_from_python_storage<python_file_buffer> rvalue_t;
      void *storage = ((rvalue_t *) data)->storage.bytes;
      object python_file(handle<>(borrowed(obj_ptr)));
      new (storage) python_file_buffer(python_file);
      data->convertible = storage;
    }
  };


  struct python_file_buffer_wrapper
  {
    typedef python_file_buffer wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("buffer", no_init)
        .def_readwrite("size", wt::buffer_size,
                       "The size of the buffer sitting "
                       "between a Python file object and a C++ stream.")
      ;
    }
  };


}} // boost_adaptbx::file_conversions

BOOST_PYTHON_MODULE(python_file_ext)
{
  using namespace boost_adaptbx::file_conversion;
  python_file_to_stream_buffer::register_conversion();
  python_file_buffer_wrapper::wrap();
}
