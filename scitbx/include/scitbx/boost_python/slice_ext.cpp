#include <boost/python/module.hpp>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include <scitbx/boost_python/slice.h>

namespace scitbx { namespace boost_python {

namespace {

  slice_item
  extract_slice_item(boost::python::object const& obj)
  {
    slice_item result;
    if (obj.ptr() == Py_None) {
      result.is_valid = false;
      result.value = 0;
    }
    else {
      result.is_valid = true;
      boost::python::extract<long> proxy(obj);
      result.value = proxy();
    }
    return result;
  }

  struct slice_from_python
  {
    slice_from_python()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<slice>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      if (!PySlice_Check(obj_ptr)) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<slice>*)
          data)->storage.bytes;
      new (storage) slice();
      data->convertible = storage;
      boost::python::handle<> slice_hdl(boost::python::borrowed(obj_ptr));
      boost::python::object slice_obj = boost::python::object(slice_hdl);
      slice* self = (slice*) storage;
      self->start = extract_slice_item(slice_obj.attr("start"));
      self->stop = extract_slice_item(slice_obj.attr("stop"));
      self->step = extract_slice_item(slice_obj.attr("step"));
    }
  };

  void init_module()
  {
    using namespace boost::python;
    slice_from_python();
  }

}}} // namespace scitbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_boost_python_slice_ext)
{
  scitbx::boost_python::init_module();
}
