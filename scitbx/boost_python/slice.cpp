/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/error.h>

namespace scitbx { namespace boost_python {

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
      boost::python::extract<int> proxy(obj);
      result.value = proxy();
    }
    return result;
  }

  slice_from_python::slice_from_python()
  {
    boost::python::converter::registry::push_back(
      &convertible,
      &construct,
      boost::python::type_id<slice>());
  }

  void*
  slice_from_python::convertible(PyObject* obj_ptr)
  {
    if (!PySlice_Check(obj_ptr)) return 0;
    return obj_ptr;
  }

  void
  slice_from_python::construct(
    PyObject* obj_ptr,
    boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    void* storage = (
      (boost::python::converter::rvalue_from_python_storage<slice>*)
        data)->storage.bytes;
    new (storage) slice();
    data->convertible = storage;
    slice& result = *((slice*)storage);
    boost::python::object slice_obj = boost::python::object(
      boost::python::handle<>(
        boost::python::borrowed(obj_ptr)));
    result.start = extract_slice_item(slice_obj.attr("start"));
    result.stop = extract_slice_item(slice_obj.attr("stop"));
    result.step = extract_slice_item(slice_obj.attr("step"));
  }

  adapted_slice::adapted_slice(slice const& sl, std::size_t sz)
  :
    step(1),
    size(0)
  {
    if (sl.step.is_valid) {
      step = sl.step.value;
    }
    SCITBX_ASSERT(step != 0);
    int sign;
    if (step > 0) {
      sign = 1;
      start = 0;
      stop = static_cast<long>(sz);
    }
    else {
      sign = -1;
      start = static_cast<long>(sz) - 1;
      stop = -1;
    }
    if (sl.start.is_valid) {
      start = positive_getitem_index(sl.start.value, sz);
    }
    if (sl.stop.is_valid) {
      stop = positive_getitem_index(sl.stop.value, sz);
    }
    if (sign * stop > sign * start) {
      size = (stop - start) / step;
      if (stop != start + step * size) size++;
    }
    stop = start + step * size;
  }

}} // namespace scitbx::boost_python
