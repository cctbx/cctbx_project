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
    long signed_sz = static_cast<long>(sz);
    if (sl.step.is_valid) {
      step = sl.step.value;
    }
    if (!sl.start.is_valid) {
      start = step < 0 ? signed_sz-1 : 0;
    }
    else {
      start = sl.start.value;
      if (start < 0) start += signed_sz;
    }
    if (!sl.stop.is_valid) {
      stop = step < 0 ? -1 : signed_sz;
    }
    else {
      stop = sl.stop.value;
      if (stop < 0) stop += signed_sz;
    }
    if (start > signed_sz-1) start = signed_sz;
    if (start < 0) start = 0;
    if      (stop < -1) stop = -1;
    else if (stop > signed_sz) stop = signed_sz;
    SCITBX_ASSERT(step != 0 || stop == start);
    long signed_size = stop - start + step;
    if      (step < 0) signed_size++;
    else if (step > 0) signed_size--;
    else {
      signed_size = 0;
      step = 1;
    }
    signed_size /= step;
    if (signed_size < 0) signed_size = 0;
    size = static_cast<std::size_t>(signed_size);
    stop = start + step * size;
  }

}} // namespace scitbx::boost_python
