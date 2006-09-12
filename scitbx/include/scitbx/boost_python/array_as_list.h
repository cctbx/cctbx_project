#ifndef SCITBX_BOOST_PYTHON_ARRAY_AS_LIST_H
#define SCITBX_BOOST_PYTHON_ARRAY_AS_LIST_H

#include <boost/python/object.hpp>

namespace scitbx { namespace boost_python {

  template <typename ElementType>
  boost::python::object
  array_as_list(
    const ElementType* array,
    std::size_t array_size)
  {
    boost::python::object result(boost::python::handle<>(
      PyList_New(array_size)));
    PyObject* r = result.ptr();
    for(std::size_t i=0;i<array_size;i++) {
      boost::python::object ri(array[i]);
      PyList_SET_ITEM(r, i, boost::python::incref(ri.ptr()));
    }
    return result;
  }

}} // scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_ARRAY_AS_LIST_H
