// $Id$

#include <boost/python/class_builder.hpp>

namespace bpl_utils {

  boost::python::tuple tuple_from_python_list_or_tuple(PyObject* p) {
    using namespace boost::python;
    if (PyList_Check(p))
      return tuple(ref(PyList_AsTuple(p)));
    else
      return tuple(ref(p, ref::increment_count));
  }
}
