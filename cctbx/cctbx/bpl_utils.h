// $Id$

#ifndef CCTBX_BPL_UTILS_H
#define CCTBX_BPL_UTILS_H

#include <boost/python/class_builder.hpp>

namespace bpl_utils {

  boost::python::tuple tuple_from_python_list_or_tuple(PyObject* p);
}

#endif // CCTBX_BPL_UTILS_H
