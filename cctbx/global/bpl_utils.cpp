// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

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
