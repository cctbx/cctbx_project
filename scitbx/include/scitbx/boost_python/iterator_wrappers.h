/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_BOOST_PYTHON_ITERATOR_WRAPPERS_H
#define SCITBX_BOOST_PYTHON_ITERATOR_WRAPPERS_H

#include <boost/python/class.hpp>

namespace scitbx { namespace boost_python {

  inline
  boost::python::object
  pass_through(boost::python::object const& o) { return o; }

  template <typename TableElementType,
            typename TableIteratorType>
  struct iterator_wrappers
  {
    typedef TableIteratorType w_t;

    static TableElementType
    next(w_t& o)
    {
      TableElementType result = o.next();
      if (!result.is_valid()) {
        PyErr_SetString(PyExc_StopIteration, "At end of table.");
        boost::python::throw_error_already_set();
      }
      return result;
    }

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name)
        .def("next", next)
        .def("__iter__", pass_through)
      ;
    }
  };

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_ITERATOR_WRAPPERS_H
