#ifndef SCITBX_BOOST_PYTHON_ITERATOR_WRAPPERS_H
#define SCITBX_BOOST_PYTHON_ITERATOR_WRAPPERS_H

#include <boost_adaptbx/iterator_range.h>

namespace scitbx { namespace boost_python {

  using boost_adaptbx::pass_through;

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
