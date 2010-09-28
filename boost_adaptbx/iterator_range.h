#ifndef BOOST_ADAPTBX_ITERATOR_RANGE_H
#define BOOST_ADAPTBX_ITERATOR_RANGE_H

#include <boost/python/return_internal_reference.hpp>

namespace boost_adaptbx {

  template <class IteratorRangeType>
  struct iterator_range_wrapper
  {
    typedef IteratorRangeType wt;

    static wt identity(wt r) { return r; }

    static typename wt::value_type next(wt &r) {
      if (r.begin() == r.end()) {
        PyErr_SetString(PyExc_StopIteration, "Exhausted range");
        boost::python::throw_error_already_set();
      }
      typename wt::value_type result = *r.begin();
      r.advance_begin(1);
      return result;
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def("__iter__", identity)
        .def("next", next, return_internal_reference<>())
        ;
    }
  };


}

#endif // GUARD
