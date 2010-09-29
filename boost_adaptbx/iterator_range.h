#ifndef BOOST_ADAPTBX_ITERATOR_RANGE_H
#define BOOST_ADAPTBX_ITERATOR_RANGE_H

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace boost_adaptbx {

  inline
  boost::python::object
  pass_through(boost::python::object const& o) { return o; }

  template <class IteratorRangeType>
  struct iterator_range_wrapper
  {
    typedef IteratorRangeType wt;

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
        .def("__iter__", pass_through)
        .def("next", next, return_internal_reference<>())
        ;
    }
  };


  template <class IteratorRangeType>
  boost::python::tuple iterator_range_as_tuple(IteratorRangeType const &p) {
    using namespace boost::python;
    boost::python::tuple result;
    for (typename IteratorRangeType::iterator i=p.begin(); i!=p.end(); ++i) {
      typename IteratorRangeType::value_type x = *i;
      result += boost::python::make_tuple(x);
    }
    return result;
  }
}

#endif // GUARD
