#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/tuple.hpp>

#include <scitbx/sparse/vector.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename T>
struct vector_wrapper
{
  typedef vector<T> wt;
  typedef typename wt::index_type index_type;
  typedef typename wt::value_type value_type;
  typedef typename wt::iterator iterator;

  static void setitem(wt& self, index_type i, T x) {
    self[i] = x;
  }

  static value_type getitem(wt& self, index_type i) {
    return self[i];
  }

  struct element_iterator
  {
    iterator cur, end;
    element_iterator(iterator first, iterator last) : cur(first), end(last)
    {}

    boost::python::tuple next() {
      if (cur == end) {
        PyErr_SetNone(PyExc_StopIteration);
        boost::python::throw_error_already_set();
      }
      index_type i = cur.index();
      value_type x = *cur++;
      return boost::python::make_tuple(i,x);
    }

    element_iterator iter() {
      return *this;
    }
  };

  struct element_iterator_wrapper
  {
    typedef element_iterator wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt>("element_iterator", no_init)
        .def("next", &wt::next)
        .def("__iter__", &wt::iter)
        ;
    }
  };

  static element_iterator iter(wt& self) {
    return element_iterator(self.begin(), self.end());
  }

  static void wrap() {
    using namespace boost::python;
    class_<wt>("vector", no_init)
      .def(init<index_type>())
      .add_property("size", &wt::size)
      .def("__setitem__", setitem)
      .def("__getitem__", getitem)
      .def("__iter__", iter)
      .def("deep_copy", &wt::deep_copy)
      .def("sort_indices", &wt::sort_indices, return_self<>())
      .def("permute",
           static_cast<wt& (wt::*)(af::const_ref<index_type> const&)>(
                                                                  &wt::permute),
           return_self<>())
      .def("as_dense_vector", &wt::as_dense_vector)
      .def("is_structurally_zero", &wt::is_structurally_zero)
    ;
  }
};

void wrap_vector() {
  vector_wrapper<double>::element_iterator_wrapper::wrap();
  vector_wrapper<double>::wrap();
}

}}}
