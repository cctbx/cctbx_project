#ifndef SCITBX_BOOST_PYTHON_STL_SET_WRAPPER_H
#define SCITBX_BOOST_PYTHON_STL_SET_WRAPPER_H

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/boost_python/utils.h>
#include <set>

namespace scitbx { namespace stl { namespace boost_python {

  template <typename ElementType>
  struct set_wrapper
  {
    typedef std::set<ElementType> w_t;
    typedef ElementType e_t;

    static void
    insert_element(w_t& self, e_t const& x) { self.insert(x); }

    static void
    insert_set(w_t& self, w_t const& other)
    {
      self.insert(other.begin(), other.end());
    }

    static bool
    contains(w_t const& self, e_t const& x)
    {
      return self.find(x) != self.end();
    }

    static e_t
    getitem(w_t const& self, std::size_t i)
    {
      if (i >= self.size()) scitbx::boost_python::raise_index_error();
      typename w_t::const_iterator p = self.begin();
      while (i > 0) { p++; i--; }
      return *p;
    }

    static boost::python::tuple
    getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(boost::python::tuple(self));
    }

    static void
    wrap(std::string const& python_name)
    {
      namespace bp = boost::python;
      bp::class_<w_t, std::auto_ptr<w_t> >(python_name.c_str())
        .def(bp::init<w_t const&>())
        .def("size", &w_t::size)
        .def("__len__", &w_t::size)
        .def("insert", insert_element)
        .def("append", insert_element)
        .def("insert", insert_set)
        .def("extend", insert_set)
        .def("erase", (std::size_t(w_t::*)(e_t const&)) &w_t::erase)
        .def("clear", &w_t::clear)
        .def("__contains__", contains)
        .def("__getitem__", getitem)
        .enable_pickling()
        .def("__getinitargs__", getinitargs)
      ;

      scitbx::boost_python::container_conversions::from_python_sequence<
        w_t, scitbx::boost_python::container_conversions::set_policy>();
    }
  };

}}} // namespace scitbx::stl::boost_python

#endif // SCITBX_BOOST_PYTHON_STL_SET_WRAPPER_H
