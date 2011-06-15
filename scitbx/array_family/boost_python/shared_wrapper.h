#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/slice.hpp>
#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/ref_from_array.h>
#include <scitbx/misc/positive_getitem_index.h>
#include <scitbx/error.h>

namespace scitbx { namespace af { namespace boost_python {

  using scitbx::positive_getitem_index;

  template <typename ElementType>
  struct shared_wrapper_default_element
  {
    static ElementType
    get() { return ElementType(); }
  };

  template <typename ElementType,
            typename GetitemReturnValuePolicy
              = boost::python::return_value_policy<
                  boost::python::copy_non_const_reference> >
  struct shared_wrapper
  {
    typedef shared<ElementType> w_t;
    typedef ElementType e_t;

    static w_t*
    init_with_default_value(
      std::size_t size)
    {
      return new w_t(size, shared_wrapper_default_element<e_t>::get());
    }

    static e_t&
    getitem_1d(w_t& self, long i)
    {
      return self[positive_getitem_index(i, self.size())];
    }

    static void
    setitem_1d(w_t& self, long i, e_t const& x)
    {
      self[positive_getitem_index(i, self.size())] = x;
    }

    static void
    delitem_1d(w_t& self, long i)
    {
      self.erase(&self[positive_getitem_index(i, self.size())]);
    }

    static w_t
    getitem_1d_slice(w_t const& self, boost::python::slice const& slice)
    {
      scitbx::boost_python::adapted_slice a_sl(slice, self.size());
      w_t result((af::reserve(a_sl.size)));
      for(long i=a_sl.start;i!=a_sl.stop;i+=a_sl.step) {
        result.push_back(self[i]);
      }
      return result;
    }

    static void
    delitem_1d_slice(w_t& self, boost::python::slice const& slice)
    {
      scitbx::boost_python::adapted_slice a_sl(slice, self.size());
      SCITBX_ASSERT(a_sl.step == 1);
      self.erase(&self[a_sl.start], &self[a_sl.stop]);
    }

    static void
    insert(w_t& self, long i, e_t const& x)
    {
      self.insert(&self[positive_getitem_index(i, self.size())], x);
    }

    static void
    extend(w_t& self, w_t const& other)
    {
      self.extend(other.begin(), other.end());
    }

    static void
    reserve(w_t& self, std::size_t size)
    {
      self.reserve(size);
    }

    static
    boost::python::class_<w_t>
    wrap(std::string const& python_name)
    {
      using namespace boost::python;
      class_<w_t> result(python_name.c_str());
      result
        .def(init<w_t const&>())
        .def(init<std::size_t const&, e_t const&>((
          boost::python::arg("size"),
          boost::python::arg("value"))))
        .def("__init__", make_constructor(
          init_with_default_value,
            default_call_policies(), (
              boost::python::arg("size"))))
        .def("size", &w_t::size)
        .def("__len__", &w_t::size)
        .def("__getitem__", getitem_1d, GetitemReturnValuePolicy())
        .def("__setitem__", setitem_1d)
        .def("__delitem__", delitem_1d)
        .def("__getitem__", getitem_1d_slice)
        .def("__delitem__", delitem_1d_slice)
        .def("deep_copy", &w_t::deep_copy)
        .def("clear", &w_t::clear)
        .def("insert", insert)
        .def("append", &w_t::push_back)
        .def("extend", extend)
        .def("reserve", reserve)
      ;

      scitbx::boost_python::container_conversions::from_python_sequence<
        w_t,
        scitbx::boost_python::container_conversions
          ::variable_capacity_policy>();

      using scitbx::array_family::boost_python::ref_from_array;
      ref_from_array<w_t, const_ref<ElementType> >();
      ref_from_array<w_t, ref<ElementType> >();

      return result;
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H
