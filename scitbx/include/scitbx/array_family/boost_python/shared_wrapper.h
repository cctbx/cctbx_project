#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/error.h>

namespace scitbx { namespace af { namespace boost_python {

  using scitbx::boost_python::positive_getitem_index;

  template <typename RefType>
  struct ref_from_shared
  {
    typedef typename RefType::value_type element_type;
    typedef shared<element_type> shared_type;

    ref_from_shared()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefType>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object none;
      if (obj_ptr != none.ptr()) {
        boost::python::object obj(boost::python::borrowed(obj_ptr));
        boost::python::extract<shared_type&> shared_proxy(obj);
        if (!shared_proxy.check()) return 0;
      }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::object none;
      element_type* bg = 0;
      std::size_t sz = 0;
      if (obj_ptr != none.ptr()) {
        boost::python::object obj(boost::python::borrowed(obj_ptr));
        shared_type& a = boost::python::extract<shared_type&>(obj)();
        bg = a.begin();
        sz = a.size();
      }
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefType>*)
          data)->storage.bytes;
      new (storage) RefType(bg, sz);
      data->convertible = storage;
    }
  };

  template <typename ElementType,
            typename GetitemReturnValuePolicy
              = boost::python::return_value_policy<
                  boost::python::copy_non_const_reference> >
  struct shared_wrapper
  {
    typedef shared<ElementType> w_t;
    typedef ElementType e_t;

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
    getitem_1d_slice(w_t const& self, scitbx::boost_python::slice const& slice)
    {
      scitbx::boost_python::adapted_slice a_sl(slice, self.size());
      w_t result((af::reserve(a_sl.size)));
      for(long i=a_sl.start;i!=a_sl.stop;i+=a_sl.step) {
        result.push_back(self[i]);
      }
      return result;
    }

    static void
    delitem_1d_slice(w_t& self, scitbx::boost_python::slice const& slice)
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
      self.insert(self.end(), other.begin(), other.end());
    }

    static void
    wrap(std::string const& python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name.c_str())
        .def(init<w_t const&>())
        .def(init<std::size_t const&, optional<e_t const&> >())
        .def("size", &w_t::size)
        .def("__len__", &w_t::size)
        .def("__getitem__", getitem_1d, GetitemReturnValuePolicy())
        .def("__setitem__", setitem_1d)
        .def("__delitem__", delitem_1d)
        .def("__getitem__", getitem_1d_slice)
        .def("__delitem__", delitem_1d_slice)
        .def("clear", &w_t::clear)
        .def("insert", insert)
        .def("append", &w_t::append)
        .def("extend", extend)
      ;

      scitbx::boost_python::container_conversions::from_python_sequence<
        w_t,
        scitbx::boost_python::container_conversions
          ::variable_capacity_policy>();

      ref_from_shared<const_ref<ElementType> >();
      ref_from_shared<ref<ElementType> >();
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_WRAPPER_H
