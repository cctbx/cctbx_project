#ifndef SCITBX_BOOST_PYTHON_STL_MAP_WRAPPER_H
#define SCITBX_BOOST_PYTHON_STL_MAP_WRAPPER_H

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace scitbx { namespace stl { namespace boost_python {

  template <typename MapType>
  struct from_python_dict
  {
    typedef typename MapType::key_type k_t;
    typedef typename MapType::mapped_type m_t;

    from_python_dict()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<MapType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::wrap_pytype<&PyDict_Type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      return PyDict_Check(obj_ptr) ? obj_ptr : 0;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::handle<> obj_hdl(boost::python::borrowed(obj_ptr));
      boost::python::object obj_obj(obj_hdl);
      boost::python::extract<boost::python::dict> obj_proxy(obj_obj);
      boost::python::dict other = obj_proxy();
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<MapType>*)
          data)->storage.bytes;
      new (storage) MapType();
      data->convertible = storage;
      MapType& self = *((MapType*)storage);
      boost::python::list keys = other.keys();
      int len_keys = boost::python::len(keys);
      for(int i=0;i<len_keys;i++) {
        boost::python::object key_obj = keys[i];
        boost::python::extract<k_t> key_proxy(key_obj);
        if (!key_proxy.check()) {
          PyErr_SetString(PyExc_KeyError, "Unsuitable type.");
          boost::python::throw_error_already_set();
        }
        boost::python::object value_obj = other[key_obj];
        boost::python::extract<m_t> value_proxy(value_obj);
        if (!value_proxy.check()) {
          PyErr_SetString(PyExc_ValueError, "Unsuitable type.");
          boost::python::throw_error_already_set();
        }
        k_t key = key_proxy();
        m_t value = value_proxy();
        self[key] = value;
      }
    }
  };

  template <typename MapType,
            typename GetitemReturnValuePolicy
              = boost::python::return_value_policy<
                  boost::python::copy_non_const_reference> >
  struct map_wrapper
  {
    typedef MapType w_t;
    typedef typename MapType::key_type k_t;
    typedef typename MapType::mapped_type m_t;

    static bool
    contains(w_t const& self, k_t const& k)
    {
      return self.find(k) != self.end();
    }

    static boost::python::object
    get(
      boost::python::object self,
      boost::python::object k,
      boost::python::object d)
    {
      namespace bp = boost::python;
      bp::extract<w_t const&> self_proxy(self);
      w_t const& c_self = self_proxy();
      bp::extract<k_t const&> k_const_proxy(k);
      if (k_const_proxy.check()) {
        k_t c_k = k_const_proxy();
        if (c_self.find(c_k) == c_self.end()) return d;
      }
      bp::extract<k_t> k_value_proxy(k);
      k_t c_k = k_value_proxy();
      if (c_self.find(c_k) == c_self.end()) return d;
      return self[k]; // call through Python to use correct return value policy
    }

    static m_t&
    setdefault_2(w_t& self, k_t const& k, m_t const& d)
    {
      if (self.find(k) == self.end()) {
        self[k] = d;
      }
      return self[k];
    }

    static m_t&
    setdefault_1(w_t& self, k_t const& k)
    {
      if (self.find(k) == self.end()) {
        self[k];
      }
      return self[k];
    }

    static m_t&
    getitem(w_t& self, k_t const& k)
    {
      if (self.find(k)==self.end()) {
        PyErr_SetString(PyExc_KeyError, "Key not in C++ map.");
        boost::python::throw_error_already_set();
      }
      return self[k];
    }

    static void
    setitem(w_t& self, k_t const& k, m_t const& x) { self[k] = x; }

    static void
    delitem(w_t& self, k_t const& k)
    {
      typename w_t::iterator pos = self.find(k);
      if (pos == self.end()) {
        PyErr_SetString(PyExc_KeyError, "Key not in C++ map.");
        boost::python::throw_error_already_set();
      }
      self.erase(pos);
    }

    static boost::python::list
    keys(w_t const& self)
    {
      boost::python::list result;
      typename w_t::const_iterator i;
      for(i=self.begin();i!=self.end();i++) {
        result.append(i->first);
      }
      return result;
    }

    static boost::python::list
    values(
      boost::python::object self)
    {
      namespace bp = boost::python;
      bp::list result;
      bp::extract<w_t const&> self_proxy(self);
      w_t const& c_self = self_proxy();
      typename w_t::const_iterator i;
      for(i=c_self.begin();i!=c_self.end();i++) {
        result.append(
          // call through Python to use correct return value policy
          self[i->first]);
      }
      return result;
    }

    static boost::python::list
    items(
      boost::python::object self)
    {
      namespace bp = boost::python;
      bp::list result;
      bp::extract<w_t const&> self_proxy(self);
      w_t const& c_self = self_proxy();
      typename w_t::const_iterator i;
      for(i=c_self.begin();i!=c_self.end();i++) {
        result.append(bp::make_tuple(
          i->first,
          // call through Python to use correct return value policy
          self[i->first]));
      }
      return result;
    }

    static void
    update(w_t& self, w_t const& other)
    {
      typename w_t::const_iterator i;
      for(i=other.begin();i!=other.end();i++) {
        self[i->first] = i->second;
      }
    }

    static boost::python::tuple
    popitem(w_t& self)
    {
      typename w_t::iterator i = self.begin();
      if (i == self.end()) {
        PyErr_SetString(PyExc_KeyError, "popitem(): C++ map is empty");
        boost::python::throw_error_already_set();
      }
      boost::python::tuple result = boost::python::make_tuple(
        i->first,
        i->second); // value is copied; cannot use return value policy here
      self.erase(i);
      return result;
    }

    static boost::python::object
    iter(w_t const& self)
    {
      boost::python::handle<> key_iter(PyObject_GetIter(keys(self).ptr()));
      return boost::python::object(key_iter);
    }

    static boost::python::tuple
    getinitargs(
      boost::python::object self)
    {
      namespace bp = boost::python;
      return bp::make_tuple(bp::dict(items(self)));
    }

    static void
    wrap(std::string const& python_name)
    {
      using namespace boost::python;
      class_<w_t, boost::shared_ptr<w_t> >(python_name.c_str())
        .def(init<w_t const&>())
        .def("size", &w_t::size)
        .def("__len__", &w_t::size)
        .def("erase", (std::size_t(w_t::*)(k_t const&)) &w_t::erase)
        .def("clear", &w_t::clear)
        .def("__contains__", contains)
        .def("has_key", contains)
        .def("get", get, (arg("k"), arg("d")=object()))
        .def("setdefault", setdefault_2, GetitemReturnValuePolicy())
        .def("setdefault", setdefault_1, GetitemReturnValuePolicy())
        .def("__getitem__", getitem, GetitemReturnValuePolicy())
        .def("__setitem__", setitem)
        .def("__delitem__", delitem)
        .def("keys", keys)
        .def("values", values)
        .def("items", items)
        .def("update", update)
        .def("popitem", popitem)
        .def("__iter__", iter)
        .enable_pickling()
        .def("__getinitargs__", getinitargs)
      ;

      from_python_dict<w_t>();
    }
  };

}}} // namespace scitbx::stl::boost_python

#endif // SCITBX_BOOST_PYTHON_STL_MAP_WRAPPER_H
