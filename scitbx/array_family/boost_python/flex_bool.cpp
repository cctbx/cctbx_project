#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/numpy_bridge.hpp>
#include <boost_adaptbx/type_id_eq.h>
#include <boost/python/args.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  template <typename UnsignedType>
  struct from_iselection
  {
    static
    flex<bool>::type*
    get(
      std::size_t size,
      af::const_ref<UnsignedType> const& iselection)
    {
      af::shared<bool> result(size, false);
      for(std::size_t i=0;i<iselection.size();i++) {
        SCITBX_ASSERT(iselection[i] < size);
        result[iselection[i]] = true;
      }
      return new flex<bool>::type(result, result.size());
    }
  };

  static const char* type_error_message_2
    = "Type of argument must be a Python bool or flex.bool.";
  static const char* type_error_message_3
    = "Type of argument must be a Python bool, flex.bool, or None.";

  boost::python::object
  eq(flex_bool const& a, boost::python::object const& b)
  {
    if (b.ptr() == boost::python::object().ptr()) { // if b is None
      return boost::python::object(false);
    }
    {
      boost::python::extract<flex_bool> b_proxy(b);
      if (b_proxy.check()) {
        return boost::python::object(a == b_proxy());
      }
    }
    {
      boost::python::extract<bool> b_proxy(b);
      if (b_proxy.check()) {
        return boost::python::object(a == b_proxy());
      }
    }
    PyErr_SetString(PyExc_TypeError, type_error_message_3);
    boost::python::throw_error_already_set();
    return boost::python::object(); // never reached, avoids warning
  }

  boost::python::object
  ne(flex_bool const& a, boost::python::object const& b)
  {
    if (b.ptr() == boost::python::object().ptr()) { // if b is None
      return boost::python::object(true);
    }
    {
      boost::python::extract<flex_bool> b_proxy(b);
      if (b_proxy.check()) {
        return boost::python::object(a != b_proxy());
      }
    }
    {
      boost::python::extract<bool> b_proxy(b);
      if (b_proxy.check()) {
        return boost::python::object(a != b_proxy());
      }
    }
    PyErr_SetString(PyExc_TypeError, type_error_message_3);
    boost::python::throw_error_already_set();
    return boost::python::object(); // never reached, avoids warning
  }

  boost::python::object
  all_eq(flex_bool const& a, boost::python::object const& b)
  {
    if (b.ptr() != boost::python::object().ptr()) { // if b is not None
      {
        boost::python::extract<flex_bool> b_proxy(b);
        if (b_proxy.check()) {
          return boost::python::object(a.all_eq(b_proxy()));
        }
      }
      {
        boost::python::extract<bool> b_proxy(b);
        if (b_proxy.check()) {
          return boost::python::object(a.all_eq(b_proxy()));
        }
      }
    }
    PyErr_SetString(PyExc_TypeError, type_error_message_2);
    boost::python::throw_error_already_set();
    return boost::python::object(); // never reached, avoids warning
  }

  boost::python::object
  all_ne(flex_bool const& a, boost::python::object const& b)
  {
    if (b.ptr() != boost::python::object().ptr()) { // if b is not None
      {
        boost::python::extract<flex_bool> b_proxy(b);
        if (b_proxy.check()) {
          return boost::python::object(a.all_ne(b_proxy()));
        }
      }
      {
        boost::python::extract<bool> b_proxy(b);
        if (b_proxy.check()) {
          return boost::python::object(a.all_ne(b_proxy()));
        }
      }
    }
    PyErr_SetString(PyExc_TypeError, type_error_message_2);
    boost::python::throw_error_already_set();
    return boost::python::object(); // never reached, avoids warning
  }

  flex_bool
  invert_a(flex_bool const& a) { return !a; }

  flex_bool
  and_a_a(flex_bool const& a1, flex_bool const& a2) { return a1 && a2; }

  flex_bool
  or_a_a(flex_bool const& a1, flex_bool const& a2) { return a1 || a2; }

  flex_bool
  iand_a_a(flex_bool a1, flex_bool const& a2)
  {
    if (a1.accessor() != a2.accessor()) {
      raise_incompatible_arrays();
    }
    for(std::size_t i=0;i<a1.size();i++) if(!a2[i]) a1[i] = false;
    return a1;
  }

  flex_bool
  ior_a_a(flex_bool a1, flex_bool const& a2)
  {
    if (a1.accessor() != a2.accessor()) {
      raise_incompatible_arrays();
    }
    for(std::size_t i=0;i<a1.size();i++) if(a2[i]) a1[i] = true;
    return a1;
  }

  flex_bool
  iand_a_s(flex_bool a1, bool a2)
  {
    if (!a2) std::fill(a1.begin(), a1.end(), false);
    return a1;
  }

  flex_bool
  ior_a_s(flex_bool a1, bool a2)
  {
    if (a2) std::fill(a1.begin(), a1.end(), true);
    return a1;
  }

  bool
  exclusive_or(bool lhs, bool rhs)
  {
    return lhs ? !rhs : rhs;
  }

  bool
  is_super_set (flex_bool const& a1, flex_bool const& a2)
  {
    SCITBX_ASSERT(a2.size() == a1.size());
    bool result = true;
    for (std::size_t i = 0; i < a1.size(); i++) {
      if (! (a2[i] == false || a1[i] == a2[i])) {
        result = false;
        break;
      }
    }
    return result;
  }

  static flex_bool
  exclusive_or_a_a(flex_bool const& a1, flex_bool const& a2)
  {
    SCITBX_ASSERT(a2.size() == a1.size());
    flex_bool result(a1.accessor(), af::init_functor_null<bool>());
    bool* res = result.begin();
    bool* res_end = result.end();
    const bool* lhs = a1.begin();
    const bool* rhs = a2.begin();
    while (res != res_end) {
      *res++ = exclusive_or(*lhs++, *rhs++);
    }
    return result;
  }

  af::shared<int>
  as_int(af::const_ref<bool> const& self)
  {
    af::shared<int> result((af::reserve(self.size())));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i] ? 1 : 0);
    }
    return result;
  }

  af::shared<double>
  as_double(af::const_ref<bool> const& self)
  {
    af::shared<double> result((af::reserve(self.size())));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i] ? 1 : 0);
    }
    return result;
  }

  af::shared<std::size_t>
  iselection(
    af::const_ref<bool, flex_grid<> > const& a,
    bool test_value=true)
  {
    af::shared<std::size_t> result;
    for(std::size_t i=0;i<a.size();i++) {
      if (a[i] == test_value) result.push_back(i);
    }
    return result;
  }

  template <typename UnsignedType>
  struct union_core
  {
    union_core(
      boost::python::object const& iselection,
      af::ref<bool> result)
    :
      ok(false)
    {
      boost::python::extract<af::const_ref<UnsignedType> > proxy(iselection);
      if (proxy.check()) {
        ok = true;
        af::const_ref<UnsignedType> iselection = proxy();
        for(std::size_t i=0;i<iselection.size();i++) {
          SCITBX_ASSERT(iselection[i] < result.size());
          result[iselection[i]] = true;
        }
      }
    }

    bool ok;
  };

  af::shared<bool>
  union_(
    std::size_t size,
    boost::python::list const& iselections)
  {
    af::shared<bool> result(size, false);
    af::ref<bool> r = result.ref();
    std::size_t n_iselections = boost::python::len(iselections);
    for(std::size_t i=0;i<n_iselections;i++) {
      bool ok = union_core<unsigned>(iselections[i], r).ok;
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
      if (!ok) ok = union_core<std::size_t>(iselections[i], r).ok;
#endif
      if (!ok) {
        throw error("iselections must be arrays of unsigned or size_t.");
      }
    }
    return result;
  }

  template <typename UnsignedType>
  struct intersection_core
  {
    intersection_core(
      boost::python::object const& iselection,
      af::ref<bool> result,
      af::ref<bool> tmp)
    :
      ok(false)
    {
      SCITBX_ASSERT(tmp.size() == result.size());
      boost::python::extract<af::const_ref<UnsignedType> > proxy(iselection);
      if (proxy.check()) {
        ok = true;
        af::const_ref<UnsignedType> iselection = proxy();
        for(std::size_t i=0;i<iselection.size();i++) {
          SCITBX_ASSERT(iselection[i] < result.size());
          tmp[iselection[i]] = true;
        }
        for(std::size_t i=0;i<result.size();i++) {
          if (tmp[i]) tmp[i] = false;
          else        result[i] = false;
        }
      }
    }

    bool ok;
  };

  af::shared<bool>
  intersection(
    std::size_t size,
    boost::python::list const& iselections)
  {
    af::shared<bool> result(size, true);
    af::shared<bool> tmp(size, false);
    af::ref<bool> r = result.ref();
    af::ref<bool> t = tmp.ref();
    std::size_t n_iselections = boost::python::len(iselections);
    for(std::size_t i=0;i<n_iselections;i++) {
      bool ok = intersection_core<unsigned>(iselections[i], r, t).ok;
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
      if (!ok) ok = intersection_core<std::size_t>(iselections[i], r, t).ok;
#endif
      if (!ok) {
        throw error("iselections must be arrays of unsigned or size_t.");
      }
    }
    return result;
  }

  template <typename UnsignedType>
  af::shared<UnsignedType>
  filter_indices(
    af::const_ref<bool> const& self,
    af::const_ref<UnsignedType> const& indices)
  {
    af::shared<UnsignedType> result;
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < self.size());
      if (self[indices[i]]) {
        result.push_back(indices[i]);
      }
    }
    return result;
  }

} // namespace <anonymous>

  void wrap_flex_bool()
  {
    using namespace boost::python;
    using boost::python::arg;

    typedef flex_wrapper<bool> f_w;

    f_w::plain("bool")
      .def_pickle(flex_pickle_single_buffered<bool>())
      .def("__init__", make_constructor(
        &from_iselection<unsigned>::get,
        default_call_policies(),
        (arg("size"), arg("iselection"))))
      .def("__init__", make_constructor(
        &from_iselection<std::size_t>::get,
        default_call_policies(),
        (arg("size"), arg("iselection"))))
      .def("__init__", make_constructor(
        flex_bool_from_numpy_array, default_call_policies()))
      .def("__eq__", eq)
      .def("__ne__", ne)
      .def("__eq__", eq)
      .def("__ne__", ne)
      .def("all_eq", all_eq, (arg("other")))
      .def("all_ne", all_ne, (arg("other")))
      .def("all_eq", all_eq, (arg("other")))
      .def("all_ne", all_ne, (arg("other")))
      .def("__invert__", invert_a)
      .def("__and__", and_a_a)
      .def("__or__", or_a_a)
      .def("__iand__", iand_a_a)
      .def("__ior__", ior_a_a)
      .def("__iand__", iand_a_s)
      .def("__ior__", ior_a_s)
      .def("exclusive_or", exclusive_or_a_a, (arg("other")))
      .def("is_super_set", is_super_set, (arg("other")))
      .def("count", f_w::count, (arg("value")))
      .def("as_int", as_int)
      .def("as_double", as_double)
      .def("iselection", iselection, (arg("test_value")=true))
      .def("filter_indices",
        (af::shared<std::size_t>(*)(
           af::const_ref<bool> const&,
           af::const_ref<std::size_t> const&)) filter_indices, (
             arg("indices")))
      .def("as_numpy_array", flex_bool_as_numpy_array)
    ;
    def("order", f_w::order_a_a, (arg("other")));
    def("union", union_, (arg("size"), arg("iselections")));
    def("intersection", intersection, (arg("size"), arg("iselections")));
    def("first_index", f_w::first_index_a_s, (arg("value")));
    def("last_index", f_w::last_index_a_s, (arg("value")));
  }

}}} // namespace scitbx::af::boost_python
