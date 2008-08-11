#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_VECTOR_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_VECTOR_WRAPPER_H

#include <boost/python/converter/registry.hpp>
#include <scitbx/array_family/ref.h>

namespace scitbx { namespace array_family { namespace boost_python {

  template <typename ArrayType, typename RefType>
  struct ref_from_array
  {
    typedef typename RefType::value_type element_type;

    ref_from_array()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::expected_pytype_for_arg<
          ArrayType>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object none;
      if (obj_ptr != none.ptr()) {
        boost::python::object obj(boost::python::borrowed(obj_ptr));
        boost::python::extract<ArrayType&> array_proxy(obj);
        if (!array_proxy.check()) return 0;
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
        ArrayType& a = boost::python::extract<ArrayType&>(obj)();
        sz = a.size();
        if (sz != 0) bg = &*a.begin();
      }
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefType>*)
          data)->storage.bytes;
      new (storage) RefType(bg, sz);
      data->convertible = storage;
    }
  };

}}} // namespace scitbx::array_family::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_VECTOR_WRAPPER_H
