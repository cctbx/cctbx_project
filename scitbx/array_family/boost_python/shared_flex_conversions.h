#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename SharedType>
  struct shared_to_flex
  {
    static PyObject* convert(SharedType const& a)
    {
      typedef typename SharedType::value_type value_type;
      versa<value_type, flex_grid<> > result(a, flex_grid<>(a.size()));
      return boost::python::incref(boost::python::object(result).ptr());
    }

    static const PyTypeObject* get_pytype()
    {
      typedef typename SharedType::value_type value_type;
      return boost::python::converter::registered<
        versa<value_type, flex_grid<> > >::converters.to_python_target_type();
    }
  };

  template <typename SharedType>
  struct shared_from_flex
  {
    typedef typename SharedType::value_type element_type;
    typedef versa<element_type, flex_grid<> > flex_type;

    shared_from_flex()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<SharedType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::expected_pytype_for_arg<
          flex_type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      boost::python::extract<flex_type const&> flex_proxy(obj);
      if (!flex_proxy.check()) return 0;
      flex_type const& a = flex_proxy();
      if (!a.accessor().is_trivial_1d()) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      flex_type const& a = boost::python::extract<flex_type const&>(obj)();
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      assert(a.accessor().is_trivial_1d());
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<SharedType>*)
          data)->storage.bytes;
      new (storage) SharedType(a);
      data->convertible = storage;
    }
  };

  template <typename ElementType>
  struct shared_flex_conversions
  {
    shared_flex_conversions()
    {
      boost::python::to_python_converter<
        shared_plain<ElementType>,
        shared_to_flex<shared_plain<ElementType> >
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
#endif
        >();
      boost::python::to_python_converter<
        shared<ElementType>,
        shared_to_flex<shared<ElementType> >
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
#endif
        >();
      shared_from_flex<shared_plain<ElementType> >();
      shared_from_flex<shared<ElementType> >();
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H
