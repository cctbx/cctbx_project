#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_C_GRID_FLEX_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_C_GRID_FLEX_CONVERSIONS_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid_periodic.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/to_python_converter.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType,
            typename CgridType>
  struct versa_c_grid_to_flex
  {
    static PyObject* convert(versa<ElementType, CgridType> const& a)
    {
      versa<ElementType, flex_grid<> > result(
        a,
        a.accessor().as_flex_grid());
      return boost::python::incref(boost::python::object(result).ptr());
    }

    static const PyTypeObject* get_pytype()
    {
      return boost::python::converter::registered<
        versa<ElementType, flex_grid<> > >::converters.to_python_target_type();
    }
  };

  template <typename RefCGridType>
  struct ref_c_grid_from_flex
  {
    typedef typename RefCGridType::value_type element_type;
    typedef typename RefCGridType::accessor_type c_grid_type;
    typedef versa<element_type, flex_grid<> > flex_type;

    ref_c_grid_from_flex()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefCGridType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::expected_pytype_for_arg<
          flex_type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      boost::python::extract<flex_type&> flex_proxy(obj);
      if (!flex_proxy.check()) return 0;
      flex_type& a = flex_proxy();
      try { c_grid_type(a.accessor()); }
      catch (...) { return 0; }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      flex_type& a = boost::python::extract<flex_type&>(obj)();
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      c_grid_type c_grid(a.accessor());
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefCGridType>*)
          data)->storage.bytes;
      new (storage) RefCGridType(a.begin(), c_grid);
      data->convertible = storage;
    }
  };

  template <typename ElementType,
            typename CGridType>
  struct c_grid_flex_conversions
  {
    c_grid_flex_conversions()
    {
      boost::python::to_python_converter<
        versa<ElementType, CGridType>,
        versa_c_grid_to_flex<ElementType, CGridType>
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
#endif
        >();
      ref_c_grid_from_flex<const_ref<ElementType, CGridType> >();
      ref_c_grid_from_flex<ref<ElementType, CGridType> >();
    }
  };

  template <typename ElementType>
  struct default_c_grid_flex_conversions
  {
    default_c_grid_flex_conversions()
    {
      c_grid_flex_conversions<ElementType, c_grid<2> >();
      c_grid_flex_conversions<ElementType, c_grid<3> >();
      c_grid_flex_conversions<ElementType, c_grid_padded<2> >();
      c_grid_flex_conversions<ElementType, c_grid_padded<3> >();
      c_grid_flex_conversions<ElementType, c_grid_periodic<3> >();
      c_grid_flex_conversions<ElementType, c_grid_padded_periodic<3> >();
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_C_GRID_FLEX_CONVERSIONS_H
