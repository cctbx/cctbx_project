#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_REF_FLEX_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_REF_FLEX_CONVERSIONS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace af { namespace boost_python {

  struct trivial_size_functor
  {
    static std::size_t get(std::size_t sz) {
      return sz;
    }
  };

  struct packed_u_size_functor
  {
    static std::size_t get(std::size_t sz) {
      return dimension_from_packed_size(sz);
    }
  };

  template <typename RefType, class SizeFunctor=trivial_size_functor>
  struct ref_from_flex
  {
    typedef typename RefType::value_type element_type;
    typedef versa<element_type, flex_grid<> > flex_type;

    ref_from_flex()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::expected_pytype_for_arg<
          flex_type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object none;
      if (obj_ptr != none.ptr()) {
        boost::python::object obj(boost::python::borrowed(obj_ptr));
        boost::python::extract<flex_type&> flex_proxy(obj);
        if (!flex_proxy.check()) return 0;
        flex_type& a = flex_proxy();
        if (!a.accessor().is_trivial_1d()) return 0;
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
        flex_type& a = boost::python::extract<flex_type&>(obj)();
        if (!a.check_shared_size()) raise_shared_size_mismatch();
        assert(a.accessor().is_trivial_1d());
        bg = a.begin();
        sz = SizeFunctor::get(a.size());
      }
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefType>*)
          data)->storage.bytes;
      new (storage) RefType(bg, sz);
      data->convertible = storage;
    }
  };


  template <typename RefType>
  struct ref_flex_grid_from_flex
  {
    typedef typename RefType::value_type element_type;
    typedef versa<element_type, flex_grid<> > flex_type;

    ref_flex_grid_from_flex()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefType>()
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
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      flex_type& a = boost::python::extract<flex_type&>(obj)();
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefType>*)
          data)->storage.bytes;
      new (storage) RefType(a.begin(), a.accessor());
      data->convertible = storage;
    }
  };

  template <typename ElementType>
  struct ref_flex_conversions
  {
    ref_flex_conversions()
    {
      ref_from_flex<const_ref<ElementType> >();
      ref_from_flex<ref<ElementType> >();
      ref_flex_grid_from_flex<const_ref<ElementType, flex_grid<> > >();
      ref_flex_grid_from_flex<ref<ElementType, flex_grid<> > >();
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_REF_FLEX_CONVERSIONS_H
