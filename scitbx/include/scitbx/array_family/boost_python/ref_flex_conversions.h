/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_REF_FLEX_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_REF_FLEX_CONVERSIONS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename RefType>
  struct ref_from_flex
  {
    typedef typename RefType::value_type element_type;
    typedef versa<element_type, flex_grid<> > flex_type;

    ref_from_flex()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<RefType>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      boost::python::extract<flex_type&> flex_proxy(obj);
      if (!flex_proxy.check()) return 0;
      flex_type& a = flex_proxy();
      if (!a.accessor().is_trivial_1d()) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::python::object obj(boost::python::borrowed(obj_ptr));
      flex_type& a = boost::python::extract<flex_type&>(obj)();
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      assert(a.accessor().is_trivial_1d());
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<RefType>*)
          data)->storage.bytes;
      new (storage) RefType(a.begin(), a.size());
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
        boost::python::type_id<RefType>());
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
