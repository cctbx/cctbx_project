/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType>
  struct shared_to_flex
  {
    static PyObject* convert(shared_plain<ElementType> const& a)
    {
      using namespace boost::python;
      using boost::python::incref; // works around gcc 2.96 bug
      versa<ElementType, flex_grid<> > result(a, flex_grid<>(a.size()));
      return incref(object(result).ptr());
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
        boost::python::type_id<SharedType>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      using namespace boost::python;
      using boost::python::borrowed; // works around gcc 2.96 bug
      object obj(borrowed(obj_ptr));
      extract<flex_type const&> flex_proxy(obj);
      if (!flex_proxy.check()) return 0;
      flex_type const& a = flex_proxy();
      if (a.accessor().nd() != 1) return 0;
      if (!a.accessor().is_0_based()) return 0;
      if (a.accessor().is_padded()) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      using namespace boost::python;
      using boost::python::borrowed; // works around gcc 2.96 bug
      using boost::python::converter::rvalue_from_python_storage; // dito
      object obj(borrowed(obj_ptr));
      flex_type const& a = extract<flex_type const&>(obj)();
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      assert(a.accessor().nd() == 1);
      assert(a.accessor().is_0_based());
      assert(!a.accessor().is_padded());
      void* storage = (
        (rvalue_from_python_storage<SharedType>*)
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
        shared_to_flex<ElementType> >();
      boost::python::to_python_converter<
        shared<ElementType>,
        shared_to_flex<ElementType> >();
      shared_from_flex<shared_plain<ElementType> >();
      shared_from_flex<shared<ElementType> >();
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLEX_CONVERSIONS_H
