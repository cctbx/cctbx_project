/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/arraytbx/flex_picklers.cpp (rwgk)
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2002 Aug: Created, based on shared_picklers.cpp (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H

#include <scitbx/array_family/type_holder.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace scitbx { namespace af { namespace boost_python {

  namespace detail {

    struct getstate_manager
    {
      getstate_manager(std::size_t a_size, std::size_t size_per_element)
      {
        str_capacity = a_size * size_per_element + 50;// extra space for a_size
        str_obj = PyString_FromStringAndSize(
          0, static_cast<int>(str_capacity + 100)); // extra space for safety
        str_begin = PyString_AS_STRING(str_obj);
        str_end = scitbx::boost_python::pickle_single_buffered::to_string(
          str_begin, a_size);
      };

      void advance(char* str_ptr)
      {
        str_end = str_ptr;
        SCITBX_ASSERT(str_end - str_begin <= str_capacity);
      }

      PyObject* finalize()
      {
        SCITBX_ASSERT(_PyString_Resize(&str_obj,
          static_cast<int>(str_end - str_begin)) == 0);
        return str_obj;
      }

      std::size_t str_capacity;
      PyObject* str_obj;
      char* str_begin;
      char* str_end;
    };

    struct setstate_manager
    {
      setstate_manager(std::size_t a_size, PyObject* state)
      {
        SCITBX_ASSERT(a_size == 0);
        str_ptr = PyString_AsString(state);
        SCITBX_ASSERT(str_ptr != 0);
        a_capacity = get_value(type_holder<std::size_t>());
      };

      template <typename ValueType>
      ValueType get_value(type_holder<ValueType>)
      {
        scitbx::boost_python::pickle_single_buffered::from_string<ValueType>
          proxy(str_ptr);
        str_ptr = proxy.end;
        return proxy.value;
      }

      void assert_end()
      {
        SCITBX_ASSERT(*str_ptr == 0);
      }

      const char* str_ptr;
      std::size_t a_capacity;
    };

  } // namespace detail

  template <typename ElementType, std::size_t SizePerElement>
  struct flex_pickle_single_buffered : boost::python::pickle_suite
  {
    static
    boost::python::tuple
    getstate(versa<ElementType, flex_grid<> > const& a)
    {
      detail::getstate_manager mgr(a.size(), SizePerElement);
      for(std::size_t i=0;i<a.size();i++) {
        mgr.advance(
          scitbx::boost_python::pickle_single_buffered::to_string(
            mgr.str_end, a[i]));
      }
      return boost::python::make_tuple(
        a.accessor(), boost::python::handle<>(mgr.finalize()));
    }

    static
    void
    setstate(
      versa<ElementType, flex_grid<> >& a,
      boost::python::tuple state)
    {
      SCITBX_ASSERT(boost::python::len(state) == 2);
      flex_grid<> a_accessor = boost::python::extract<flex_grid<> >(
        state[0])();
      detail::setstate_manager
        mgr(a.size(), boost::python::object(state[1]).ptr());
      shared_plain<ElementType> b = a.as_base_array();
      b.reserve(mgr.a_capacity);
      for(std::size_t i=0;i<mgr.a_capacity;i++) {
        b.push_back(mgr.get_value(type_holder<ElementType>()));
      }
      mgr.assert_end();
      SCITBX_ASSERT(b.size() == a_accessor.size1d());
      a.resize(a_accessor);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H
