// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/class_builder.hpp>
#include <cctbx/error.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace af {

  namespace {

    struct getstate_manager
    {
      getstate_manager(std::size_t a_size, std::size_t size_per_element)
      {
        str_capacity = a_size * size_per_element;
        str_obj = PyString_FromStringAndSize(
          0, static_cast<int>(str_capacity + 100));
        str_begin = PyString_AS_STRING(str_obj);
        str_end = str_begin;
        sprintf(str_end, "%lu", static_cast<unsigned long>(a_size));
        while (*str_end) str_end++;
        *str_end++ = ' ';
      };

      void advance()
      {
        while (*str_end) str_end++;
        *str_end++ = ' ';
        cctbx_assert(str_end - str_begin <= str_capacity);
      }

      PyObject* finalize()
      {
        str_capacity = str_end - str_begin;
        cctbx_assert(
          _PyString_Resize(&str_obj, static_cast<int>(str_capacity)) == 0);
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
        cctbx_assert(a_size == 0);
        str_ptr = PyString_AsString(state);
        cctbx_assert(str_ptr != 0);
        cctbx_assert(sscanf(str_ptr, "%lu", &a_capacity) == 1);
        while (*str_ptr != ' ') str_ptr++;
        str_ptr++;
      };

      void advance()
      {
        while (*str_ptr != ' ') str_ptr++;
        str_ptr++;
      }

      void finalize()
      {
        cctbx_assert(*str_ptr == 0);
      }

      char* str_ptr;
      std::size_t a_capacity;
    };

    template <typename ElementType>
    struct num_picklers
    {
      static
      boost::python::ref
      getstate(
        shared<ElementType> const& a,
        std::size_t size_per_element,
        const char* fmt)
      {
        getstate_manager mgr(a.size(), size_per_element);
        for(std::size_t i=0;i<a.size();i++) {
          sprintf(mgr.str_end, fmt, a[i]);
          mgr.advance();
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(
        shared<ElementType>& a,
        boost::python::ref state,
        const char* fmt)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          ElementType val;
          cctbx_assert(sscanf(mgr.str_ptr, fmt, &val) == 1);
          mgr.advance();
          a.push_back(val);
        }
        mgr.finalize();
      }
    };

  } // namespace <anonymous>

  boost::python::ref shared_int_getstate(shared<int> const& a)
  {
    return num_picklers<int>::getstate(a, 12, "%d");
  }

  void shared_int_setstate(shared<int>& a, boost::python::ref state)
  {
    num_picklers<int>::setstate(a, state, "%d");
  }

  boost::python::ref shared_long_getstate(shared<long> const& a)
  {
    return num_picklers<long>::getstate(a, 21, "%ld");
  }

  void shared_long_setstate(shared<long>& a, boost::python::ref state)
  {
    num_picklers<long>::setstate(a, state, "%ld");
  }

  boost::python::ref shared_float_getstate(shared<float> const& a)
  {
    return num_picklers<float>::getstate(a, 12, "%.6g");
  }

  void shared_float_setstate(shared<float>& a, boost::python::ref state)
  {
    num_picklers<float>::setstate(a, state, "%g");
  }

  boost::python::ref shared_double_getstate(shared<double> const& a)
  {
    return num_picklers<double>::getstate(a, 18, "%.12g");
  }

  void shared_double_setstate(shared<double>& a, boost::python::ref state)
  {
    num_picklers<double>::setstate(a, state, "%lg");
  }

}} // namespace cctbx::af
