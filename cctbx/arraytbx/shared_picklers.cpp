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

    template <typename ElementType>
    struct num_picklers
    {
      static
      boost::python::ref
      getstate(shared<ElementType> const& a, const char* fmt)
      {
        std::size_t str_capacity = a.size() * 18;
        PyObject* str_obj = PyString_FromStringAndSize(
          0, static_cast<int>(str_capacity+100));
        char* str_begin = PyString_AS_STRING(str_obj);
        char* str_end = str_begin;
        sprintf(str_end, "%lu", static_cast<unsigned long>(a.size()));
        while (*str_end) str_end++;
        *str_end++ = ' ';
        for(std::size_t i=0;i<a.size();i++) {
          sprintf(str_end, fmt, a[i]);
          while (*str_end) str_end++;
          *str_end++ = ' ';
          cctbx_assert(str_end - str_begin <= str_capacity);
        }
        str_capacity = str_end - str_begin;
        cctbx_assert(
          _PyString_Resize(&str_obj, static_cast<int>(str_capacity)) == 0);
        return boost::python::ref(str_obj);
      }

      static
      void
      setstate(
        shared<ElementType>& a,
        boost::python::ref state,
        const char* fmt)
      {
        cctbx_assert(a.size() == 0);
        char* str_ptr = PyString_AsString(state.get());
        cctbx_assert(str_ptr != 0);
        std::size_t str_capacity;
        cctbx_assert(sscanf(str_ptr, "%lu", &str_capacity) == 1);
        while (*str_ptr != ' ') str_ptr++;
        str_ptr++;
        a.reserve(str_capacity);
        for(std::size_t i=0;i<str_capacity;i++) {
          ElementType val;
          cctbx_assert(sscanf(str_ptr, fmt, &val) == 1);
          while (*str_ptr != ' ') str_ptr++;
          str_ptr++;
          a.push_back(val);
        }
        cctbx_assert(*str_ptr == 0);
      }
    };

  } // namespace <anonymous>

  boost::python::ref shared_int_getstate(shared<int> const& a)
  {
    return num_picklers<int>::getstate(a, "%d");
  }

  void shared_int_setstate(shared<int>& a, boost::python::ref state)
  {
    num_picklers<int>::setstate(a, state, "%d");
  }

  boost::python::ref shared_long_getstate(shared<long> const& a)
  {
    return num_picklers<long>::getstate(a, "%ld");
  }

  void shared_long_setstate(shared<long>& a, boost::python::ref state)
  {
    num_picklers<long>::setstate(a, state, "%ld");
  }

  boost::python::ref shared_float_getstate(shared<float> const& a)
  {
    return num_picklers<float>::getstate(a, "%.6g");
  }

  void shared_float_setstate(shared<float>& a, boost::python::ref state)
  {
    num_picklers<float>::setstate(a, state, "%g");
  }

  boost::python::ref shared_double_getstate(shared<double> const& a)
  {
    return num_picklers<double>::getstate(a, "%.12g");
  }

  void shared_double_setstate(shared<double>& a, boost::python::ref state)
  {
    num_picklers<double>::setstate(a, state, "%lg");
  }

}} // namespace cctbx::af
