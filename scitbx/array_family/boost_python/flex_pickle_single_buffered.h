#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H

#include <stdint.h>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <scitbx/serialization/single_buffered.h>
#include <scitbx/type_holder.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace af { namespace boost_python {

  namespace detail {

    struct getstate_manager
    {
      getstate_manager(std::size_t a_size, std::size_t size_per_element)
      {
        str_capacity = a_size * size_per_element + 50;// extra space for a_size
#ifdef IS_PY3K
        str_obj = PyBytes_FromStringAndSize(
#else
        str_obj = PyString_FromStringAndSize(
#endif
          0, static_cast<boost::python::ssize_t>(str_capacity + 100));
            // extra space for safety
#ifdef IS_PY3K
        str_begin = PyBytes_AS_STRING(str_obj);
#else
        str_begin = PyString_AS_STRING(str_obj);
#endif
        str_end = scitbx::serialization::single_buffered::to_string(
          str_begin, a_size);
      };

      void advance(char* str_ptr)
      {
        str_end = str_ptr;
        SCITBX_ASSERT(str_end - str_begin <= str_capacity);
      }

      PyObject* finalize()
      {
#ifdef IS_PY3K
        if (_PyBytes_Resize(&str_obj,
#else
        if (_PyString_Resize(&str_obj,
#endif
              static_cast<boost::python::ssize_t>(str_end - str_begin)) != 0) {
          boost::python::throw_error_already_set();
        }
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
#ifdef IS_PY3K
        if(PyUnicode_Check(state)) {
          str_ptr = PyUnicode_AsUTF8(state);
        }
        else {
          SCITBX_ASSERT(PyBytes_Check(state));
          str_ptr = PyBytes_AsString(state);
        }
#else
        str_ptr = PyString_AsString(state);
#endif
        SCITBX_ASSERT(str_ptr != 0);
        a_capacity = get_value(type_holder<std::size_t>());
      };

      template <typename ValueType>
      ValueType get_value(type_holder<ValueType>)
      {
        scitbx::serialization::single_buffered::from_string<ValueType>
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

  template <typename T>
  struct pickle_size_per_element;

  template <>
  struct pickle_size_per_element<bool>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = 1);
  };

  template <>
  struct pickle_size_per_element<int>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(int)+1));
  };

  template <>
  struct pickle_size_per_element<unsigned int>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(unsigned int)+1));
  };

  template <>
  struct pickle_size_per_element<long>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(long)+1));
  };

  template <>
  struct pickle_size_per_element<unsigned long>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(unsigned long)+1));
  };

  template <>
  struct pickle_size_per_element<long long>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(long long)+1));
  };

  template <>
  struct pickle_size_per_element<unsigned long long>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(unsigned long long)+1));
  };

  template <>
  struct pickle_size_per_element<float>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(float)+3));
  };

  template <>
  struct pickle_size_per_element<double>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(double)+3));
  };

  template <>
  struct pickle_size_per_element<std::complex<float> >
  {
    BOOST_STATIC_CONSTANT(std::size_t,
      value = (2*pickle_size_per_element<float>::value));
  };

  template <>
  struct pickle_size_per_element<std::complex<double> >
  {
    BOOST_STATIC_CONSTANT(std::size_t,
      value = (2*pickle_size_per_element<double>::value));
  };

  template <>
  struct pickle_size_per_element<int16_t>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(int16_t)+1));
  };

  template <>
  struct pickle_size_per_element<uint16_t>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(uint16_t)+1));
  };

  template <>
  struct pickle_size_per_element<int8_t>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(int8_t)+1));
  };

  template <>
  struct pickle_size_per_element<uint8_t>
  {
    BOOST_STATIC_CONSTANT(std::size_t, value = (sizeof(uint8_t)+1));
  };

  template <typename ElementType,
            std::size_t SizePerElement
              = pickle_size_per_element<ElementType>::value>
  struct flex_pickle_single_buffered : boost::python::pickle_suite
  {
    static
    boost::python::tuple
    getstate(versa<ElementType, flex_grid<> > const& a)
    {
      detail::getstate_manager mgr(a.size(), SizePerElement);
      for(std::size_t i=0;i<a.size();i++) {
        mgr.advance(
          scitbx::serialization::single_buffered::to_string(
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
      SCITBX_ASSERT(b.size() == a_accessor.size_1d());
      a.resize(a_accessor);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_SINGLE_BUFFERED_H
