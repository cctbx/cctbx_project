#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H

#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <scitbx/serialization/double_buffered.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType,
            typename DoubleBufferedToString
              = scitbx::serialization::double_buffered::to_string,
            typename DoubleBufferedFromString
              = scitbx::serialization::double_buffered::from_string>
  struct flex_pickle_double_buffered : boost::python::pickle_suite
  {
    static
    boost::python::tuple
    getstate(versa<ElementType, flex_grid<> > const& a)
    {
      DoubleBufferedToString accu;
      accu << a.size();
      for(std::size_t i=0;i<a.size();i++) accu << a[i];
#ifdef IS_PY3K
      return boost::python::make_tuple(a.accessor(),
        boost::python::handle<>(PyBytes_FromStringAndSize(accu.buffer.c_str(), accu.buffer.size())));
#else
      return boost::python::make_tuple(a.accessor(), accu.buffer);
#endif
    }

    static
    void
    setstate(versa<ElementType, flex_grid<> >& a, boost::python::tuple state)
    {
      SCITBX_ASSERT(boost::python::len(state) == 2);
      SCITBX_ASSERT(a.size() == 0);
      flex_grid<> a_accessor = boost::python::extract<flex_grid<> >(
        state[0])();
      PyObject* py_str = boost::python::object(state[1]).ptr();
#ifdef IS_PY3K
      DoubleBufferedFromString inp(PyBytes_AsString(py_str));
#else
      DoubleBufferedFromString inp(PyString_AsString(py_str));
#endif
      std::size_t a_capacity;
      inp >> a_capacity;
      shared_plain<ElementType> b = a.as_base_array();
      b.reserve(a_capacity);
      ElementType val(flex_default_element<ElementType>::get());
      for(std::size_t i=0;i<a_capacity;i++) {
        inp >> val;
        b.push_back(val);
      }
      inp.assert_end();
      SCITBX_ASSERT(b.size() == a_accessor.size_1d());
      a.resize(a_accessor, flex_default_element<ElementType>::get());
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H
