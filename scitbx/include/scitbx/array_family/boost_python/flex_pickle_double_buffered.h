/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H

#include <scitbx/boost_python/pickle_double_buffered.h>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType,
            typename DoubleBufferedToString
              = scitbx::boost_python::pickle_double_buffered::to_string,
            typename DoubleBufferedFromString
              = scitbx::boost_python::pickle_double_buffered::from_string>
  struct flex_pickle_double_buffered : boost::python::pickle_suite
  {
    static
    boost::python::tuple
    getstate(versa<ElementType, flex_grid<> > const& a)
    {
      DoubleBufferedToString accu;
      accu << a.size();
      for(std::size_t i=0;i<a.size();i++) accu << a[i];
      return boost::python::make_tuple(a.accessor(),
#if 0
        accu.buffer
#else
        // XXX work around bug in boost 1.29.0
        boost::python::object(boost::python::handle<>(
          PyString_FromStringAndSize(
            accu.buffer.begin(),
            static_cast<int>(accu.buffer.size()))))
#endif
        );
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
      DoubleBufferedFromString inp(py_str);
      std::size_t a_capacity;
      inp >> a_capacity;
      shared_plain<ElementType> b = a.as_base_array();
      b.reserve(a_capacity);
      ElementType val;
      for(std::size_t i=0;i<a_capacity;i++) {
        inp >> val;
        b.push_back(val);
      }
      inp.assert_end();
      SCITBX_ASSERT(b.size() == a_accessor.size_1d());
      a.resize(a_accessor);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_PICKLE_DOUBLE_BUFFERED_H
