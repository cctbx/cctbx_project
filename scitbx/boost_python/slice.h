#ifndef SCITBX_BOOST_PYTHON_SLICE_H
#define SCITBX_BOOST_PYTHON_SLICE_H

#include <cstddef>
#include <boost/python/slice.hpp>

namespace scitbx { namespace boost_python {

  struct adapted_slice
  {
    long start, stop, step;
    std::size_t size;

    adapted_slice(boost::python::slice const& sl, std::size_t sz);
  };

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_SLICE_H
