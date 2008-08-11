#ifndef SCITBX_BOOST_PYTHON_SLICE_H
#define SCITBX_BOOST_PYTHON_SLICE_H

#include <cstddef>

namespace scitbx { namespace boost_python {

  struct slice_item
  {
    bool is_valid;
    long value;
  };

  struct slice
  {
    slice_item start, stop, step;
  };

  struct adapted_slice
  {
    long start, stop, step;
    std::size_t size;

    adapted_slice(slice const& sl, std::size_t sz);
  };

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_SLICE_H
