/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_BOOST_PYTHON_SLICE_H
#define SCITBX_BOOST_PYTHON_SLICE_H

#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace boost_python {

  struct slice_item
  {
    bool is_valid;
    int value;
  };

  slice_item
  extract_slice_item(boost::python::object const& obj);

  struct slice
  {
    slice_item start, stop, step;
  };

  struct slice_from_python
  {
    slice_from_python();

    static void* convertible(PyObject* obj_ptr);

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data);
  };

  struct adapted_slice
  {
    long start, stop, step;
    std::size_t size;

    adapted_slice(slice const& sl, std::size_t sz);
  };

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_SLICE_H
