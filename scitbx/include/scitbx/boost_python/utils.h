/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/bpl_utils.h (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_BOOST_PYTHON_UTILS_H
#define SCITBX_BOOST_PYTHON_UTILS_H

#include <boost/python/object.hpp>
#include <string>

namespace scitbx { namespace boost_python {

  boost::python::object
  cvs_revision(const std::string& revision);

  boost::python::handle<>
  import_module(const char* module_name);

  void raise_index_error();

  inline
  std::size_t
  positive_getitem_index(long i, std::size_t size)
  {
    if (i >= 0) {
      if (i >= size) raise_index_error();
      return i;
    }
    if (-i > size) raise_index_error();
    return size + i;
  }

  boost::python::object
  range(long start, long len, long step=1);

  boost::python::object
  range(long len);

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_UTILS_H
