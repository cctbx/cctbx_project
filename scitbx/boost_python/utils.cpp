/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <scitbx/boost_python/utils.h>

namespace scitbx { namespace boost_python {

  boost::python::object
  cvs_revision(const std::string& revision)
  {
    using namespace boost::python;
    return object(borrowed(object(
      revision.substr(11, revision.size()-11-2)).ptr()));
  }

  boost::python::handle<>
  import_module(const char* module_name)
  {
    using namespace boost::python;
    return handle<>(PyImport_ImportModule(const_cast<char*>(module_name)));
  }

  void raise_index_error()
  {
    PyErr_SetString(PyExc_IndexError, "Index out of range.");
    boost::python::throw_error_already_set();
  }

  boost::python::object
  range(long start, long len, long step)
  {
    return boost::python::object(boost::python::handle<>(
      PyRange_New(start, len, step, 1)));
  }

  boost::python::object
  range(long len)
  {
    return range(0, len);
  }

}} // namespace scitbx::boost_python
