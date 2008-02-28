#ifndef SCITBX_BOOST_PYTHON_UTILS_H
#define SCITBX_BOOST_PYTHON_UTILS_H

#include <boost/python/object.hpp>
#include <string>

namespace scitbx { namespace boost_python {

  boost::python::handle<>
  import_module(const char* module_name);

  void raise_index_error();

  boost::python::object
  range(long start, long len, long step=1);

  boost::python::object
  range(long len);

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_UTILS_H
