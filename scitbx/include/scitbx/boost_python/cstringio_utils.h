#ifndef SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H
#define SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H

#include <cStringIO.h>

namespace scitbx { namespace boost_python {

  inline
  void
  cstringio_import()
  {
    const char* module_name = "cStringIO";
    const char* cobject_name = "cStringIO_CAPI";
    char* m = const_cast<char*>(module_name);
    char* c = const_cast<char*>(cobject_name);
    PycStringIO = reinterpret_cast<PycStringIO_CAPI*>(PyCObject_Import(m, c));
  }

}} // namespace scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H
