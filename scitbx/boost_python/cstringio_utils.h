#ifndef SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H
#define SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H

#if PY_MAJOR_VERSION < 3

#include <stdexcept>
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

  inline
  void
  cstringio_output_check(PyObject* sio)
  {
    if (!PycStringIO_OutputCheck(sio)) {
      throw std::invalid_argument(
        "cstringio argument must be a cStringIO.StringIO instance.");
    }
  }

  inline
  void
  cstringio_cwrite(PyObject* sio, const char* s, unsigned n)
  {
    PycStringIO->cwrite(
      sio,
#if PY_VERSION_HEX >= 0x02050000
      s,
#else
      const_cast<char*>(s),
#endif
      static_cast<boost::python::ssize_t>(n));
  }

}} // namespace scitbx::boost_python

#endif // python 3
#endif // SCITBX_BOOST_PYTHON_CSTRINGIO_UTILS_H
