#ifndef IOTBX_PDB_WRITE_UTILS_BPL_H
#define IOTBX_PDB_WRITE_UTILS_BPL_H

#include <iotbx/pdb/write_utils.h>

#if PY_MAJOR_VERSION < 3

#include <scitbx/boost_python/cstringio_utils.h>

namespace iotbx { namespace pdb { namespace write_utils {

  struct cstringio_write : stream_write
  {
    PyObject* sio;

    cstringio_write(PyObject* sio_)
    :
      sio(sio_)
    {
      scitbx::boost_python::cstringio_output_check(sio);
    }

    virtual void
    operator()(const char* s, unsigned n)
    {
      scitbx::boost_python::cstringio_cwrite(sio, s, n);
    }
  };

}}} // namespace iotbx::pdb::write_utils

#endif // python 3
#endif // IOTBX_PDB_WRITE_UTILS_BPL_H
