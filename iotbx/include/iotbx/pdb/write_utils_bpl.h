#ifndef IOTBX_PDB_WRITE_UTILS_BPL_H
#define IOTBX_PDB_WRITE_UTILS_BPL_H

#include <iotbx/pdb/write_utils.h>
#include <scitbx/boost_python/cstringio_utils.h>

namespace iotbx { namespace pdb { namespace write_utils {

  struct cstringio_write : stream_write
  {
    PyObject* sio;

    cstringio_write(PyObject* sio_)
    :
      sio(sio_)
    {
      if (!PycStringIO_OutputCheck(sio)) {
        throw std::invalid_argument(
          "cstringio argument must be a cStringIO.StringIO instance.");
      }
    }

    virtual void
    operator()(const char* s, unsigned n)
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
  };

}}} // namespace iotbx::pdb::write_utils

#endif // IOTBX_PDB_WRITE_UTILS_BPL_H
