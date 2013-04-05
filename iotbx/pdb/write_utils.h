#ifndef IOTBX_PDB_WRITE_UTILS_H
#define IOTBX_PDB_WRITE_UTILS_H

#include <iotbx/pdb/small_str.h>
#include <iotbx/pdb/namespace.h>
#include <boost/noncopyable.hpp>
#include <fstream>
#include <string>

namespace iotbx { namespace pdb { namespace write_utils {

  void
  model_record(
    stream_write& write,
    std::string const& model_id);

  struct fstream_write : stream_write
  {
    std::ofstream* stream;

    fstream_write(std::ofstream* stream_) : stream(stream_) {}

    virtual void
    operator()(const char* s, unsigned n)
    {
      stream->write(s, n);
    }
  };

  struct fstream_open_close : boost::noncopyable
  {
    std::string file_name;
    std::ofstream out;

    fstream_open_close(
      const char* file_name_,
      bool open_append);

    ~fstream_open_close();
  };

  struct sstream_write : stream_write
  {
    std::stringstream* stream;

    sstream_write(std::stringstream* stream_) : stream(stream_) {}

    virtual void
    operator()(const char* s, unsigned n)
    {
      stream->write(s, n);
    }
  };


  struct sstream_open_close : boost::noncopyable
  {
    std::stringstream out;

    sstream_open_close();

    ~sstream_open_close();
  };


}}} // namespace iotbx::pdb::write_utils

#endif // IOTBX_PDB_WRITE_UTILS_H
