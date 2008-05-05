#include <iotbx/pdb/write_utils.h>
#include <ctype.h>

namespace iotbx { namespace pdb { namespace write_utils {

  void
  model_record(
    stream_write& write,
    str8 const& model_id)
  {
    write("MODEL", 5U);
    unsigned n = model_id.rstripped_size();
    if (n != 0) {
      write(" ", 1U);
      for(unsigned i=n;i<8U;i++) write(" ", 1U);
      write(model_id.elems, n);
    }
    write("\n", 1U);
  }

  fstream_open_close::fstream_open_close(
    const char* file_name_,
    bool open_append)
  :
    file_name(file_name_ ? file_name_ : "")
  {
    SCITBX_ASSERT(file_name.size() != 0);
    std::ios::openmode mode = std::ios::out | std::ios::binary;
    if (open_append) mode |= std::ios::app;
    out.open(file_name.c_str(), mode);
    if (out.fail()) {
      throw std::runtime_error(
        "Cannot open file for writing: \"" + file_name + "\"");
    }
  }

  fstream_open_close::~fstream_open_close()
  {
    if (out.fail()) {
      throw std::runtime_error(
        "Failure writing to file: \"" + file_name + "\"");
    }
    out.close();
  }

}}} // namespace iotbx::pdb::write_utils
