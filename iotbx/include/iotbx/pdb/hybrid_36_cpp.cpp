#include <iotbx/pdb/hybrid_36_cpp.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <stdexcept>
#include <cstring>

namespace iotbx { namespace pdb { namespace hybrid_36 {

  std::string
  encode(unsigned width, int value)
  {
    char result[16];
    result[0] = ' ';
    const char* errmsg = hy36encode(width, value, result);
    if (errmsg) {
      for(unsigned i=0;i<width;i++) {
        if (result[i] != '*') {
          throw std::runtime_error("internal error: result not reset.");
        }
      }
      if (std::strcmp(errmsg, "value out of range.") == 0) {
        throw std::invalid_argument(errmsg);
      }
      throw std::runtime_error(errmsg);
    }
    return std::string(result);
  }

  int
  decode(unsigned width, const char* s, unsigned s_size)
  {
    int result = -1;
    const char* errmsg = hy36decode(width, s, s_size, &result);
    if (errmsg) {
      if (result != 0) {
        throw std::runtime_error("internal error: result not reset.");
      }
      if (std::strcmp(errmsg, "invalid number literal.") == 0) {
        throw std::invalid_argument(errmsg);
      }
      throw std::runtime_error(errmsg);
    }
    return result;
  }

  int
  decode(unsigned width, std::string const& s)
  {
    return decode(width, s.c_str(), static_cast<unsigned>(s.size()));
  }

}}} // namespace iotbx::pdb::hybrid_36
