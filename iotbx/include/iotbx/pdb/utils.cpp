#include <iotbx/pdb/utils.h>

namespace iotbx { namespace pdb { namespace utils {

  boost::int64_t
  base_256_ordinal(const char* s)
  {
    static const boost::int64_t zero = static_cast<boost::int64_t>(
      *reinterpret_cast<const unsigned char*>("0"));

    if (s == 0) return zero;
    while (*s == ' ') s++; // ignore leading spaces
    if (*s == '\0') return zero;
    bool negative;
    if (*s == '-') {
      negative = true;
      s++;
    }
    else {
      negative = false;
    }
    boost::int64_t result = static_cast<boost::int64_t>(
      *reinterpret_cast<const unsigned char*>(s++));
    while (*s) {
      result *= 256;
      result += static_cast<boost::int64_t>(
        *reinterpret_cast<const unsigned char*>(s++));
    }
    if (negative) return -result;
    return result;
  }

}}} // namespace iotbx::pdb::utils
