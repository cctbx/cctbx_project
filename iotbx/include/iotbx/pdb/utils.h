#ifndef IOTBX_PDB_UTILS_H
#define IOTBX_PDB_UTILS_H

namespace iotbx { namespace pdb {

//! Miscellaneous algorithms.
namespace utils {

  //! Integer equivalent of base-256 pseudo-numbers.
  /*! Similar to interpretation of base-10 numbers.
      Leading spaces are ignored.
      If the first non-space is a minus sign, it is interpreted as such.
      All following characters are interpreted as base-256 digits,
      spaces and minus signs included.
      NOTE: trailing spaces are significant!
   */
  int
  base_256_ordinal(const char* s)
  {
    static const int zero = static_cast<int>(
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
    int result = static_cast<int>(
      *reinterpret_cast<const unsigned char*>(s++));
    while (*s) {
      result *= 256;
      result += static_cast<int>(
        *reinterpret_cast<const unsigned char*>(s++));
    }
    if (negative) return -result;
    return result;
  }

}}} // namespace iotbx::pdb::utils

#endif // IOTBX_PDB_UTILS_H
