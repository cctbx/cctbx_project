#ifndef IOTBX_PDB_UTILS_H
#define IOTBX_PDB_UTILS_H

#include <boost/cstdint.hpp>

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
  boost::int64_t
  base_256_ordinal(const char* s);

}}} // namespace iotbx::pdb::utils

#endif // IOTBX_PDB_UTILS_H
