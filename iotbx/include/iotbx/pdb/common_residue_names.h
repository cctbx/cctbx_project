#ifndef IOTBX_PDB_COMMON_RESIDUE_NAMES_H
#define IOTBX_PDB_COMMON_RESIDUE_NAMES_H

#include <iotbx/pdb/small_str.h>
#include <set>

namespace iotbx { namespace pdb { namespace common_residue_names {

  static const char* water[] = {
    "HOH ",
    "H2O ",
    "OH2 ",
    "DOD ",
    "D2O ",
    "OD2 ",
    "WAT ",
    "TIP ",
    "TIP3",
    0
  };

  inline
  const std::set<str4>&
  water_set()
  {
    static std::set<str4> result;
    if (result.size() == 0) {
      for(const char** n=water; *n; n++) result.insert(str4(*n));
    }
    return result;
  }

}}} // namespace iotbx::pdb::common_residue_names

#endif // IOTBX_PDB_COMMON_RESIDUE_NAMES_H
