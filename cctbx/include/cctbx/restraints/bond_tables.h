#ifndef CCTBX_RESTRAINTS_BOND_TABLES_H
#define CCTBX_RESTRAINTS_BOND_TABLES_H

#include <cctbx/restraints/bond.h>
#include <map>

namespace cctbx { namespace restraints {

  typedef sgtbx::rt_mx bond_sym_op;
  typedef std::vector<sgtbx::rt_mx> bond_sym_ops;
  typedef std::map<unsigned, bond_sym_ops> bond_sym_dict;
  typedef af::shared<bond_sym_dict> bond_sym_table;

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_BOND_TABLES_H
