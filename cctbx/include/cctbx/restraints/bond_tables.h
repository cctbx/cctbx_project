#ifndef CCTBX_RESTRAINTS_BOND_TABLES_H
#define CCTBX_RESTRAINTS_BOND_TABLES_H

#include <cctbx/restraints/bond.h>
#include <map>

namespace cctbx { namespace restraints {

  typedef std::map<unsigned, bond_params> bond_params_dict;
  typedef af::shared<bond_params_dict> bond_params_table;

  typedef sgtbx::rt_mx bond_sym_op;
  typedef std::vector<sgtbx::rt_mx> bond_sym_ops;
  typedef std::map<unsigned, bond_sym_ops> bond_sym_dict;
  typedef af::shared<bond_sym_dict> bond_sym_table;

  typedef std::set<unsigned> bond_asu_j_sym_group;
  typedef std::vector<bond_asu_j_sym_group> bond_asu_j_sym_groups;
  typedef std::map<unsigned, bond_asu_j_sym_groups> bond_asu_dict;
  typedef af::shared<bond_asu_dict> bond_asu_table;

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_BOND_TABLES_H
