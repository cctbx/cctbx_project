#ifndef CCTBX_CRYSTAL_PAIR_TABLES_H
#define CCTBX_CRYSTAL_PAIR_TABLES_H

#include <cctbx/sgtbx/rt_mx.h>
#include <map>

namespace cctbx { namespace crystal {

  typedef sgtbx::rt_mx pair_sym_op;
  typedef std::vector<sgtbx::rt_mx> pair_sym_ops;
  typedef std::map<unsigned, pair_sym_ops> pair_sym_dict;
  typedef af::shared<pair_sym_dict> pair_sym_table;

  typedef std::set<unsigned> pair_asu_j_sym_group;
  typedef std::vector<pair_asu_j_sym_group> pair_asu_j_sym_groups;
  typedef std::map<unsigned, pair_asu_j_sym_groups> pair_asu_dict;
  typedef af::shared<pair_asu_dict> pair_asu_table_table;

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_PAIR_TABLES_H
