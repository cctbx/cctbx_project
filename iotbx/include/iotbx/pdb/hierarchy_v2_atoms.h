#ifndef IOTBX_PDB_HIERARCHY_V2_ATOMS_H
#define IOTBX_PDB_HIERARCHY_V2_ATOMS_H

#include <iotbx/pdb/hierarchy_v2.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 { namespace atoms {

  af::shared<vec3>
  extract_xyz(
    af::const_ref<atom> const& atoms);

  af::shared<vec3>
  extract_sigxyz(
    af::const_ref<atom> const& atoms);

  af::shared<double>
  extract_occ(
    af::const_ref<atom> const& atoms);

  af::shared<double>
  extract_sigocc(
    af::const_ref<atom> const& atoms);

  af::shared<double>
  extract_b(
    af::const_ref<atom> const& atoms);

  af::shared<double>
  extract_sigb(
    af::const_ref<atom> const& atoms);

  af::shared<sym_mat3>
  extract_uij(
    af::const_ref<atom> const& atoms);

  af::shared<sym_mat3>
  extract_siguij(
    af::const_ref<atom> const& atoms);

  af::shared<std::size_t>
  extract_hetero(
    af::const_ref<atom> const& atoms);

  void
  set_xyz(
    af::ref<atom> const& atoms,
    af::const_ref<vec3> const& new_xyz);

  void
  set_sigxyz(
    af::ref<atom> const& atoms,
    af::const_ref<vec3> const& new_sigxyz);

  void
  set_occ(
    af::ref<atom> const& atoms,
    af::const_ref<double> const& new_occ);

  void
  set_sigocc(
    af::ref<atom> const& atoms,
    af::const_ref<double> const& new_sigocc);

  void
  set_b(
    af::ref<atom> const& atoms,
    af::const_ref<double> const& new_b);

  void
  set_sigb(
    af::ref<atom> const& atoms,
    af::const_ref<double> const& new_sigb);

  void
  set_uij(
    af::ref<atom> const& atoms,
    af::const_ref<sym_mat3> const& new_uij);

  void
  set_siguij(
    af::ref<atom> const& atoms,
    af::const_ref<sym_mat3> const& new_siguij);

  void
  reset_tmp(
    af::const_ref<atom> const& atoms,
    int first_value=0,
    int increment=1);

  void
  reset_tmp_for_occupancy_groups_simple(
    af::const_ref<atom> const& atoms);

}}}} // namespace iotbx::pdb::hierarchy_v2::atoms

#endif // IOTBX_PDB_HIERARCHY_V2_ATOMS_H
