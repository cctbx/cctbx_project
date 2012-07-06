#ifndef IOTBX_PDB_HIERARCHY_ATOMS_H
#define IOTBX_PDB_HIERARCHY_ATOMS_H

#include <iotbx/pdb/hierarchy.h>

namespace iotbx { namespace pdb { namespace hierarchy { namespace atoms {

  af::shared<std::string>
  extract_serial(
    af::const_ref<atom> const& atoms);

  af::shared<std::string>
  extract_name(
    af::const_ref<atom> const& atoms);

  af::shared<std::string>
  extract_segid(
    af::const_ref<atom> const& atoms);

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

#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
  af::shared<sym_mat3>
  extract_siguij(
    af::const_ref<atom> const& atoms);
#endif

  af::shared<std::size_t>
  extract_hetero(
    af::const_ref<atom> const& atoms);

  af::shared<std::string>
  extract_element(
    af::const_ref<atom> const& atoms,
    bool strip=false);

  af::shared<std::size_t>
  extract_i_seq(
    af::const_ref<atom> const& atoms);

  af::shared<std::size_t>
  extract_tmp_as_size_t(
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

#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
  void
  set_siguij(
    af::ref<atom> const& atoms,
    af::const_ref<sym_mat3> const& new_siguij);
#endif

  void
  reset_serial(
    af::const_ref<atom> const& atoms,
    int first_value=1);

  void
  reset_i_seq(
    af::const_ref<atom> const& atoms);

  std::size_t
  set_chemical_element_simple_if_necessary(
    af::ref<atom> const& atoms,
    bool tidy_existing=true);

  class atom_tmp_sentinel : boost::noncopyable
  {
    protected:
      std::vector<atom> atoms_;

    public:
      atom_tmp_sentinel(
        af::const_ref<atom> const& atoms);

      ~atom_tmp_sentinel();
  };

  std::auto_ptr<atom_tmp_sentinel>
  reset_tmp(
    af::const_ref<atom> const& atoms,
    int first_value=0,
    int increment=1);

  std::auto_ptr<atom_tmp_sentinel>
  reset_tmp_for_occupancy_groups_simple(
    af::const_ref<atom> const& atoms);

}}}} // namespace iotbx::pdb::hierarchy::atoms

#endif // IOTBX_PDB_HIERARCHY_ATOMS_H
