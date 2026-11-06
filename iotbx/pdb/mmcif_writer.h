#ifndef IOTBX_PDB_MMCIF_WRITER_H
#define IOTBX_PDB_MMCIF_WRITER_H

#include <iotbx/pdb/hierarchy.h>
#include <scitbx/array_family/shared.h>
#include <string>

namespace iotbx { namespace pdb { namespace mmcif {

  namespace af = scitbx::af;

  //! The label values of the mmCIF writer that Python resolves.
  /*! One value per chain, residue group or atom group of the hierarchy,
      over all models in hierarchy order (root._atom_site_labels).
   */
  struct atom_site_labels
  {
    af::shared<std::string>
      chain_auth_asym_id, rg_label_asym_id, rg_auth_seq_id, ag_label_seq_id;
  };

  //! The _atom_site and _atom_site_anisotrop columns of a hierarchy.
  /*! The _atom_site members have one value per atom; auth_segid and
      auth_break are empty when not requested. The aniso members have one
      value per atom with U. chem_comp_ids are the residue names, sorted,
      each once; struct_asym_ids the label_asym_ids in order of first
      appearance.
   */
  struct atom_site_output
  {
    af::shared<std::string>
      group_pdb, id, label_atom_id, auth_atom_id, label_alt_id,
      label_comp_id, auth_asym_id, auth_seq_id, ins_code,
      cartn_x, cartn_y, cartn_z, occupancy, b_iso, type_symbol,
      formal_charge, fp, fdp, resolution, label_asym_id, label_entity_id,
      label_seq_id, model_num, auth_segid, auth_break;
    af::shared<std::string>
      aniso_id, aniso_auth_atom_id, aniso_label_alt_id, aniso_auth_comp_id,
      aniso_auth_asym_id, aniso_auth_seq_id, aniso_ins_code,
      aniso_u11, aniso_u22, aniso_u33, aniso_u12, aniso_u13, aniso_u23;
    af::shared<std::string> chem_comp_ids, struct_asym_ids;
  };

  //! The columns of root.as_cif_block, in one pass over the atoms.
  /*! A port of the atom loop of iotbx.pdb.hierarchy.root.
      _as_cif_block_python, which stays as the reference of
      iotbx/regression/tst_mmcif_writer_cpp.py: every value is the same
      string. Python resolves the per-chain and per-residue values (labels)
      and whether the hierarchy has break records (output_break_records is
      used as given).

      One difference: label_alt_id is the altloc of the atom group that
      holds the atom; the reference read the altloc of atom.parent(). The
      two differ when parent() is another atom group or none: atoms added
      with atom_group::append_atom_with_other_parent, the atoms of a
      select() without copy_atoms, whose parents are the original's atom
      groups, and the atoms such a selection shares with its original after
      remove_alt_confs() or remove_hd() on either. Where the reference
      raised AttributeError for a missing parent, this writes the value.

      Text is stripped as Python's str.strip() strips it: the ASCII blanks
      and the Unicode ones, as UTF-8. A double is written as Python's
      "%.<precision>f" % value writes it: snprintf, with "nan", "inf" and
      "-inf" for the values that are not finite, and '.' as the decimal
      point whatever the C locale. On an x86 build with -ffast-math, the
      test of fp and fdp against zero takes a NaN for zero, so an atom
      whose other value is 0 or NaN gets "." in both columns. No CI job
      that runs the iotbx tests, and no released build, uses that flag.
      Throws std::invalid_argument when a
      labels array does not have one value per chain, residue group or
      atom group, when a value is written with a negative precision, and
      when a serial is not a hybrid-36 number in a single-model hierarchy.
   */
  atom_site_output
  atom_site_columns_from_hierarchy(
    hierarchy::root const& root,
    atom_site_labels const& labels,
    int coordinate_precision,
    int occupancy_precision,
    int b_iso_precision,
    int u_aniso_precision,
    int resolution_precision,
    bool segid_as_auth_segid,
    bool output_break_records);

}}} // namespace iotbx::pdb::mmcif

#endif // IOTBX_PDB_MMCIF_WRITER_H
