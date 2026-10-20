#ifndef IOTBX_PDB_MMCIF_HIERARCHY_H
#define IOTBX_PDB_MMCIF_HIERARCHY_H

#include <iotbx/pdb/hierarchy.h>
#include <scitbx/array_family/shared.h>
#include <string>

namespace iotbx { namespace pdb { namespace mmcif {

  namespace af = scitbx::af;

  //! The resolved columns of an _atom_site loop (and _atom_site_anisotrop).
  /*! Python resolves the tags and fallbacks and fills this struct with
      flex.std_string handles. Mandatory members have the row count n;
      optional members are empty when the column is absent. The aniso
      members have the row count of _atom_site_anisotrop or are empty.
   */
  struct atom_site_columns
  {
    af::shared<std::string>
      type_symbol, atom_name, alt_id, label_asym_id, auth_asym_id,
      comp_id, seq_id, cartn_x, cartn_y, cartn_z, b_iso, occupancy;
    af::shared<std::string>
      model_num, atom_id, ins_code, group_pdb, segid, auth_break,
      formal_charge, resolution, fp, fdp;
    af::shared<std::string>
      aniso_id, aniso_11, aniso_22, aniso_33, aniso_12, aniso_13, aniso_23;
    bool aniso_is_b;

    atom_site_columns() : aniso_is_b(false) {}
  };

  //! Parse a numeric cell with the grammar of flex.double(flex.std_string).
  /*! Throws std::invalid_argument on anything flex.double rejects: empty
      cell, leading whitespace, trailing characters, over- or underflow,
      malformed esd bracket.
   */
  double
  parse_cell_double(std::string const& cell);

  //! Integer cell: optional surrounding whitespace, optional sign, digits.
  bool
  try_parse_cell_int(std::string const& cell, long& value);

  //! As try_parse_cell_int, throwing std::invalid_argument on failure.
  long
  parse_cell_int(std::string const& cell);

  //! The hierarchy of an _atom_site loop, built in one pass over the rows;
  //! differences from the Python builder: F1-F11 below.
  /*! A bad cell raises std::invalid_argument (ValueError in Python) naming
      the tag, the row (1-based within the loop) and the value.

      The Python builder is iotbx.pdb.mmcif.pdb_hierarchy_builder_python;
      tst_mmcif_hierarchy_cpp.py exercises each difference. F5 for missing
      columns and F10 come from the column resolution in
      iotbx.pdb.mmcif.atom_site_columns_from_block.

      F1  anisotrop join on the exact _atom_site.id string (Python used the
          stripped hybrid-36 serial: lost the ANISOU of ids above 99999 or
          with leading zeros, joined atom 01 to anisotrop 1)
      F2  B[i][j] converted to U (Python: NameError)
      F3  no pdbx_PDB_model_num: one model, id "" (Python: TypeError)
      F4  no _atom_site.id: serial from the row order, 1-based; ValueError
          when there are atom rows and anisotrop rows (Python: TypeError)
      F5  bad numeric, charge, id, seq id or anisotrop cell, or missing
          type_symbol or numeric column: ValueError with tag, row and value
          (Python: assorted exceptions without the location)
      F6  no pdbx_PDB_ins_code: icode " ", as for "?" (Python: "")
      F7  non-integer seq id: right-justified into four if shorter,
          ValueError if longer; ".", "?" and blank: ValueError (Python:
          AssertionError unless the id had exactly four characters; a blank
          four-character id was accepted as resseq "    ")
      F8  atom name empty after stripping: four spaces (Python: IndexError)
      F9  every chain boundary starts a residue group (Python put the atoms
          of a label_asym_id chain break with unchanged seq id and icode
          into the previous chain and left the new chain empty)
      F10 single-item _atom_site: every item read as a one-element column
          (Python used the first character of group_PDB, id, model number,
          insertion code and phenix_scat_dispersion_real/_imag)
      F11 cell grammar: resolution, fp and fdp follow the flex.double
          grammar (esd bracket accepted, overflow such as 1e500 a
          ValueError); integer cells (id, seq id, charge) a strict integer
          grammar (optional sign and surrounding blanks, digits only).
          Python used float() and int(): no esd bracket, inf on overflow,
          underscores accepted (seq id 1_0 gave resseq "  10")
   */
  hierarchy::root
  build_hierarchy(atom_site_columns const& columns);

}}} // namespace iotbx::pdb::mmcif

#endif // IOTBX_PDB_MMCIF_HIERARCHY_H
