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
  /*! Throws std::invalid_argument on anything flex.double rejects on every
      platform: empty cell, leading whitespace, trailing characters, a
      malformed esd bracket, overflow; and on what flex.double accepts only
      on some platforms: hex floats, underflow and subnormal values.
      Whitespace is the ASCII set, whatever the locale.
   */
  double
  parse_cell_double(std::string const& cell);

  //! Integer cell: optional ASCII whitespace around, optional sign, digits.
  bool
  try_parse_cell_int(std::string const& cell, long long& value);

  //! As try_parse_cell_int, throwing std::invalid_argument on failure.
  long long
  parse_cell_int(std::string const& cell);

  //! The hierarchy of an _atom_site loop, built in one pass over the rows.
  /*! A port of iotbx.pdb.mmcif._pdb_hierarchy_builder_python, which stays
      as the reference of tst_mmcif_hierarchy_parity.py; the hierarchies are
      identical except for the two differences below, which
      tst_mmcif_hierarchy_cpp.py asserts on both builders.

      Error reporting. A bad numeric, charge, id, seq id or anisotrop cell,
      a null chain id (., ? or a blank), an _atom_site_anisotrop column
      whose row count differs from that of
      _atom_site_anisotrop.id, or a text cell with a byte outside the CIF
      1.1 character set (printable ASCII, tab, line feed, carriage return)
      raises std::invalid_argument (ValueError in Python, a CifBuilderError
      through hierarchy_from_atom_site) naming the tag, the row (1-based
      within the loop) and the value, with bytes outside printable ASCII
      as \xHH. The text rule applies to every text column read,
      label_asym_id included. A missing type_symbol or numeric column is a
      CifBuilderError from atom_site_columns_from_block, before this
      function runs. The Python builder raised assorted exceptions, mostly
      without the location; it joined anisotrop columns of unequal length
      silently; a byte it could not decode as UTF-8 raised
      UnicodeDecodeError in any column it read, label_asym_id included; and
      it stored valid UTF-8 text (a chain id written as a Greek letter) in
      the std::string fields: chain id, resname, altloc, model id. Such
      text is not CIF 1.1, and only a non-standard file or a direct
      assignment puts it into a hierarchy; since the writer emits whatever
      a hierarchy holds, such a hierarchy is rejected when read back, by
      pickle and deepcopy too when it does not fit PDB format.

      Cell grammar. Every numeric cell (Cartn, B, occupancy, anisotrop,
      resolution, fp, fdp) follows parse_cell_double: an esd bracket and an
      exponent without its E (1.5+03, read as 1500) are accepted; blanks
      around the number, underscores, hex, overflow such as 1e500,
      underflow such as 1e-400 and subnormal values such as 4.9e-324 or
      1e-310 are rejected. Integer cells (id, seq id, charge) follow a
      strict ASCII integer grammar: optional sign, digits, ASCII blanks
      around. The Python builder read Cartn, B, occupancy and anisotrop
      cells with flex.double, whose grammar differs only by hex, underflow
      and subnormals, accepted where the platform's strtod accepts them;
      it read resolution, fp and fdp with float() and the integer cells
      with int(): blanks around a float, underscores, Unicode digits and
      Unicode blanks accepted, neither an esd bracket nor a bare exponent,
      inf on overflow, 0.0 on underflow and the subnormal itself (seq id
      1_0 read as 10).
   */
  hierarchy::root
  build_hierarchy(atom_site_columns const& columns);

}}} // namespace iotbx::pdb::mmcif

#endif // IOTBX_PDB_MMCIF_HIERARCHY_H
