#include <iotbx/pdb/mmcif_hierarchy.h>
#include <iotbx/pdb/hybrid_36_cpp.h>
#include <iotbx/pdb/common_residue_names.h>
#include <iotbx/error.h>
#include <cctbx/adptbx.h>
#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cctype>
#include <climits>
#include <locale.h>
// strtod_l / newlocale (POSIX) and _strtod_l / _create_locale (MSVC) are
// declared through <cstdlib> and <locale.h>; this is the include set
// xcif/numeric.cpp builds with on macOS, Linux and Windows.

namespace iotbx { namespace pdb { namespace mmcif {

  namespace {

    void
    throw_cell_error(
      const char* tag,
      std::size_t row,
      const char* reason,
      std::string const& value)
    {
      std::ostringstream o;
      o << tag << " row " << (row + 1) << ": " << reason << ": '" << value << "'";
      throw std::invalid_argument(o.str());
    }

    void
    check_length(
      const char* name,
      af::shared<std::string> const& column,
      std::size_t n,
      bool optional)
    {
      if (optional && column.size() == 0) return;
      if (column.size() != n) {
        std::ostringstream o;
        o << "atom_site_columns." << name << " has " << column.size()
          << " rows, expected " << n;
        throw std::invalid_argument(o.str());
      }
    }

    void
    check_lengths(atom_site_columns const& c, std::size_t n)
    {
      check_length("type_symbol", c.type_symbol, n, false);
      check_length("atom_name", c.atom_name, n, false);
      check_length("alt_id", c.alt_id, n, false);
      check_length("label_asym_id", c.label_asym_id, n, false);
      check_length("auth_asym_id", c.auth_asym_id, n, false);
      check_length("comp_id", c.comp_id, n, false);
      check_length("seq_id", c.seq_id, n, false);
      check_length("cartn_x", c.cartn_x, n, false);
      check_length("cartn_y", c.cartn_y, n, false);
      check_length("cartn_z", c.cartn_z, n, false);
      check_length("b_iso", c.b_iso, n, false);
      check_length("occupancy", c.occupancy, n, false);
      check_length("model_num", c.model_num, n, true);
      check_length("atom_id", c.atom_id, n, true);
      check_length("ins_code", c.ins_code, n, true);
      check_length("group_pdb", c.group_pdb, n, true);
      check_length("segid", c.segid, n, true);
      check_length("auth_break", c.auth_break, n, true);
      check_length("formal_charge", c.formal_charge, n, true);
      check_length("resolution", c.resolution, n, true);
      check_length("fp", c.fp, n, true);
      check_length("fdp", c.fdp, n, true);
      std::size_t m = c.aniso_id.size();
      check_length("aniso_11", c.aniso_11, m, false);
      check_length("aniso_22", c.aniso_22, m, false);
      check_length("aniso_33", c.aniso_33, m, false);
      check_length("aniso_12", c.aniso_12, m, false);
      check_length("aniso_13", c.aniso_13, m, false);
      check_length("aniso_23", c.aniso_23, m, false);
    }

    // strtod is locale-sensitive; CIF always uses '.'. Same approach as
    // xcif/numeric.cpp: a C-locale strtod.
#ifdef _WIN32
    double strtod_c(const char* s, char** end)
    {
      static _locale_t loc = _create_locale(LC_ALL, "C");
      return _strtod_l(s, end, loc);
    }
#else
    double strtod_c(const char* s, char** end)
    {
      static locale_t loc = ::newlocale(LC_ALL_MASK, "C", (locale_t)0);
      return ::strtod_l(s, end, loc);
    }
#endif

    bool
    all_digits_with_optional_sign(std::string const& s)
    {
      std::size_t i = 0;
      if (i < s.size() && (s[i] == '+' || s[i] == '-')) i++;
      if (i == s.size()) return false;
      for(; i < s.size(); i++) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
      }
      return true;
    }

    // flex.double validates the esd with boost::lexical_cast<int>: an
    // optional sign, digits, and a value that fits in an int.
    bool
    is_int_esd(std::string const& s)
    {
      if (!all_digits_with_optional_sign(s)) return false;
      errno = 0;
      char* end = 0;
      long v = std::strtol(s.c_str(), &end, 10);
      if (*end != '\0' || errno == ERANGE) return false;
      return v >= INT_MIN && v <= INT_MAX;
    }

  } // namespace <anonymous>

  double
  parse_cell_double(std::string const& s)
  {
    // The grammar of scitbx flex_double.cpp from_std_string.
    std::string::size_type ob = s.find_first_of('(');
    std::string::size_type cb = s.find_last_of(')');
    std::string value;
    std::string esd;
    if (ob == std::string::npos) {
      value = s;
      if (value.size() == 0) {
        throw std::invalid_argument("Empty string (floating-point value expected).");
      }
    }
    else if (cb == std::string::npos) {
      throw std::invalid_argument("Missing closing parenthesis: \"" + s + "\"");
    }
    else {
      if (cb != s.size() - 1) {
        throw std::invalid_argument(
          "Unexpected trailing characters after \")\": \"" + s + "\"");
      }
      if (ob == 0) {
        throw std::invalid_argument("Empty value part: \"" + s + "\"");
      }
      if (cb - ob - 1 == 0) {
        throw std::invalid_argument("Empty esd part: \"" + s + "\"");
      }
      value = s.substr(0, ob);
      esd = s.substr(ob + 1, cb - ob - 1);
    }
    { // insert E if necessary (e.g. .34+05 -> .34E+05)
      std::string::size_type ip = value.find_last_of('+');
      std::string::size_type in = value.find_last_of('-');
      if (ip == std::string::npos) ip = 0;
      if (in == std::string::npos) in = 0;
      std::string::size_type is = std::max(ip, in);
      if (is != 0) {
        char c = value[is - 1];
        if (c == '.' || std::isdigit(static_cast<unsigned char>(c))) {
          value = value.substr(0, is) + "E" + value.substr(is);
        }
      }
    }
    if (std::isspace(static_cast<unsigned char>(value[0]))) {
      throw std::invalid_argument("Invalid floating-point value: \"" + s + "\"");
    }
    errno = 0;
    char* end = 0;
    double result = strtod_c(value.c_str(), &end);
    if (   end == value.c_str() || end != value.c_str() + value.size()
        || errno == ERANGE) {
      throw std::invalid_argument("Invalid floating-point value: \"" + s + "\"");
    }
    if (esd.size() != 0 && !is_int_esd(esd)) {
      throw std::invalid_argument("Invalid esd part: \"" + s + "\"");
    }
    return result;
  }

  bool
  try_parse_cell_int(std::string const& cell, long& value)
  {
    std::size_t b = 0;
    std::size_t e = cell.size();
    while (b < e && std::isspace(static_cast<unsigned char>(cell[b]))) b++;
    while (e > b && std::isspace(static_cast<unsigned char>(cell[e - 1]))) e--;
    if (b == e) return false;
    std::string core = cell.substr(b, e - b);
    if (!all_digits_with_optional_sign(core)) return false;
    errno = 0;
    char* end = 0;
    long v = std::strtol(core.c_str(), &end, 10);
    if (*end != '\0' || errno == ERANGE) return false;
    value = v;
    return true;
  }

  long
  parse_cell_int(std::string const& cell)
  {
    long value;
    if (!try_parse_cell_int(cell, value)) {
      throw std::invalid_argument("not an integer: '" + cell + "'");
    }
    return value;
  }

  namespace {

    std::string
    strip_spaces(std::string const& s)
    {
      std::size_t b = s.find_first_not_of(" \t\n\r\f\v");
      if (b == std::string::npos) return std::string();
      std::size_t e = s.find_last_not_of(" \t\n\r\f\v");
      return s.substr(b, e - b + 1);
    }

    bool
    is_polymer_residue(std::string const& raw_comp_id)
    {
      // Python: common_residue_names_get_class(resname.strip()) in the four
      // classes; strip() also removes the newlines of a text-field value
      std::string name = strip_spaces(raw_comp_id);
      if (name.size() == 0) return false;
      std::string const& cls = common_residue_names::get_class(name);
      return (   cls == "common_amino_acid" || cls == "modified_amino_acid"
              || cls == "common_rna_dna" || cls == "modified_rna_dna");
    }

    //! Port of iotbx.pdb.mmcif.format_pdb_atom_name, plus F8 for the empty name.
    std::string
    format_atom_name(std::string const& name, std::string const& element)
    {
      if (name.size() >= 4) return name;
      std::string s = strip_spaces(name);
      if (s.size() == 0) return std::string("    ");
      if (element.size() == 1) {
        if (std::toupper(static_cast<unsigned char>(s[0]))
            == std::toupper(static_cast<unsigned char>(element[0]))) {
          s = " " + s;
        }
      }
      while (s.size() < 4) s += ' ';
      return s;
    }

    std::string
    pad_resname(std::string const& comp_id)
    {
      // Python "%3s" % resname: right-justify to 3, longer unchanged
      if (comp_id.size() >= 3) return comp_id;
      return std::string(3 - comp_id.size(), ' ') + comp_id;
    }

    std::string
    resseq_for(std::string const& seq_id, std::size_t row)
    {
      static const char* tag = "_atom_site.auth_seq_id (or label_seq_id)";
      long v;
      if (try_parse_cell_int(seq_id, v)) {
        if (v >= -999 && v <= 2436111) {
          return hybrid_36::encode(4U, static_cast<int>(v));
        }
      }
      // F7: a blank seq id is rejected like the nulls; kept as four spaces
      // it would merge neighbouring residues into one residue group
      if (seq_id == "." || seq_id == "?" || strip_spaces(seq_id).size() == 0) {
        throw_cell_error(tag, row, "residue sequence id is a null value", seq_id);
      }
      if (seq_id.size() == 4) return seq_id;
      if (seq_id.size() < 4) {
        return std::string(4 - seq_id.size(), ' ') + seq_id;   // F7
      }
      throw_cell_error(tag, row,
        "not an integer and longer than four characters", seq_id);
      return std::string(); // not reached
    }

    std::string
    serial_for(atom_site_columns const& c, std::size_t row)
    {
      if (c.atom_id.size() == 0) {
        return hybrid_36::encode(5U, static_cast<int>(row + 1));   // F4
      }
      std::string const& cell = c.atom_id[row];
      long v;
      if (!try_parse_cell_int(cell, v) || v < -9999 || v > 87440031) {
        throw_cell_error("_atom_site.id", row, "not an integer in hybrid-36 range", cell);
      }
      return hybrid_36::encode(5U, static_cast<int>(v));
    }

    double
    parse_or_throw(const char* tag, af::shared<std::string> const& col, std::size_t row)
    {
      try {
        return parse_cell_double(col[row]);
      }
      catch (std::invalid_argument const& e) {
        throw_cell_error(tag, row, e.what(), col[row]);
      }
      return 0; // not reached
    }

    double
    parse_b_or_occ(const char* tag, af::shared<std::string> const& col, std::size_t row)
    {
      std::string const& cell = col[row];
      // integrative models write . or ? for occupancy and B
      if (cell == "." || cell == "?") return 0.0;
      return parse_or_throw(tag, col, row);
    }

    std::string
    charge_for(std::string const& cell, std::size_t row)
    {
      // Python: sign from a leading or trailing '-', strip " -+", int, 0 -> no sign
      std::string sign = "+";
      if (cell.size() != 0 && (cell[0] == '-' || cell[cell.size() - 1] == '-')) sign = "-";
      std::size_t b = cell.find_first_not_of(" -+");
      std::string digits;
      if (b != std::string::npos) {
        std::size_t e = cell.find_last_not_of(" -+");
        digits = cell.substr(b, e - b + 1);
      }
      long v;
      if (!try_parse_cell_int(digits, v)) {
        throw_cell_error("_atom_site.pdbx_formal_charge", row, "not an integer charge", cell);
      }
      if (v == 0) sign = "";
      std::ostringstream o;
      o << v << sign;
      std::string result = o.str();
      if (result.size() > 2) {
        throw_cell_error("_atom_site.pdbx_formal_charge", row,
          "charge longer than two characters", cell);
      }
      return result;
    }

    typedef std::map<std::pair<std::string, std::string>, hierarchy::atom_group>
      atom_group_map;

  } // namespace <anonymous>

  hierarchy::root
  build_hierarchy(atom_site_columns const& c)
  {
    const std::size_t n = c.atom_name.size();
    check_lengths(c, n);
    const bool have_model = c.model_num.size() != 0;
    const bool have_ins = c.ins_code.size() != 0;
    const bool have_group = c.group_pdb.size() != 0;
    const bool have_segid = c.segid.size() != 0;
    const bool have_break = c.auth_break.size() != 0;
    const bool have_charge = c.formal_charge.size() != 0;
    const bool have_res = c.resolution.size() != 0;
    const bool have_fp = c.fp.size() != 0;
    const bool have_fdp = c.fdp.size() != 0;
    const bool have_aniso = c.aniso_id.size() != 0;
    if (have_aniso && n != 0 && c.atom_id.size() == 0) {
      throw std::invalid_argument(
        "_atom_site_anisotrop present but _atom_site.id missing");
    }
    std::map<std::string, std::size_t> aniso_index;
    std::vector<scitbx::sym_mat3<double> > aniso_uij;
    if (have_aniso) {
      const std::size_t m = c.aniso_id.size();
      aniso_uij.reserve(m);
      af::shared<std::string> const* cols[6] = {
        &c.aniso_11, &c.aniso_22, &c.aniso_33, &c.aniso_12, &c.aniso_13, &c.aniso_23};
      static const char* tags_u[6] = {
        "_atom_site_anisotrop.U[1][1]", "_atom_site_anisotrop.U[2][2]",
        "_atom_site_anisotrop.U[3][3]", "_atom_site_anisotrop.U[1][2]",
        "_atom_site_anisotrop.U[1][3]", "_atom_site_anisotrop.U[2][3]"};
      static const char* tags_b[6] = {
        "_atom_site_anisotrop.B[1][1]", "_atom_site_anisotrop.B[2][2]",
        "_atom_site_anisotrop.B[3][3]", "_atom_site_anisotrop.B[1][2]",
        "_atom_site_anisotrop.B[1][3]", "_atom_site_anisotrop.B[2][3]"};
      for(std::size_t r=0;r<m;r++) {
        scitbx::sym_mat3<double> u;
        for(unsigned k=0;k<6;k++) {
          u[k] = parse_or_throw(c.aniso_is_b ? tags_b[k] : tags_u[k], *cols[k], r);
        }
        if (c.aniso_is_b) u = cctbx::adptbx::b_as_u(u);
        aniso_uij.push_back(u);
        aniso_index[c.aniso_id[r]] = r;   // last occurrence wins
      }
    }
    hierarchy::root result;
    hierarchy::model model;
    hierarchy::chain chain;
    hierarchy::residue_group rg;
    atom_group_map groups;
    bool rg_is_first_in_chain = true;
    std::string prev_model, prev_auth, prev_label, prev_seq, prev_icode;
    for(std::size_t i=0;i<n;i++) {
      std::string const model_num = have_model ? c.model_num[i] : std::string();
      std::string const& auth = c.auth_asym_id[i];
      std::string const& label = c.label_asym_id[i];
      std::string icode(" ");
      if (have_ins) {
        icode = c.ins_code[i];
        if (icode == "?" || icode == ".") icode = " ";
        if (icode.size() > 1) {
          throw_cell_error("_atom_site.pdbx_PDB_ins_code", i,
            "insertion code longer than one character", c.ins_code[i]);
        }
      }
      bool new_model = (i == 0) || model_num != prev_model;
      bool new_chain = new_model || auth != prev_auth
        || (label != prev_label && i > 0 && is_polymer_residue(c.comp_id[i-1]));
      bool new_rg = new_chain || c.seq_id[i] != prev_seq || icode != prev_icode;   // F9: new_chain
      if (new_model) {
        model = hierarchy::model(model_num.c_str());
        result.append_model(model);
      }
      if (new_chain) {
        chain = hierarchy::chain(auth.c_str());
        model.append_chain(chain);
      }
      if (new_rg) {
        rg = hierarchy::residue_group(resseq_for(c.seq_id[i], i).c_str(), icode.c_str(), true);
        chain.append_residue_group(rg);
        rg_is_first_in_chain = new_chain;
        groups.clear();
      }
      std::string altloc = c.alt_id[i];
      if (altloc == "." || altloc == "?") altloc.clear();
      std::pair<std::string, std::string> key(altloc, c.comp_id[i]);
      atom_group_map::iterator g = groups.find(key);
      if (g == groups.end()) {
        hierarchy::atom_group ag(altloc, pad_resname(c.comp_id[i]));
        if (altloc.size() == 0) rg.insert_atom_group(0, ag);
        else rg.append_atom_group(ag);
        g = groups.insert(std::make_pair(key, ag)).first;
      }
      hierarchy::atom atom;
      g->second.append_atom(atom);
      if (c.type_symbol[i].size() > 2) {
        throw_cell_error("_atom_site.type_symbol", i,
          "element longer than two characters", c.type_symbol[i]);
      }
      atom.set_element(c.type_symbol[i].c_str());
      std::string const name = format_atom_name(c.atom_name[i], c.type_symbol[i]);
      if (name.size() > 4) {
        throw_cell_error("_atom_site.label_atom_id (or auth_atom_id)", i,
          "atom name longer than four characters", c.atom_name[i]);
      }
      atom.set_name(name.c_str());
      atom.set_xyz(scitbx::vec3<double>(
        parse_or_throw("_atom_site.Cartn_x", c.cartn_x, i),
        parse_or_throw("_atom_site.Cartn_y", c.cartn_y, i),
        parse_or_throw("_atom_site.Cartn_z", c.cartn_z, i)));
      atom.set_b(parse_b_or_occ("_atom_site.B_iso_or_equiv", c.b_iso, i));
      atom.set_occ(parse_b_or_occ("_atom_site.occupancy", c.occupancy, i));
      atom.set_serial(serial_for(c, i).c_str());
      if (have_segid) {
        std::string s = c.segid[i].substr(0, 4);
        while (s.size() < 4) s += ' ';
        atom.set_segid(s.c_str());
      }
      else {
        atom.set_segid("    ");
      }
      if (have_break && !rg_is_first_in_chain && c.auth_break[i] == "1") {
        rg.data->link_to_previous = false;
      }
      if (have_group && c.group_pdb[i] == "HETATM") atom.set_hetero(true);
      if (have_charge) {
        std::string const& cell = c.formal_charge[i];
        if (cell != "?" && cell != ".") atom.set_charge(charge_for(cell, i).c_str());
      }
      if (have_res) {
        std::string const& cell = c.resolution[i];
        if (cell != "?" && cell != ".") {
          atom.set_resolution(parse_or_throw("_atom_site.phenix_resolution", c.resolution, i));
        }
      }
      if (have_fp) {
        std::string const& cell = c.fp[i];
        if (cell != "?" && cell != ".") {
          atom.set_fp(parse_or_throw("_atom_site.phenix_scat_dispersion_real", c.fp, i));
        }
      }
      if (have_fdp) {
        std::string const& cell = c.fdp[i];
        if (cell != "?" && cell != ".") {
          atom.set_fdp(parse_or_throw("_atom_site.phenix_scat_dispersion_imag", c.fdp, i));
        }
      }
      if (have_aniso) {
        std::map<std::string, std::size_t>::const_iterator a = aniso_index.find(c.atom_id[i]);
        if (a != aniso_index.end()) atom.set_uij(aniso_uij[a->second]);
      }
      prev_model = model_num;
      prev_auth = auth;
      prev_label = label;
      prev_seq = c.seq_id[i];
      prev_icode = icode;
    }
    if (result.models_size() == 1) {
      result.models()[0].data->id = "";
    }
    return result;
  }

}}} // namespace iotbx::pdb::mmcif
