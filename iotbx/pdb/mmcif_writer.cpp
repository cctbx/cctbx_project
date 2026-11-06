#include <iotbx/pdb/mmcif_writer.h>
#include <clocale>
#include <cstdio>
#include <cstring>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace iotbx { namespace pdb { namespace mmcif {

  namespace {

    // The ASCII characters of Python's str.isspace().
    bool
    is_ascii_space(char c)
    {
      return c == ' ' || (c >= '\t' && c <= '\r')
        || (c >= '\x1c' && c <= '\x1f');
    }

    // The other characters of Python's str.isspace(), as UTF-8: U+0085,
    // U+00A0 (two bytes); U+1680, U+2000..U+200A, U+2028, U+2029, U+202F,
    // U+205F, U+3000 (three bytes).
    bool
    is_utf8_space(unsigned char a, unsigned char b)
    {
      return a == 0xc2 && (b == 0x85 || b == 0xa0);
    }

    bool
    is_utf8_space(unsigned char a, unsigned char b, unsigned char c)
    {
      if (a == 0xe1) return b == 0x9a && c == 0x80;
      if (a == 0xe2 && b == 0x80) {
        return (c >= 0x80 && c <= 0x8a) || c == 0xa8 || c == 0xa9
          || c == 0xaf;
      }
      if (a == 0xe2) return b == 0x81 && c == 0x9f;
      return a == 0xe3 && b == 0x80 && c == 0x80;
    }

    // The number of bytes of the str.isspace() character at the start of
    // [begin, end), or 0.
    std::size_t
    space_at_start(const char* begin, const char* end)
    {
      const unsigned char* s = reinterpret_cast<const unsigned char*>(begin);
      std::size_t n = end - begin;
      if (n >= 1 && is_ascii_space(begin[0])) return 1;
      if (n >= 2 && is_utf8_space(s[0], s[1])) return 2;
      if (n >= 3 && is_utf8_space(s[0], s[1], s[2])) return 3;
      return 0;
    }

    // The same for the character that ends at end.
    std::size_t
    space_at_end(const char* begin, const char* end)
    {
      const unsigned char* s = reinterpret_cast<const unsigned char*>(end);
      std::size_t n = end - begin;
      if (n >= 1 && is_ascii_space(end[-1])) return 1;
      if (n >= 2 && is_utf8_space(s[-2], s[-1])) return 2;
      if (n >= 3 && is_utf8_space(s[-3], s[-2], s[-1])) return 3;
      return 0;
    }

    // Python's str.strip() of the UTF-8 text [begin, end).
    std::string
    python_strip(const char* begin, const char* end)
    {
      while (true) {
        std::size_t n = space_at_start(begin, end);
        if (n == 0) break;
        begin += n;
      }
      while (true) {
        std::size_t n = space_at_end(begin, end);
        if (n == 0) break;
        end -= n;
      }
      return std::string(begin, end);
    }

    std::string
    python_strip(std::string const& s)
    {
      return python_strip(s.data(), s.data() + s.size());
    }

    std::string
    python_strip(const char* s)
    {
      return python_strip(s, s + std::strlen(s));
    }

    // Python's "%.<precision>f" % value. snprintf rounds as Python does,
    // but spells a NaN with its sign on some platforms (-nan, -nan(ind))
    // and writes the decimal point of LC_NUMERIC, a comma in many locales
    // once something calls setlocale(LC_ALL, "") (wxtbx.plots does).
    // radix is that decimal point.
    std::string
    format_fixed(double value, int precision, std::string const& radix)
    {
      if (precision < 0) {
        std::ostringstream o;
        o << "negative precision: " << precision;
        throw std::invalid_argument(o.str());
      }
      char buffer[64];
      int n = std::snprintf(buffer, sizeof(buffer), "%.*f", precision, value);
      if (n < 0) throw std::runtime_error("snprintf failed");
      std::string result;
      if (static_cast<std::size_t>(n) < sizeof(buffer)) {
        result.assign(buffer, n);
      }
      else {
        std::vector<char> big(n + 1);
        std::snprintf(&big[0], big.size(), "%.*f", precision, value);
        result.assign(&big[0], n);
      }
      // A value that is not finite has no digit first: the C library
      // writes nan, -nan, -nan(ind) or nan(snan), inf or -inf; Python
      // writes nan, inf or -inf. The text is tested, not std::isnan, which
      // -ffast-math folds to false.
      std::size_t first = result[0] == '-' ? 1 : 0;
      if (result[first] < '0' || result[first] > '9') {
        if (result[first] == 'n' || result[first] == 'N') return "nan";
        return first ? "-inf" : "inf";
      }
      if (radix != ".") {
        std::size_t i = result.find(radix);
        if (i != std::string::npos) result.replace(i, radix.size(), ".");
      }
      return result;
    }

    std::string
    decimal(long value)
    {
      char buffer[32];
      int n = std::snprintf(buffer, sizeof(buffer), "%ld", value);
      return std::string(buffer, n);
    }

    // _atom_site.pdbx_formal_charge as the Python reference writes it:
    // charge_tidy, then the sign in front ("2-" is -2, "2+" is 2), ? for
    // none.
    std::string
    formal_charge(hierarchy::atom const& atom)
    {
      boost::optional<std::string> tidy = atom.charge_tidy();
      if (!tidy) return "?";
      std::string charge = python_strip(*tidy);
      if (charge.empty()) return "?";
      if (charge.size() == 2) {
        if (charge[1] == '-') return std::string("-") + charge[0];
        return charge.substr(0, 1);
      }
      return charge;
    }

    void
    check_size(
      const char* name,
      std::size_t size,
      const char* level,
      std::size_t expected)
    {
      if (size != expected) {
        std::ostringstream o;
        o << "atom_site_labels." << name << " has " << size
          << " values, expected one per " << level << ": " << expected;
        throw std::invalid_argument(o.str());
      }
    }

  } // namespace <anonymous>

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
    bool output_break_records)
  {
    std::vector<hierarchy::model> const& models = root.models();
    std::size_t n_chains = 0;
    std::size_t n_rgs = 0;
    std::size_t n_ags = 0;
    std::size_t n_atoms = 0;
    for (std::size_t i_model = 0; i_model < models.size(); i_model++) {
      std::vector<hierarchy::chain> const& chains = models[i_model].chains();
      n_chains += chains.size();
      for (std::size_t i_chain = 0; i_chain < chains.size(); i_chain++) {
        std::vector<hierarchy::residue_group> const&
          rgs = chains[i_chain].residue_groups();
        n_rgs += rgs.size();
        for (std::size_t i_rg = 0; i_rg < rgs.size(); i_rg++) {
          std::vector<hierarchy::atom_group> const&
            ags = rgs[i_rg].atom_groups();
          n_ags += ags.size();
          for (std::size_t i_ag = 0; i_ag < ags.size(); i_ag++) {
            n_atoms += ags[i_ag].atoms_size();
          }
        }
      }
    }
    check_size("chain_auth_asym_id", labels.chain_auth_asym_id.size(),
      "chain", n_chains);
    check_size("rg_label_asym_id", labels.rg_label_asym_id.size(),
      "residue group", n_rgs);
    check_size("rg_auth_seq_id", labels.rg_auth_seq_id.size(),
      "residue group", n_rgs);
    check_size("ag_label_seq_id", labels.ag_label_seq_id.size(),
      "atom group", n_ags);

    atom_site_output out;
    af::shared<std::string>* per_atom[] = {
      &out.group_pdb, &out.id, &out.label_atom_id, &out.auth_atom_id,
      &out.label_alt_id, &out.label_comp_id, &out.auth_asym_id,
      &out.auth_seq_id, &out.ins_code, &out.cartn_x, &out.cartn_y,
      &out.cartn_z, &out.occupancy, &out.b_iso, &out.type_symbol,
      &out.formal_charge, &out.fp, &out.fdp, &out.resolution,
      &out.label_asym_id, &out.label_entity_id, &out.label_seq_id,
      &out.model_num};
    for (std::size_t i = 0; i < sizeof(per_atom) / sizeof(per_atom[0]); i++) {
      per_atom[i]->reserve(n_atoms);
    }
    if (segid_as_auth_segid) out.auth_segid.reserve(n_atoms);
    if (output_break_records) out.auth_break.reserve(n_atoms);

    std::string const radix(std::localeconv()->decimal_point);
    static const std::string atom_record("ATOM");
    static const std::string hetatm_record("HETATM");
    static const std::string null_value(".");
    static const std::string unknown_value("?");
    // _atom_site.id is the key of the loop. The serial restarts in every
    // model, so in a multi-model hierarchy the id is the row number; the
    // anisotrop id is the same value.
    bool multi_model = models.size() > 1;
    long i_row = 0;
    std::size_t i_chain_all = 0;
    std::size_t i_rg_all = 0;
    std::size_t i_ag_all = 0;
    std::set<std::string> chem_comp_ids;
    std::set<std::string> struct_asym_ids;
    for (std::size_t i_model = 0; i_model < models.size(); i_model++) {
      hierarchy::model const& model = models[i_model];
      std::string model_num = model.data->id.empty()
        ? std::string("1") : python_strip(model.data->id);
      // Set once per model, not per chain: the first residue group of
      // every chain but a model's first gets auth_break 1 when it is not
      // linked to the previous one (the Python reference does this).
      bool is_first_in_chain = true;
      std::vector<hierarchy::chain> const& chains = model.chains();
      for (std::size_t i_chain = 0; i_chain < chains.size(); i_chain++) {
        std::string const&
          auth_asym_id = labels.chain_auth_asym_id[i_chain_all++];
        std::vector<hierarchy::residue_group> const&
          rgs = chains[i_chain].residue_groups();
        for (std::size_t i_rg = 0; i_rg < rgs.size(); i_rg++) {
          hierarchy::residue_group const& rg = rgs[i_rg];
          bool is_first_after_break = !(
            is_first_in_chain || rg.data->link_to_previous);
          std::string const& label_asym_id = labels.rg_label_asym_id[i_rg_all];
          std::string const& auth_seq_id = labels.rg_auth_seq_id[i_rg_all];
          i_rg_all++;
          std::string icode(rg.data->icode.elems);
          if (icode == " " || icode.empty()) icode = unknown_value;
          std::vector<hierarchy::atom_group> const& ags = rg.atom_groups();
          for (std::size_t i_ag = 0; i_ag < ags.size(); i_ag++) {
            hierarchy::atom_group const& ag = ags[i_ag];
            std::string comp_id = python_strip(ag.data->resname);
            std::string const& label_seq_id = labels.ag_label_seq_id[i_ag_all];
            i_ag_all++;
            std::string const& label_alt_id = ag.data->altloc.empty()
              ? null_value : ag.data->altloc;
            std::vector<hierarchy::atom> const& atoms = ag.atoms();
            if (!atoms.empty()) {
              chem_comp_ids.insert(comp_id);
              if (struct_asym_ids.insert(label_asym_id).second) {
                out.struct_asym_ids.push_back(label_asym_id);
              }
            }
            for (std::size_t i_atom = 0; i_atom < atoms.size(); i_atom++) {
              hierarchy::atom const& atom = atoms[i_atom];
              hierarchy::atom_data const& a = *atom.data;
              i_row++;
              std::string atom_id = decimal(
                multi_model ? i_row : static_cast<long>(atom.serial_as_int()));
              std::string name = python_strip(a.name.elems);
              out.group_pdb.push_back(a.hetero ? hetatm_record : atom_record);
              out.id.push_back(atom_id);
              out.label_atom_id.push_back(name);
              out.auth_atom_id.push_back(name);
              out.label_alt_id.push_back(label_alt_id);
              out.label_comp_id.push_back(comp_id);
              out.auth_asym_id.push_back(auth_asym_id);
              out.auth_seq_id.push_back(auth_seq_id);
              out.ins_code.push_back(icode);
              out.cartn_x.push_back(
                format_fixed(a.xyz[0], coordinate_precision, radix));
              out.cartn_y.push_back(
                format_fixed(a.xyz[1], coordinate_precision, radix));
              out.cartn_z.push_back(
                format_fixed(a.xyz[2], coordinate_precision, radix));
              out.occupancy.push_back(
                format_fixed(a.occ, occupancy_precision, radix));
              out.b_iso.push_back(format_fixed(a.b, b_iso_precision, radix));
              out.type_symbol.push_back(python_strip(a.element.elems));
              out.formal_charge.push_back(formal_charge(atom));
              if (a.fp == 0 && a.fdp == 0) {
                out.fp.push_back(null_value);
                out.fdp.push_back(null_value);
              }
              else {
                out.fp.push_back(format_fixed(a.fp, 4, radix));
                out.fdp.push_back(format_fixed(a.fdp, 4, radix));
              }
              out.resolution.push_back(
                format_fixed(a.resolution, resolution_precision, radix));
              out.label_asym_id.push_back(label_asym_id);
              out.label_entity_id.push_back(unknown_value);
              out.label_seq_id.push_back(label_seq_id);
              out.model_num.push_back(model_num);
              if (segid_as_auth_segid) {
                out.auth_segid.push_back(std::string(a.segid.elems));
              }
              if (output_break_records) {
                out.auth_break.push_back(is_first_after_break ? "1" : "0");
              }
              if (atom.uij_is_defined()) {
                out.aniso_id.push_back(atom_id);
                out.aniso_auth_atom_id.push_back(name);
                out.aniso_label_alt_id.push_back(label_alt_id);
                out.aniso_auth_comp_id.push_back(comp_id);
                out.aniso_auth_asym_id.push_back(auth_asym_id);
                out.aniso_auth_seq_id.push_back(auth_seq_id);
                out.aniso_ins_code.push_back(icode);
                out.aniso_u11.push_back(
                  format_fixed(a.uij[0], u_aniso_precision, radix));
                out.aniso_u22.push_back(
                  format_fixed(a.uij[1], u_aniso_precision, radix));
                out.aniso_u33.push_back(
                  format_fixed(a.uij[2], u_aniso_precision, radix));
                out.aniso_u12.push_back(
                  format_fixed(a.uij[3], u_aniso_precision, radix));
                out.aniso_u13.push_back(
                  format_fixed(a.uij[4], u_aniso_precision, radix));
                out.aniso_u23.push_back(
                  format_fixed(a.uij[5], u_aniso_precision, radix));
              }
              is_first_in_chain = false;
              is_first_after_break = false;
            }
          }
        }
      }
    }
    for (std::set<std::string>::const_iterator i = chem_comp_ids.begin();
         i != chem_comp_ids.end(); i++) {
      out.chem_comp_ids.push_back(*i);
    }
    return out;
  }

}}} // namespace iotbx::pdb::mmcif
