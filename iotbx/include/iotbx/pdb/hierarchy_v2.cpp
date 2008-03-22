#include <iotbx/pdb/hierarchy_v2.h>
#include <iotbx/pdb/common_residue_names.h>
#include <cctbx/eltbx/chemical_elements.h>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb { namespace hierarchy_v2 {

namespace {

  template <typename ParentType, typename ChildType>
  void
  detach_copy_children(
    ParentType const& new_parent,
    std::vector<ChildType>& new_children,
    std::vector<ChildType> const& old_children)
  {
    unsigned n = static_cast<unsigned>(old_children.size());
    if (n != 0) {
      new_children.reserve(n);
      const ChildType* s = &*old_children.begin();
      for(unsigned i=0;i<n;i++) {
        new_children.push_back(ChildType(new_parent, *s++));
      }
    }
  }

} // namespace <anonymous>

#define IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_SET(P, T) \
  T& \
  T::set_parent(P const& parent) \
  { \
    if (data->parent.lock().get() != 0) { \
      throw std::runtime_error(#T " has another parent " #P " already."); \
    } \
    data->parent = parent.data; \
    return *this; \
  }

#define IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET(P, T) \
  boost::optional<P> \
  T::parent() const \
  { \
    shared_ptr<P##_data> parent = data->parent.lock(); \
    if (parent.get() == 0) return boost::optional<P>(); \
    return boost::optional<P>(P(parent, true)); \
  }

#define IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET_SET(P, T) \
  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_SET(P, T) \
  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET(P, T)

#define IOTBX_PDB_HIERARCHY_V2_CPP_APPEND_ETC(T, C) \
  unsigned \
  T::C##s_size() const \
  { \
    return static_cast<unsigned>(data->C##s.size()); \
  } \
\
  std::vector<C> const& \
  T::C##s() const { return data->C##s; } \
\
  long \
  T::find_##C##_index( \
    hierarchy_v2::C const& C, \
    bool must_be_present) const \
  { \
    long n = static_cast<long>(data->C##s.size()); \
    for(long i=0;i<n;i++) { \
      if (data->C##s[i].data.get() == C.data.get()) return i; \
    } \
    if (must_be_present) { \
      throw std::runtime_error(#C " not in " #T "."); \
    } \
    return -1; \
  } \
\
  void \
  T::pre_allocate_##C##s(unsigned number_of_additional_##C##s) \
  { \
    data->C##s.reserve(data->C##s.size()+number_of_additional_##C##s); \
  } \
\
  void \
  T::insert_##C(long i, hierarchy_v2::C& C) \
  { \
    data->C##s.insert( \
      data->C##s.begin() \
        + positive_getitem_index(i, data->C##s.size(), true), \
      C.set_parent(*this)); \
  } \
\
  void \
  T::append_##C(hierarchy_v2::C& C) \
  { \
    data->C##s.push_back(C.set_parent(*this)); \
  } \
\
  void \
  T::remove_##C(long i) \
  { \
    std::size_t j =  positive_getitem_index(i, data->C##s.size()); \
    data->C##s[j].clear_parent(); \
    data->C##s.erase(data->C##s.begin() + j); \
  } \
\
  void \
  T::remove_##C(hierarchy_v2::C& C) \
  { \
    data->C##s.erase(data->C##s.begin() + find_##C##_index(C, true)); \
    C.clear_parent(); \
  }

#define IOTBX_PDB_HIERARCHY_V2_CPP_DETACHED_COPY_ETC(P, T, C) \
  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET_SET(P, T) \
  T::T( \
    P const& parent, \
    T const& other) \
  : \
    data(new T##_data(parent.data, *other.data)) \
  { \
    detach_copy_children(*this, data->C##s, other.data->C##s); \
  } \
\
  T \
  T::detached_copy() const \
  { \
    T result(data.get()); \
    detach_copy_children(result, result.data->C##s, data->C##s); \
    return result; \
  } \
  IOTBX_PDB_HIERARCHY_V2_CPP_APPEND_ETC(T, C)

  IOTBX_PDB_HIERARCHY_V2_CPP_APPEND_ETC(root, model)
  IOTBX_PDB_HIERARCHY_V2_CPP_DETACHED_COPY_ETC(root, model, chain)
  IOTBX_PDB_HIERARCHY_V2_CPP_DETACHED_COPY_ETC(model, chain, residue_group)
  IOTBX_PDB_HIERARCHY_V2_CPP_DETACHED_COPY_ETC(chain, residue_group, atom_group)
  IOTBX_PDB_HIERARCHY_V2_CPP_DETACHED_COPY_ETC(residue_group, atom_group, atom)
  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET_SET(atom_group, atom)

  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET(chain, conformer)
  IOTBX_PDB_HIERARCHY_V2_CPP_PARENT_GET(conformer, residue)

  root
  root::deep_copy() const
  {
    root result;
    result.data->info = data->info.deep_copy();
    detach_copy_children(result, result.data->models, data->models);
    return result;
  }

  atom_data::atom_data(
    weak_ptr<atom_group_data> const& parent_,
    atom_data const& other)
  :
    parent(parent_),
    name(other.name),
    segid(other.segid),
    element(other.element),
    charge(other.charge),
    serial(other.serial),
    xyz(other.xyz), sigxyz(other.sigxyz),
    occ(other.occ), sigocc(other.sigocc),
    b(other.b), sigb(other.sigb),
    uij(other.uij), siguij(other.siguij),
    hetero(other.hetero),
    tmp(0)
  {}

  atom
  atom::detached_copy() const
  {
    return atom(
      data->name.elems, data->segid.elems,
      data->element.elems, data->charge.elems, data->serial.elems,
      data->xyz, data->sigxyz,
      data->occ, data->sigocc,
      data->b, data->sigb,
      data->uij, data->siguij,
      data->hetero);
  }

  unsigned
  root::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_V2_CPP_ROOT_ATOM_GROUPS_LOOPS \
    std::vector<model> const& models = this->models(); \
    unsigned n_mds = models_size(); \
    for(unsigned i_md=0;i_md<n_mds;i_md++) { \
      unsigned n_chs = models[i_md].chains_size(); \
      std::vector<chain> const& chains = models[i_md].chains(); \
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) { \
        unsigned n_rgs = chains[i_ch].residue_groups_size(); \
        std::vector<residue_group> const& rgs = chains[i_ch].residue_groups(); \
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) { \
          unsigned n_ags = rgs[i_rg].atom_groups_size(); \
          std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
            unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_V2_CPP_ROOT_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}}}
    return result;
  }

  unsigned
  model::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_V2_CPP_MODEL_ATOM_GROUPS_LOOPS \
    unsigned n_chs = chains_size(); \
    std::vector<chain> const& chains = this->chains(); \
    for(unsigned i_ch=0;i_ch<n_chs;i_ch++) { \
      unsigned n_rgs = chains[i_ch].residue_groups_size(); \
      std::vector<residue_group> const& rgs = chains[i_ch].residue_groups(); \
      for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) { \
        unsigned n_ags = rgs[i_rg].atom_groups_size(); \
        std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
        for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
          unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_V2_CPP_MODEL_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}}
    return result;
  }

  unsigned
  chain::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_V2_CPP_CHAIN_ATOM_GROUPS_LOOPS \
    unsigned n_rgs = residue_groups_size(); \
    std::vector<residue_group> const& rgs = residue_groups(); \
    for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) { \
      unsigned n_ags = rgs[i_rg].atom_groups_size(); \
      std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
      for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
        unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_V2_CPP_CHAIN_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}
    return result;
  }

  unsigned
  residue_group::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_V2_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS \
    unsigned n_ags = atom_groups_size(); \
    std::vector<atom_group> const& ags = atom_groups(); \
    for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
      unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_V2_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS
      result += n_ats;
    }
    return result;
  }

  unsigned
  conformer::atoms_size() const
  {
    unsigned result = 0;
    unsigned n_rds = residues_size();
    std::vector<residue> const& rds = residues();
    for(unsigned i_rd=0;i_rd<n_rds;i_rd++) {
      unsigned n_ats = rds[i_rd].atoms_size();
      result += n_ats;
    }
    return result;
  }

  af::shared<atom>
  root::atoms() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_V2_CPP_ROOT_ATOM_GROUPS_LOOPS
#define IOTBX_PDB_HIERARCHY_V2_CPP_ATOM_GROUP_PUSH_BACK_LOOP \
      std::vector<atom> const& ats = ags[i_ag].atoms(); \
      for(unsigned i_at=0;i_at<n_ats;i_at++) { \
        result.push_back(ats[i_at]); \
      }
      IOTBX_PDB_HIERARCHY_V2_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}}}
    return result;
  }

  af::shared<atom>
  model::atoms() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_V2_CPP_MODEL_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_V2_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}}
    return result;
  }

  af::shared<atom>
  chain::atoms() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_V2_CPP_CHAIN_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_V2_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}
    return result;
  }

  af::shared<atom>
  residue_group::atoms() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_V2_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_V2_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }
    return result;
  }

  af::shared<atom>
  conformer::atoms() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    unsigned n_rds = residues_size();
    std::vector<residue> const& rds = residues();
    for(unsigned i_rd=0;i_rd<n_rds;i_rd++) {
      unsigned n_ats = rds[i_rd].atoms_size();
      std::vector<atom> const& ats = rds[i_rd].atoms();
      for(unsigned i_at=0;i_at<n_ats;i_at++) {
        result.push_back(ats[i_at]);
      }
    }
    return result;
  }

  void
  atom_label_columns_formatter::format(
    char* result) const
  {
    char blank = ' ';
    copy_left_justified(result, 4U, name, 4U, blank);
    copy_left_justified(result+4, 1U, altloc, 1U, blank);
    copy_right_justified(result+5, 3U, resname, 3U, blank);
    copy_right_justified(result+8, 2U, chain_id, 2U, blank);
    copy_right_justified(result+10, 4U, resseq, 4U, blank);
    copy_left_justified(result+14, 1U, icode, 1U, blank);
  }

  void
  atom_label_columns_formatter::format(
    char* result,
    hierarchy_v2::atom const& atom)
  {
    name = atom.data->name.elems;
    shared_ptr<atom_group_data> ag_lock = atom.data->parent.lock();
    const atom_group_data* ag = ag_lock.get();
    if (ag == 0) {
      altloc = resname = resseq = icode = chain_id = 0;
      format(result);
    }
    else {
      altloc = ag->altloc.elems;
      resname = ag->resname.elems;
      shared_ptr<residue_group_data> rg_lock = ag->parent.lock();
      const residue_group_data* rg = rg_lock.get();
      if (rg == 0) {
        resseq = icode = chain_id = 0;
        format(result);
      }
      else {
        resseq = rg->resseq.elems;
        icode = rg->icode.elems;
        shared_ptr<chain_data> c_lock = rg->parent.lock();
        chain_data const* c = c_lock.get();
        chain_id = (c == 0 ? 0 : c->id.c_str());
        format(result);
      }
    }
  }

  void
  atom::format_atom_record_serial_label_columns(
    char* result,
    atom_label_columns_formatter* label_formatter) const
  {
    char blank = ' ';
    data->serial.copy_right_justified(result+6, 5U, blank);
    result[11] = blank;
    if (label_formatter == 0) {
      atom_label_columns_formatter().format(result+12, *this);
    }
    else {
      label_formatter->name = data->name.elems;
      label_formatter->format(result+12);
    }
  }

  unsigned
  atom::format_atom_record_segid_element_charge_columns(
    char* result,
    unsigned blanks_start_at) const
  {
    char blank = ' ';
    data->segid.copy_left_justified(result+72, 4U, blank);
    data->element.copy_right_justified(result+76, 2U, blank);
    data->charge.copy_left_justified(result+78, 2U, blank);
    for(unsigned i=79;i!=71;i--) {
      if (result[i] != blank) {
        copy_left_justified(
          result+blanks_start_at, 72U-blanks_start_at, 0, 0U, blank);
        result[i+1] = '\0';
        return i+1;
      }
    }
    result[blanks_start_at] = '\0';
    return blanks_start_at;
  }

  void
  atom::format_pdb_element_charge_columns(
    char* result) const
  {
    char blank = ' ';
    data->element.copy_right_justified(result, 2U, blank);
    data->charge.copy_left_justified(result+2, 2U, blank);
  }

  std::string
  atom::pdb_label_columns() const
  {
    char result[15];
    atom_label_columns_formatter().format(result, *this);
    return std::string(result, 15U);
  }

  small_str<15>
  atom::pdb_label_columns_small_str() const
  {
    small_str<15> result;
    atom_label_columns_formatter().format(result.elems, *this);
    result.elems[15] = '\0';
    return result;
  }

  std::string
  atom::pdb_element_charge_columns() const
  {
    char result[4];
    format_pdb_element_charge_columns(result);
    return std::string(result, 4U);
  }

  unsigned
  atom::format_atom_record(
    char* result,
    atom_label_columns_formatter* label_formatter,
    bool cut_after_label_columns) const
  {
    char blank = ' ';
    std::memcpy(result, (data->hetero ? "HETATM" : "ATOM  "), 6U);
    format_atom_record_serial_label_columns(result, label_formatter);
    if (cut_after_label_columns) {
      result[27U] = '\0';
      return 27U;
    }
    copy_left_justified(result+27, 3U, 0, 0U, blank);
    char *r = result + 30;
    for(unsigned i=0;i<3;i++) {
      std::sprintf(r, "%8.3f", data->xyz[i]);
      r += 8;
      if (*r != '\0') {
        throw std::runtime_error(
          std::string("atom ") + "XYZ"[i] + " coordinate value"
          " does not fit into F8.3 format.");
      }
    }
    std::sprintf(r, "%6.2f", data->occ);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom ") + "occupancy factor"
        " does not fit into F6.2 format.");
    }
    std::sprintf(r, "%6.2f", data->b);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom ") + "B-factor"
        " does not fit into F6.2 format.");
    }
    return format_atom_record_segid_element_charge_columns(result, 66U);
  }

  unsigned
  atom::format_sigatm_record(
    char* result,
    atom_label_columns_formatter* label_formatter) const
  {
    char blank = ' ';
    std::memcpy(result, "SIGATM", 6U);
    format_atom_record_serial_label_columns(result, label_formatter);
    copy_left_justified(result+27, 3U, 0, 0U, blank);
    char *r = result + 30;
    for(unsigned i=0;i<3;i++) {
      std::sprintf(r, "%8.3f", data->sigxyz[i]);
      r += 8;
      if (*r != '\0') {
        throw std::runtime_error(
          std::string("atom sigma ") + "XYZ"[i] + " coordinate value"
          " does not fit into F8.3 format.");
      }
    }
    std::sprintf(r, "%6.2f", data->sigocc);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom sigma ") + "occupancy factor"
        " does not fit into F6.2 format.");
    }
    std::sprintf(r, "%6.2f", data->sigb);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom sigma ") + "B-factor"
        " does not fit into F6.2 format.");
    }
    return format_atom_record_segid_element_charge_columns(result, 66U);
  }

  unsigned
  atom::format_anisou_record(
    char* result,
    atom_label_columns_formatter* label_formatter) const
  {
    char blank = ' ';
    std::memcpy(result, "ANISOU", 6U);
    format_atom_record_serial_label_columns(result, label_formatter);
    result[27] = blank;
    char *r = result + 28;
    for(unsigned i=0;i<6;i++) {
      std::sprintf(r, "%7.0f", data->uij[i]*10000.);
      r += 7;
      if (*r != '\0') {
        static const char* uij_labels[] = {
          "U11", "U22", "U33", "U12", "U13", "U23"};
        throw std::runtime_error(
          std::string("atom ") + uij_labels[i] + " value * 10000"
          " does not fit into F7.0 format.");
      }
    }
    return format_atom_record_segid_element_charge_columns(result, 70U);
  }

  unsigned
  atom::format_siguij_record(
    char* result,
    atom_label_columns_formatter* label_formatter) const
  {
    char blank = ' ';
    std::memcpy(result, "SIGUIJ", 6U);
    format_atom_record_serial_label_columns(result, label_formatter);
    result[27] = blank;
    char *r = result + 28;
    for(unsigned i=0;i<6;i++) {
      std::sprintf(r, "%7.0f", data->siguij[i]*10000.);
      r += 7;
      if (*r != '\0') {
        static const char* uij_labels[] = {
          "U11", "U22", "U33", "U12", "U13", "U23"};
        throw std::runtime_error(
          std::string("atom sigma ") + uij_labels[i] + " value * 10000"
          " does not fit into F7.0 format.");
      }
    }
    return format_atom_record_segid_element_charge_columns(result, 70U);
  }

  unsigned
  atom::format_atom_record_group(
    char* result,
    atom_label_columns_formatter* label_formatter,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij) const
  {
    char newline = '\n';
    unsigned str_len = 0;
    if (atom_hetatm) {
      str_len += format_atom_record(result, label_formatter);
    }
    if (  sigatm
        && (  !data->sigxyz.const_ref().all_le(0)
            || data->sigocc > 0
            || data->sigb > 0)) {
      if (str_len != 0) result[str_len++] = newline;
      str_len += format_sigatm_record(result+str_len, label_formatter);
    }
    if (anisou && uij_is_defined()) {
      if (str_len != 0) result[str_len++] = newline;
      str_len += format_anisou_record(result+str_len, label_formatter);
    }
    if (siguij && siguij_is_defined()) {
      if (str_len != 0) result[str_len++] = newline;
      str_len += format_siguij_record(result+str_len, label_formatter);
    }
    if (str_len == 0) result[0] = '\0';
    return str_len;
  }

  bool
  atom::element_is_hydrogen() const
  {
    const char* e = data->element.elems;
    return (
         (e[0] == ' ' && (e[1] == 'H' || e[1] == 'D'))
      || ((e[0] == 'H' || e[0] == 'D') && (e[1] == ' ' || e[1] == '\0')));
  }

  boost::optional<std::string>
  atom::determine_chemical_element_simple() const
  {
    std::set<std::string> const& chemical_elements
      = cctbx::eltbx::chemical_elements::proper_and_isotopes_upper_set();
    std::string e(data->element.elems);
    std::string l;
    if (e[0] == ' ') {
      l = e[1];
    }
    else {
      e[1] = toupper(e[1]); // first character must be upper case, but
                            // second character may be upper or lower case
      l = e;
    }
    if (chemical_elements.find(l) != chemical_elements.end()) {
      return boost::optional<std::string>(e);
    }
    if (e == "  ") {
      std::string n(data->name.elems, 2);
      if (isdigit(n[0])) n[0] = ' ';
      if (n[0] == ' ') l = n[1];
      else             l = n;
      if (chemical_elements.find(l) != chemical_elements.end()) {
        return boost::optional<std::string>(n);
      }
    }
    return boost::optional<std::string>();
  }

  std::string
  atom_group::confid() const
  {
    char blank = ' ';
    char result[5];
    data->altloc.copy_left_justified(result, 1U, blank);
    data->resname.copy_right_justified(result+1, 3U, blank);
    result[4] = '\0';
    return std::string(result);
  }

namespace {

  template <typename DataType>
  std::string
  make_resid(DataType const& data)
  {
    char blank = ' ';
    char result[6];
    data->resseq.copy_right_justified(result, 4U, blank);
    data->icode.copy_left_justified(result+4, 1U, blank);
    result[5] = '\0';
    return std::string(result);
  }

} // namespace <anonymous>

  std::string
  residue_group::resid() const { return make_resid(data); }

  str5
  residue_group::resid_small_str() const
  {
    char blank = ' ';
    str5 result;
    data->resseq.copy_right_justified(result.elems, 4U, blank);
    data->icode.copy_left_justified(result.elems+4, 1U, blank);
    result.elems[5] = '\0';
    return result;
  }

  std::string
  residue::resid() const { return make_resid(data); }

  bool
  residue_group::have_conformers() const
  {
    typedef std::vector<atom_group>::const_iterator agi_t;
    agi_t ag_end = data->atom_groups.end();
    for(agi_t agi=data->atom_groups.begin();agi!=ag_end;agi++) {
      char altloc = agi->data->altloc.elems[0];
      if (altloc != '\0' && altloc != blank_altloc_char) {
        return true;
      }
    }
    return false;
  }

  void
  residue_group::merge_atom_groups(
    atom_group& primary,
    atom_group& secondary)
  {
    SCITBX_ASSERT(secondary.data->altloc == primary.data->altloc);
    SCITBX_ASSERT(secondary.data->resname == primary.data->resname);
    if (primary.parent_ptr().get() != data.get()) {
      throw std::runtime_error(
        "\"primary\" atom_group has a different or no parent"
        " (this residue_group must be the parent).");
    }
    if (secondary.data.get() == primary.data.get()) {
      throw std::runtime_error(
        "\"primary\" and \"secondary\" atom_groups are identical.");
    }
    unsigned n_atoms = secondary.atoms_size();
    boost::scoped_array<atom> atom_buffer(new atom[n_atoms]);
    for(unsigned i=0;i<n_atoms;i++) { // copy to buffer ...
      atom_buffer[i] = secondary.atoms()[i];
    }
    for(unsigned i=n_atoms;i!=0;) { // ... so we can remove from end ...
      secondary.remove_atom(--i);
    }
    primary.pre_allocate_atoms(n_atoms);
    for(unsigned i=0;i<n_atoms;i++) { // ... and then append in original order.
      primary.append_atom(atom_buffer[i]);
    }
    SCITBX_ASSERT(secondary.atoms_size() == 0);
  }

  unsigned
  residue_group::move_blank_altloc_atom_groups_to_front()
  {
    unsigned n_blank_altloc_atom_groups = 0;
    unsigned n_ag = atom_groups_size();
    for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
      atom_group const& ag = data->atom_groups[i_ag];
      char altloc = ag.data->altloc.elems[0];
      if (altloc == '\0' || altloc == blank_altloc_char) {
        if (i_ag != n_blank_altloc_atom_groups) {
          atom_group ag_by_value = ag;
          remove_atom_group(i_ag);
          insert_atom_group(n_blank_altloc_atom_groups, ag_by_value);
        }
        n_blank_altloc_atom_groups++;
      }
    }
    return n_blank_altloc_atom_groups;
  }

  af::tiny<unsigned, 2>
  residue_group::edit_blank_altloc()
  {
    unsigned
      n_blank_altloc_atom_groups = move_blank_altloc_atom_groups_to_front();
    if (n_blank_altloc_atom_groups == 0) {
      return af::tiny<unsigned, 2>(0, 0);
    }
    typedef std::set<str4> ss4;
    typedef std::map<str3, ss4> ms3ss4;
    ms3ss4 blank_name_sets;
    unsigned i_ag = 0;
    for(;i_ag<n_blank_altloc_atom_groups;i_ag++) {
      atom_group& ag = data->atom_groups[i_ag];
      ag.data->altloc.elems[0] = '\0';
      ss4& blank_name_set = blank_name_sets[ag.data->resname];
      unsigned n_atoms = ag.atoms_size();
      for(unsigned i_atom=0;i_atom<n_atoms;i_atom++) {
        blank_name_set.insert(ag.atoms()[i_atom].data->name);
      }
    }
    ms3ss4 blank_but_alt_sets;
    ms3ss4::const_iterator blank_name_sets_end = blank_name_sets.end();
    unsigned n_ag = atom_groups_size();
    for(;i_ag<n_ag;i_ag++) {
      atom_group& ag = data->atom_groups[i_ag];
      ms3ss4::const_iterator blank_name_sets_iter = blank_name_sets.find(
        ag.data->resname);
      if (blank_name_sets_iter == blank_name_sets_end) continue;
      ss4 const& blank_name_set = blank_name_sets_iter->second;
      ss4::const_iterator blank_name_set_end = blank_name_set.end();
      ss4* blank_but_alt_resname = 0;
      unsigned n_atoms = ag.atoms_size();
      for(unsigned i_atom=0;i_atom<n_atoms;i_atom++) {
        str4 const& atom_name = ag.atoms()[i_atom].data->name;
        ss4::const_iterator blank_name_set_iter = blank_name_set.find(
          atom_name);
        if (blank_name_set_iter == blank_name_set_end) continue;
        if (blank_but_alt_resname == 0) {
          blank_but_alt_resname = &blank_but_alt_sets[ag.data->resname];
        }
        blank_but_alt_resname->insert(atom_name);
      }
    }
    unsigned n_blank_but_alt_atom_groups = 0;
    if (blank_but_alt_sets.size() != 0) {
      ms3ss4::const_iterator blank_but_alt_sets_end = blank_but_alt_sets.end();
      for(unsigned i_ag=0;i_ag<n_blank_altloc_atom_groups;) {
        atom_group ag = data->atom_groups[i_ag];
        ms3ss4::const_iterator
          blank_but_alt_sets_iter = blank_but_alt_sets.find(ag.data->resname);
        if (blank_but_alt_sets_iter == blank_but_alt_sets_end) {
          i_ag++;
          continue;
        }
        ss4::const_iterator
          blank_but_alt_set_end = blank_but_alt_sets_iter->second.end();
        atom_group* new_atom_group = 0;
        unsigned n_atoms = ag.atoms_size();
        for(unsigned i_atom=0;i_atom<n_atoms;) {
          hierarchy_v2::atom atom = ag.atoms()[i_atom];
          ss4::const_iterator
            blank_but_alt_set_iter = blank_but_alt_sets_iter->second.find(
              atom.data->name);
          if (blank_but_alt_set_iter == blank_but_alt_set_end) {
            i_atom++;
            continue;
          }
          if (new_atom_group == 0) {
            unsigned i = n_blank_altloc_atom_groups
                       + n_blank_but_alt_atom_groups;
            atom_group new_ag(blank_altloc_cstr, ag.data->resname.elems);
            insert_atom_group(i, new_ag);
            new_atom_group = &data->atom_groups[i];
            n_blank_but_alt_atom_groups++;
          }
          ag.remove_atom(i_atom);
          n_atoms--;
          new_atom_group->append_atom(atom);
        }
        if (ag.atoms_size() == 0) {
          remove_atom_group(i_ag);
          n_blank_altloc_atom_groups--;
        }
        else {
          i_ag++;
        }
      }
    }
    return af::tiny<unsigned, 2>(
      n_blank_altloc_atom_groups,
      n_blank_but_alt_atom_groups);
  }

  void
  chain::merge_residue_groups(
    residue_group& primary,
    residue_group& secondary)
  {
    SCITBX_ASSERT(secondary.data->resseq == primary.data->resseq);
    SCITBX_ASSERT(secondary.data->icode == primary.data->icode);
    const chain_data* data_get = data.get();
    if (primary.parent_ptr().get() != data_get) {
      throw std::runtime_error(
        "\"primary\" residue_group has a different or no parent"
        " (this chain must be the parent).");
    }
    if (secondary.parent_ptr().get() != data_get) {
      throw std::runtime_error(
        "\"secondary\" residue_group has a different or no parent"
        " (this chain must be the parent).");
    }
    typedef std::map<str3, atom_group> s3ag;
    typedef std::map<str1, s3ag> s1s3ag;
    s1s3ag altloc_resname_dict;
    unsigned n_ag = primary.atom_groups_size();
    for(unsigned i=0;i<n_ag;i++) {
      atom_group const& ag = primary.atom_groups()[i];
      altloc_resname_dict[ag.data->altloc][ag.data->resname] = ag;
    }
    n_ag = secondary.atom_groups_size();
    std::vector<atom_group> append_buffer;
    append_buffer.reserve(n_ag);// allocate potential max to avoid reallocation
    for(unsigned i=n_ag;i!=0;) {
      i--;
      atom_group ag = secondary.atom_groups()[i];
      bool have_primary = false;
      s1s3ag::const_iterator altloc_iter = altloc_resname_dict.find(
        ag.data->altloc);
      if (altloc_iter != altloc_resname_dict.end()) {
        s3ag::const_iterator resname_iter = altloc_iter->second.find(
          ag.data->resname);
        if (resname_iter != altloc_iter->second.end()) {
          secondary.remove_atom_group(i);
          atom_group primary_ag = resname_iter->second;
          primary.merge_atom_groups(primary_ag, ag);
          have_primary = true;
        }
      }
      if (!have_primary) {
        secondary.remove_atom_group(i);
        append_buffer.push_back(ag);
      }
    }
    n_ag = static_cast<unsigned>(append_buffer.size());
    for(unsigned i=n_ag;i!=0;) {
      i--;
      primary.append_atom_group(append_buffer[i]);
    }
    SCITBX_ASSERT(secondary.atom_groups_size() == 0);
    remove_residue_group(secondary);
  }

  af::shared<std::size_t>
  chain::merge_disconnected_residue_groups_with_pure_altloc()
  {
    af::shared<std::size_t> result;
    typedef std::map<str1, std::vector<unsigned> > s1i;
    typedef std::map<str4, s1i> s4s1i;
    s4s1i matching_resseq_icode;
    unsigned n_rg = residue_groups_size();
    for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
      residue_group const& rg = data->residue_groups[i_rg];
      matching_resseq_icode[rg.data->resseq][rg.data->icode].push_back(i_rg);
    }
    std::vector<residue_group> residue_groups_copy;
    std::vector<unsigned> flags; // 0: not modified, 1: primary, 2: secondary
    unsigned result_size = 0;
    s4s1i::const_iterator mri_end = matching_resseq_icode.end();
    for(s4s1i::const_iterator
          mri=matching_resseq_icode.begin();mri!=mri_end;mri++) {
      s1i::const_iterator mrii_end = mri->second.end();
      for(s1i::const_iterator
            mrii=mri->second.begin();mrii!=mrii_end;mrii++) {
        std::vector<unsigned> const& i_rgs = mrii->second;
        unsigned n_i_rgs = static_cast<unsigned>(i_rgs.size());
        if (n_i_rgs == 1U) continue;
        std::set<str1> altlocs;
        altlocs.insert('\0');
        altlocs.insert(blank_altloc_char);
        unsigned altlocs_size = 2U;
        for(unsigned i_i_rgs=0;i_i_rgs<n_i_rgs;i_i_rgs++) {
          unsigned i_rg = i_rgs[i_i_rgs];
          residue_group const& rg = (
            result_size == 0 ? data->residue_groups[i_rg]
                             : residue_groups_copy[i_rg]);
          unsigned n_ag = rg.atom_groups_size();
          std::vector<atom_group> const& ags = rg.atom_groups();
          for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
            altlocs.insert(ags[i_ag].data->altloc);
            if (altlocs.size() == altlocs_size) {
              goto next_i_rgs;
            }
            altlocs_size++;
          }
        }
        { // merge
          if (result_size == 0) {
            residue_groups_copy = residue_groups();
            flags.resize(n_rg, 0);
          }
          unsigned i_rg = i_rgs[0];
          residue_group& primary = residue_groups_copy[i_rg];
          flags[i_rg] = 1U;
          for(unsigned i_i_rgs=1U;i_i_rgs<n_i_rgs;i_i_rgs++) {
            i_rg = i_rgs[i_i_rgs];
            residue_group& secondary = residue_groups_copy[i_rg];
            flags[i_rg] = 2U;
            merge_residue_groups(primary, secondary);
            SCITBX_ASSERT(secondary.atom_groups_size() == 0);
            SCITBX_ASSERT(secondary.parent_ptr().get() == 0);
          }
          result_size++;
        }
        next_i_rgs:;
      }
    }
    if (result_size != 0) {
      result.reserve(result_size);
      unsigned j_rg = 0;
      for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
        unsigned f = flags[i_rg];
        if      (f == 0) j_rg++;
        else if (f == 1U) result.push_back(j_rg++);
      }
      SCITBX_ASSERT(result.size() == result_size);
    }
    return result;
  }

  af::shared<af::tiny<std::size_t, 2> >
  chain::find_pure_altloc_ranges(
    const char* common_residue_name_class_only) const
  {
    af::shared<af::tiny<std::size_t, 2> > result;
    unsigned n_rg = residue_groups_size();
    unsigned range_start = n_rg;
    unsigned skip = 0;
    for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
      residue_group const& rg = data->residue_groups[i_rg];
      unsigned n_ags = rg.atom_groups_size();
      std::vector<atom_group> const& ags = rg.atom_groups();
      bool is_pure_altloc = (   n_ags != 0
                             && ags[0].data->altloc.elems[0] != '\0');
      if (common_residue_name_class_only != 0) {
        skip = 1;
        for(unsigned i_ag=0;i_ag<n_ags;i_ag++) {
          if (   common_residue_names::get_class(ags[i_ag].data->resname)
              == common_residue_name_class_only) {
            skip = 0;
            break;
          }
        }
      }
      if (skip != 0 || !rg.data->link_to_previous) {
        if (range_start+1U < i_rg) {
          result.push_back(af::tiny<std::size_t, 2>(range_start, i_rg));
        }
        range_start = (is_pure_altloc ? i_rg+skip : n_rg);
      }
      else {
        if (!is_pure_altloc) {
          if (range_start+1U < i_rg) {
            result.push_back(af::tiny<std::size_t, 2>(range_start, i_rg));
          }
          range_start = n_rg;
        }
        else if (range_start == n_rg) {
          range_start = i_rg + skip;
        }
      }
    }
    if (range_start+1U < n_rg) {
      result.push_back(af::tiny<std::size_t, 2>(range_start, n_rg));
    }
    return result;
  }

  void
  atoms_reset_tmp(
    af::const_ref<atom> const& atoms,
    int first_value,
    int increment)
  {
    int value = first_value;
    for(const atom* a=atoms.begin();a!=atoms.end();a++) {
      a->data->tmp = value;
      value += increment;
    }
  }

  void
  atoms_reset_tmp_for_occupancy_groups_simple(
    af::const_ref<atom> const& atoms)
  {
    int value = 0;
    for(const atom* a=atoms.begin();a!=atoms.end();a++,value++) {
      a->data->tmp = (a->element_is_hydrogen() ? -1 : value);
    }
  }

  residue::residue(
    conformer const& parent,
    const char* resname,
    const char* resseq,
    const char* icode,
    bool link_to_previous,
    bool is_pure_main_conf,
    std::vector<atom> const& atoms)
  :
    data(new residue_data(
      parent.data,
      resname, resseq, icode, link_to_previous, is_pure_main_conf,
      atoms))
  {}

  void
  conformer::append_residue(
    const char* resname,
    const char* resseq,
    const char* icode,
    bool link_to_previous,
    bool is_pure_main_conf,
    std::vector<atom> const& atoms)
  {
    data->residues.push_back(residue(
      *this,
      resname, resseq, icode, link_to_previous, is_pure_main_conf,
      atoms));
  }

  af::shared<conformer>
  chain::conformers() const
  {
    const char nulc = '\0';
    std::vector<char> altlocs;
    typedef std::map<char, unsigned> mcu;
    mcu altloc_indices;
    unsigned n_rg = residue_groups_size();
    for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
      residue_group const& rg = data->residue_groups[i_rg];
      unsigned n_ag = rg.atom_groups_size();
      std::vector<atom_group> const& ags = rg.atom_groups();
      for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
        char altloc = ags[i_ag].data->altloc.elems[0];
        if (altloc == nulc) continue;
        if (altloc_indices.find(altloc) == altloc_indices.end()) {
          altlocs.push_back(altloc);
          unsigned i = altloc_indices.size(); // important to get .size() ...
          altloc_indices[altloc] = i; // ... before using []
        }
      }
    }
    unsigned n_cf = static_cast<unsigned>(altloc_indices.size());
    if (n_cf == 0) {
      altlocs.push_back(nulc);
      n_cf++;
    }
    af::shared<conformer> result((af::reserve(n_cf)));
    for(unsigned i_cf=0;i_cf<n_cf;i_cf++) {
      result.push_back(conformer(*this, str1(altlocs[i_cf]).elems));
    }
    std::vector<str3> resnames; // allocate once
    resnames.reserve(32U); // not critical
    std::set<str3> resnames_with_altloc; // allocate once
    std::vector<std::vector<atom_group> > altloc_ags(n_cf); // allocate once
    for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
      residue_group const& rg = data->residue_groups[i_rg];
      for(unsigned i_cf=0;i_cf<n_cf;i_cf++) {
        // preserving allocations reduces re-allocations
        altloc_ags[i_cf].clear();
      }
      unsigned n_ag = rg.atom_groups_size();
      std::vector<atom_group> const& ags = rg.atom_groups();
      for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
        atom_group const& ag = ags[i_ag];
        char altloc = ag.data->altloc.elems[0];
        if (altloc == nulc) {
          for(unsigned i_cf=0;i_cf<n_cf;i_cf++) {
            altloc_ags[i_cf].push_back(ag);
          }
        }
        else {
          altloc_ags[altloc_indices[altloc]].push_back(ag);
        }
      }
      for(unsigned i_cf=0;i_cf<n_cf;i_cf++) {
        resnames.clear();
        typedef std::map<str3, std::vector<atom> > ms3va;
        ms3va resname_atoms;
        resnames_with_altloc.clear();
        std::vector<atom_group> const& ags = altloc_ags[i_cf];
        unsigned n_ag = static_cast<unsigned>(ags.size());
        for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
          atom_group const& ag = ags[i_ag];
          ms3va::iterator rai = resname_atoms.find(ag.data->resname);
          std::vector<atom>* atoms;
          if (rai != resname_atoms.end()) {
            atoms = &rai->second;
          }
          else {
            resnames.push_back(ag.data->resname);
            atoms = &resname_atoms[ag.data->resname];
            unsigned n_reserve = 100U / n_ag; // not critical
            if (n_reserve != 0) {
              atoms->reserve(n_reserve);
            }
          }
          atoms->insert(atoms->end(), ag.atoms().begin(), ag.atoms().end());
          if (ag.data->altloc.elems[0] != nulc) {
            resnames_with_altloc.insert(ag.data->resname);
          }
        }
        unsigned n_rn = static_cast<unsigned>(resnames.size());
        for(unsigned i_rn=0;i_rn<n_rn;i_rn++) {
          str3 const& resname = resnames[i_rn];
          ms3va::const_iterator rai = resname_atoms.find(resname);
          SCITBX_ASSERT(rai != resname_atoms.end());
          result[i_cf].append_residue(
            resname.elems,
            rg.data->resseq.elems,
            rg.data->icode.elems,
            rg.data->link_to_previous,
            /* is_pure_main_conf */ (   resnames_with_altloc.find(resname)
                                     == resnames_with_altloc.end()),
            rai->second);
        }
      }
    }
    return result;
  }

}}} // namespace iotbx::pdb::hierarchy_v2
