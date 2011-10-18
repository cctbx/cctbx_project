#include <iotbx/pdb/hierarchy.h>
#include <iotbx/pdb/common_residue_names.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <iotbx/error.h>
#include <cctbx/eltbx/chemical_elements.h>
#include <boost/format.hpp>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb { namespace hierarchy {

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

#define IOTBX_PDB_HIERARCHY_CPP_PARENT_SET(P, T) \
  T& \
  T::set_parent(P const& parent) \
  { \
    if (data->parent.lock().get() != 0) { \
      throw std::runtime_error(#T " has another parent " #P " already."); \
    } \
    data->parent = parent.data; \
    return *this; \
  }

#define IOTBX_PDB_HIERARCHY_CPP_PARENT_GET(P, T) \
  boost::optional<P> \
  T::parent(bool optional) const \
  { \
    shared_ptr<P##_data> parent = data->parent.lock(); \
    if (parent.get() == 0) { \
      if (optional) return boost::optional<P>(); \
      throw std::runtime_error(#T " has no parent " #P); \
    } \
    return boost::optional<P>(P(parent, true)); \
  }

#define IOTBX_PDB_HIERARCHY_CPP_PARENT_GET_SET(P, T) \
  IOTBX_PDB_HIERARCHY_CPP_PARENT_SET(P, T) \
  IOTBX_PDB_HIERARCHY_CPP_PARENT_GET(P, T)

#define IOTBX_PDB_HIERARCHY_CPP_APPEND_ETC(T, C) \
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
    hierarchy::C const& C, \
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
  T::insert_##C(long i, hierarchy::C& C) \
  { \
    data->C##s.insert( \
      data->C##s.begin() \
        + positive_getitem_index(i, data->C##s.size(), true), \
      C.set_parent(*this)); \
  } \
\
  void \
  T::append_##C(hierarchy::C& C) \
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
  T::remove_##C(hierarchy::C& C) \
  { \
    data->C##s.erase(data->C##s.begin() + find_##C##_index(C, true)); \
    C.clear_parent(); \
  }

#define IOTBX_PDB_HIERARCHY_CPP_DETACHED_COPY_ETC(P, T, C) \
  IOTBX_PDB_HIERARCHY_CPP_PARENT_GET_SET(P, T) \
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
  IOTBX_PDB_HIERARCHY_CPP_APPEND_ETC(T, C)

  IOTBX_PDB_HIERARCHY_CPP_APPEND_ETC(root, model)
  IOTBX_PDB_HIERARCHY_CPP_DETACHED_COPY_ETC(root, model, chain)
  IOTBX_PDB_HIERARCHY_CPP_DETACHED_COPY_ETC(model, chain, residue_group)
  IOTBX_PDB_HIERARCHY_CPP_DETACHED_COPY_ETC(chain, residue_group, atom_group)
  IOTBX_PDB_HIERARCHY_CPP_DETACHED_COPY_ETC(residue_group, atom_group, atom)
  IOTBX_PDB_HIERARCHY_CPP_PARENT_GET_SET(atom_group, atom)

  IOTBX_PDB_HIERARCHY_CPP_PARENT_GET(chain, conformer)
  IOTBX_PDB_HIERARCHY_CPP_PARENT_GET(conformer, residue)

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
    xyz(other.xyz), sigxyz(other.sigxyz),
    occ(other.occ), sigocc(other.sigocc),
    b(other.b), sigb(other.sigb),
    uij(other.uij),
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
    siguij(other.siguij),
#endif
    i_seq(0), tmp(0), have_sentinel(false),
    hetero(other.hetero), serial(other.serial), name(other.name),
    segid(other.segid), element(other.element), charge(other.charge)
  {}

  atom
  atom::detached_copy() const
  {
    return atom(
      data->xyz, data->sigxyz,
      data->occ, data->sigocc,
      data->b, data->sigb,
      data->uij,
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
      data->siguij,
#else
      sym_mat3(-1,-1,-1,-1,-1,-1),
#endif
      data->hetero, data->serial.elems, data->name.elems,
      data->segid.elems, data->element.elems, data->charge.elems);
  }

  unsigned
  root::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_CPP_ROOT_RESIDUE_GROUPS_LOOPS \
    std::vector<model> const& models = this->models(); \
    unsigned n_mds = models_size(); \
    for(unsigned i_md=0;i_md<n_mds;i_md++) { \
      unsigned n_chs = models[i_md].chains_size(); \
      std::vector<chain> const& chains = models[i_md].chains(); \
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) { \
        unsigned n_rgs = chains[i_ch].residue_groups_size(); \
        std::vector<residue_group> const& rgs = chains[i_ch].residue_groups(); \
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
#define IOTBX_PDB_HIERARCHY_CPP_ROOT_ATOM_GROUPS_LOOPS \
        IOTBX_PDB_HIERARCHY_CPP_ROOT_RESIDUE_GROUPS_LOOPS \
          unsigned n_ags = rgs[i_rg].atom_groups_size(); \
          std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
            unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_CPP_ROOT_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}}}
    return result;
  }

  unsigned
  model::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_CPP_MODEL_RESIDUE_GROUPS_LOOPS \
    unsigned n_chs = chains_size(); \
    std::vector<chain> const& chains = this->chains(); \
    for(unsigned i_ch=0;i_ch<n_chs;i_ch++) { \
      unsigned n_rgs = chains[i_ch].residue_groups_size(); \
      std::vector<residue_group> const& rgs = chains[i_ch].residue_groups(); \
      for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
#define IOTBX_PDB_HIERARCHY_CPP_MODEL_ATOM_GROUPS_LOOPS \
        IOTBX_PDB_HIERARCHY_CPP_MODEL_RESIDUE_GROUPS_LOOPS \
        unsigned n_ags = rgs[i_rg].atom_groups_size(); \
        std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
        for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
          unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_CPP_MODEL_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}}
    return result;
  }

  unsigned
  chain::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_CPP_CHAIN_RESIDUE_GROUPS_LOOPS \
    unsigned n_rgs = residue_groups_size(); \
    std::vector<residue_group> const& rgs = residue_groups(); \
    for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
#define IOTBX_PDB_HIERARCHY_CPP_CHAIN_ATOM_GROUPS_LOOPS \
        IOTBX_PDB_HIERARCHY_CPP_CHAIN_RESIDUE_GROUPS_LOOPS \
      unsigned n_ags = rgs[i_rg].atom_groups_size(); \
      std::vector<atom_group> const& ags = rgs[i_rg].atom_groups(); \
      for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
        unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_CPP_CHAIN_ATOM_GROUPS_LOOPS
      result += n_ats;
    }}
    return result;
  }

  unsigned
  residue_group::atoms_size() const
  {
    unsigned result = 0;
#define IOTBX_PDB_HIERARCHY_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS \
    unsigned n_ags = atom_groups_size(); \
    std::vector<atom_group> const& ags = atom_groups(); \
    for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
      unsigned n_ats = ags[i_ag].atoms_size();
    IOTBX_PDB_HIERARCHY_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS
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
  root::atoms_sequential_conf() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_ROOT_ATOM_GROUPS_LOOPS
#define IOTBX_PDB_HIERARCHY_CPP_ATOM_GROUP_PUSH_BACK_LOOP \
      std::vector<atom> const& ats = ags[i_ag].atoms(); \
      for(unsigned i_at=0;i_at<n_ats;i_at++) { \
        result.push_back(ats[i_at]); \
      }
      IOTBX_PDB_HIERARCHY_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}}}
    return result;
  }

  af::shared<atom>
  model::atoms_sequential_conf() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_MODEL_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}}
    return result;
  }

  af::shared<atom>
  chain::atoms_sequential_conf() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_CHAIN_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_CPP_ATOM_GROUP_PUSH_BACK_LOOP
    }}
    return result;
  }

  af::shared<atom>
  residue_group::atoms_sequential_conf() const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_RESIDUE_GROUP_ATOM_GROUPS_LOOPS
      IOTBX_PDB_HIERARCHY_CPP_ATOM_GROUP_PUSH_BACK_LOOP
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
    char* result,
    bool add_model,
    bool add_segid) const
  {
    char blank = ' ';
    if (add_model) {
      if (model_id != 0) {
        std::size_t l = std::strlen(model_id);
        IOTBX_ASSERT(l <= 8);
        unsigned n = static_cast<unsigned>(l);
        unsigned m = std::max(4U, n);
        std::memcpy(result, "model=\"", 7U);
        result += 7;
        copy_right_justified(result, m, model_id, n, blank);
        result += m;
        std::memcpy(result, "\" ", 2U);
        result += 2;
      }
      if (name != 0) {
        std::memcpy(result, "pdb=\"", 5U);
        result += 5;
      }
      else {
        std::memcpy(result, "pdbres=\"", 8U);
        result += 8;
      }
    }
    if (name != 0) {
      copy_left_justified(result, 4U, name, 4U, blank);
      copy_left_justified(result+4, 1U, altloc, 1U, blank);
      result += 5;
    }
    copy_right_justified(result, 3U, resname, 3U, blank);
    copy_right_justified(result+3, 2U, chain_id, 2U, blank);
    copy_right_justified(result+5, 4U, resseq, 4U, blank);
    copy_left_justified(result+9, 1U, icode, 1U, blank);
    result += 10;
    if (add_model) {
      result[0] = '"';
      result++;
    }
    if (add_segid && segid != 0 && str4(segid).stripped_size() != 0) {
      std::memcpy(result, " segid=\"", 8U);
      copy_left_justified(result+8, 4U, segid, 4U, blank);
      result[12] = '"';
      result += 13;
    }
    if (add_model || add_segid) {
      result[0] = '\0';
    }
  }

  void
  atom_label_columns_formatter::format(
    char* result,
    shared_ptr<chain_data> const& ch_lock,
    bool add_model,
    bool add_segid)
  {
    chain_data const* ch = ch_lock.get();
    if (ch == 0) {
      chain_id = model_id = 0;
      format(result, add_model, add_segid);
    }
    else {
      chain_id = ch->id.elems;
      if (!add_model) {
        model_id = 0;
        format(result, add_model, add_segid);
      }
      else {
        shared_ptr<model_data> md_lock = ch->parent.lock();
        model_data const* md = md_lock.get();
        if (md == 0) {
          model_id = 0;
          format(result, add_model, add_segid);
        }
        else {
          model_id = (md->id.size() == 0 ? 0 : md->id.elems);
          format(result, add_model, add_segid);
        }
      }
    }
  }

  void
  atom_label_columns_formatter::format(
    char* result,
    hierarchy::atom const& atom,
    bool add_model,
    bool add_segid,
    bool pdbres)
  {
    name = (pdbres ? 0 : atom.data->name.elems);
    segid = atom.data->segid.elems;
    shared_ptr<atom_group_data> ag_lock = atom.data->parent.lock();
    const atom_group_data* ag = ag_lock.get();
    if (ag == 0) {
      altloc = resname = resseq = icode = chain_id = model_id = 0;
      format(result, add_model, add_segid);
    }
    else {
      altloc = ag->altloc.elems;
      resname = ag->resname.elems;
      shared_ptr<residue_group_data> rg_lock = ag->parent.lock();
      const residue_group_data* rg = rg_lock.get();
      if (rg == 0) {
        resseq = icode = chain_id = model_id = 0;
        format(result, add_model, add_segid);
      }
      else {
        resseq = rg->resseq.elems;
        icode = rg->icode.elems;
        format(result, rg->parent.lock(), add_model, add_segid);
      }
    }
  }

  void
  atom_label_columns_formatter::format(
    char* result,
    hierarchy::residue const& residue)
  {
    name = 0;
    altloc = 0;
    resname = residue.data->resname.elems;
    resseq = residue.data->resseq.elems;
    icode = residue.data->icode.elems;
    shared_ptr<conformer_data> cf_lock = residue.data->parent.lock();
    const conformer_data* cf = cf_lock.get();
    if (cf == 0) {
      chain_id = model_id = 0;
      format(
        result, /* add_model */ true, /* add_segid */ true);
    }
    else {
      format(
        result, cf->parent.lock(), /* add_model */ true, /* add_segid */ true);
    }
  }

namespace {

  void
  atom_with_labels_init_label_formatter(
    atom_with_labels const& self,
    atom_label_columns_formatter& label_formatter)
  {
    label_formatter.altloc = self.altloc.elems;
    label_formatter.resname = self.resname.elems;
    label_formatter.resseq = self.resseq.elems;
    label_formatter.icode = self.icode.elems;
    label_formatter.chain_id = self.chain_id.elems;
    label_formatter.model_id = (
      self.model_id.size() == 0 ? 0 : self.model_id.elems);
  }

} // namespace <anonymous>

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
    unsigned segid_start,
    unsigned blanks_start_at) const
  {
    char blank = ' ';
    data->segid.copy_left_justified(result+segid_start, 4U, blank);
    unsigned i = segid_start + 4U;
    data->element.copy_right_justified(result+i, 2U, blank);
    i += 2U;
    data->charge.copy_left_justified(result+i, 2U, blank);
    i += 2U;
    if (blanks_start_at != i) {
      for(;;) {
        i--;
        if (result[i] != blank) {
          copy_left_justified(
            result+blanks_start_at, segid_start-blanks_start_at, 0, 0U, blank);
          result[i+1] = '\0';
          return i+1;
        }
        if (i == segid_start) break;
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

  small_str<19>
  atom::pdb_label_columns_segid_small_str() const
  {
    char blank = ' ';
    small_str<19> result(small_str_no_init);
    atom_label_columns_formatter().format(result.elems, *this);
    data->segid.copy_left_justified(result.elems+15, 4U, blank);
    result.elems[19] = '\0';
    return result;
  }

  std::string
  atom::pdb_element_charge_columns() const
  {
    char result[4];
    format_pdb_element_charge_columns(result);
    return std::string(result, 4U);
  }

  std::string
  atom::id_str(bool pdbres, bool suppress_segid) const
  {
    char result[52];
    atom_label_columns_formatter().format(
      result,
      *this,
      /* add_model */ true,
      /* add_segid */ !suppress_segid,
      pdbres);
    return std::string(result);
  }

  std::string
  atom_with_labels::id_str(bool pdbres, bool suppress_segid) const
  {
    char result[52];
    atom_label_columns_formatter label_formatter;
    label_formatter.name = (pdbres ? 0 : data->name.elems);
    label_formatter.segid = data->segid.elems;
    atom_with_labels_init_label_formatter(*this, label_formatter);
    label_formatter.format(
      result, /* add_model */ true, /* add_segid */ !suppress_segid);
    return std::string(result);
  }

  std::string
  residue::id_str(int suppress_segid) const
  {
    char result[50];
    bool throw_segid_not_unique = false;
    atom_label_columns_formatter label_formatter;
    label_formatter.segid = 0;
    if (suppress_segid <= 0) {
      std::vector<atom> const& ats = atoms();
      unsigned n_ats = atoms_size();
      std::set<str4> segid_set;
      for(unsigned i_at=0;i_at<n_ats;i_at++) {
        atom const& a = ats[i_at];
        segid_set.insert(a.data->segid);
      }
      if (segid_set.size() != 0) {
        throw_segid_not_unique = (   segid_set.size() != 1U
                                  && suppress_segid == 0);
        if (throw_segid_not_unique || segid_set.size() == 1U) {
          label_formatter.segid = ats[0].data->segid.elems;
        }
      }
    }
    label_formatter.format(result, *this);
    if (throw_segid_not_unique) {
      throw std::invalid_argument(
        "residue.id_str(suppress_segid=false): segid is not unique:\n"
        + ("  " + std::string(result)));
    }
    return std::string(result);
  }

  unsigned
  atom::format_atom_record(
    char* result,
    atom_label_columns_formatter* label_formatter,
    const char* replace_floats_with) const
  {
    char blank = ' ';
    std::memcpy(result, (data->hetero ? "HETATM" : "ATOM  "), 6U);
    format_atom_record_serial_label_columns(result, label_formatter);
    unsigned segid_start;
    unsigned blanks_start_at;
    if (replace_floats_with != 0) {
      segid_start = 27U;
      unsigned i=0;
      while (replace_floats_with[i] != '\0' && segid_start != 72U) {
        result[segid_start++] = replace_floats_with[i++];
      }
      blanks_start_at = segid_start + 8U;
    }
    else {
      copy_left_justified(result+27, 3U, 0, 0U, blank);
      char *r = result + 30;
      for(unsigned i=0;i<3;i++) {
        std::sprintf(r, "%8.3f", std::min(std::max(-1.e7, data->xyz[i]), 1.e8));
        if (r[8] != '\0' && r[5] != '.' && r[6] != '.' && r[7] != '.') {
          throw std::runtime_error(
            std::string("atom ") + "XYZ"[i] + " coordinate value"
            " does not fit into F8.3 format:\n"
            + "  \"" + std::string(result, 27U) + "\"\n"
            + "  value: " + (boost::format("%.3f") % data->xyz[i]).str());
        }
        r += 8;
      }
      std::sprintf(r, "%6.2f", std::min(std::max(-1.e5, data->occ), 1.e6));
      if (r[6] != '\0' && r[4] != '.' && r[5] != '.') {
        throw std::runtime_error(
          std::string("atom occupancy factor does not fit into F6.2 format:\n")
          + "  \"" + std::string(result, 27U) + "\"\n"
          + "  occupancy factor: " + (boost::format("%.2f") % data->occ).str());
      }
      r += 6;
      std::sprintf(r, "%6.2f", std::min(std::max(-1.e5, data->b), 1.e6));
      if (r[6] != '\0' && r[4] != '.' && r[5] != '.') {
        throw std::runtime_error(
          std::string("atom B-factor does not fit into F6.2 format:\n")
          + "  \"" + std::string(result, 27U) + "\"\n"
          + "  B-factor: " + (boost::format("%.2f") % data->b).str());
      }
      segid_start = 72U;
      blanks_start_at = 66U;
    }
    return format_atom_record_segid_element_charge_columns(
      result, segid_start, blanks_start_at);
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
      std::sprintf(r, "%8.3f", std::min(std::max(-1.e7, data->sigxyz[i]),1.e8));
      if (r[8] != '\0' && r[5] != '.' && r[6] != '.' && r[7] != '.') {
        throw std::runtime_error(
          std::string("atom sigma ") + "XYZ"[i] + " coordinate value"
          " does not fit into F8.3 format:\n"
          + "  \"" + std::string(result, 27U) + "\"\n"
          + "  value: " + (boost::format("%.3f") % data->sigxyz[i]).str());
      }
      r += 8;
    }
    std::sprintf(r, "%6.2f", std::min(std::max(-1.e5, data->sigocc), 1.e6));
    if (r[6] != '\0' && r[4] != '.' && r[5] != '.') {
      throw std::runtime_error(std::string(
          "atom sigma occupancy factor does not fit into F6.2 format:\n")
        + "  \"" + std::string(result, 27U) + "\"\n"
        + "  sigma occupancy factor: "
        + (boost::format("%.2f") % data->sigocc).str());
    }
    r += 6;
    std::sprintf(r, "%6.2f", std::min(std::max(-1.e5, data->sigb), 1.e6));
    if (r[6] != '\0' && r[4] != '.' && r[5] != '.') {
      throw std::runtime_error(std::string(
          "atom sigma B-factor does not fit into F6.2 format:\n")
        + "  \"" + std::string(result, 27U) + "\"\n"
        + "  sigma B-factor: "
        + (boost::format("%.2f") % data->sigb).str());
    }
    return format_atom_record_segid_element_charge_columns(result, 72U, 66U);
  }

namespace {

  void
  throw_f70_error(
    unsigned i,
    double value,
    const char* result,
    const char* sigma)
  {
    static const char* uij_labels[] = {
      "U11", "U22", "U33", "U12", "U13", "U23"};
    throw std::runtime_error(
      std::string("atom ") + sigma + uij_labels[i]
      + " value * 10000 does not fit into F7.0 format:\n"
      + "  \"" + std::string(result, 27U) + "\"\n"
      + "  value * 10000: " + (boost::format("%.0f") % value).str());
  }

} // namespace <anonymous>

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
      double value = data->uij[i]*10000.;
      std::sprintf(r, "%7.0f", std::min(std::max(-1.e7, value), 1.e8));
      r += 7;
      if (*r != '\0') throw_f70_error(i, value, result, "");
    }
    return format_atom_record_segid_element_charge_columns(result, 72U, 70U);
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
      double value =
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
        data->siguij[i]
#else
        -1
#endif
        *10000.;
      std::sprintf(r, "%7.0f", std::min(std::max(-1.e7, value), 1.e8));
      r += 7;
      if (*r != '\0') throw_f70_error(i, value, result, "sigma ");
    }
    return format_atom_record_segid_element_charge_columns(result, 72U, 70U);
  }

  std::string
  atom_with_labels::format_atom_record(
    const char* replace_floats_with) const
  {
    char result[81];
    atom_label_columns_formatter label_formatter;
    atom_with_labels_init_label_formatter(*this, label_formatter);
    unsigned str_len = atom::format_atom_record(
      result, &label_formatter, replace_floats_with);
    return std::string(result, str_len);
  }

#define IOTBX_LOC(type) \
  std::string \
  atom_with_labels::format_##type##_record() const \
  { \
    char result[81]; \
    atom_label_columns_formatter label_formatter; \
    atom_with_labels_init_label_formatter(*this, label_formatter); \
    unsigned str_len = atom::format_##type##_record( \
      result, &label_formatter); \
    return std::string(result, str_len); \
  }

  IOTBX_LOC(sigatm)
  IOTBX_LOC(anisou)
  IOTBX_LOC(siguij)

#undef IOTBX_LOC

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

  std::string
  atom_with_labels::format_atom_record_group(
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij) const
  {
    char result[324];
    atom_label_columns_formatter label_formatter;
    atom_with_labels_init_label_formatter(*this, label_formatter);
    unsigned str_len = atom::format_atom_record_group(
      result, &label_formatter, atom_hetatm, sigatm, anisou, siguij);
    return std::string(result, str_len);
  }

  std::string
  atom::quote(bool full) const
  {
    char result[82];
    result[0] = '"';
    unsigned str_len = format_atom_record(
      result+1,
      /* label_formatter */ 0,
      /* replace_floats_with */ (full ? 0 : ".*."));
    result[++str_len] = '"';
    return std::string(result, ++str_len);
  }

  std::string
  atom_with_labels::quote(bool full) const
  {
    char result[82];
    atom_label_columns_formatter label_formatter;
    atom_with_labels_init_label_formatter(*this, label_formatter);
    result[0] = '"';
    unsigned str_len = atom::format_atom_record(
      result+1,
      &label_formatter,
      /* replace_floats_with */ (full ? 0 : ".*."));
    result[++str_len] = '"';
    return std::string(result, ++str_len);
  }

  atom_with_labels
  atom::fetch_labels() const
  {
    str8 model_id;
    str2 chain_id;
    str4 resseq;
    str1 icode;
    str1 altloc;
    str3 resname;
    boost::optional<atom_group> ag = parent();
    if (ag) {
      altloc = ag->data->altloc;
      resname = ag->data->resname;
      boost::optional<residue_group> rg = ag->parent();
      if (rg) {
        resseq = rg->data->resseq;
        icode = rg->data->icode;
        boost::optional<chain> ch = rg->parent();
        if (ch) {
          chain_id = ch->data->id;
          boost::optional<model> md = ch->parent();
          if (md) {
            model_id = md->data->id;
          }
        }
      }
    }
    return atom_with_labels(
      *this,
      model_id.elems,
      chain_id.elems,
      resseq.elems,
      icode.elems,
      altloc.elems,
      resname.elems,
      /* is_first_in_chain */ false,
      /* is_first_after_break */ false);
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
    char pad_with = ' ';
    char erj[2];
    data->element.copy_right_justified(erj, 2U, pad_with);
    std::string e(erj, 2);
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
    if (e == "  " && data->name.size() >= 2) {
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

  bool
  atom::set_chemical_element_simple_if_necessary(
    bool tidy_existing)
  {
    if (data->element.stripped_size() != 0 && !tidy_existing) return false;
    boost::optional<std::string> e = determine_chemical_element_simple();
    if (!e || *e == data->element.elems) return false;
    IOTBX_ASSERT(e->size() <= 2);
    char pad_with = ' ';
    copy_right_justified(
      data->element.elems,
      data->element.capacity(),
      e->c_str(),
      static_cast<unsigned>(e->size()),
      pad_with);
    data->element.elems[data->element.capacity()] = '\0';
    return true;
  }

  boost::optional<std::string>
  atom::charge_tidy(
    bool strip) const
  {
    char charge[2];
    std::memcpy(charge, data->charge.elems, 2);
    unsigned result_size = 2;
    while (true) {
      if (charge[0] == '\0') {
        if (strip) {
          result_size = 0;
        }
        else {
          charge[0] = charge[1] = ' ';
        }
        break;
      }
      if (charge[1] == '\0') charge[1] = ' ';
      if (   (charge[0] == ' ' || charge[0] == '0')
          && (charge[1] == ' ' || charge[1] == '0')) {
        if (strip) {
          result_size = 0;
        }
        else {
          charge[0] = charge[1] = ' ';
        }
        break;
      }
      if (charge[0] == '+' || charge[0] == '-') {
        if (charge[1] == charge[0]) {
          charge[0] = '2';
          break;
        }
        std::swap(charge[0], charge[1]);
      }
      if (charge[1] == '+' || charge[1] == '-') {
        if (charge[0] == ' ') {
          charge[0] = '1';
          break;
        }
        if (isdigit(charge[0])) {
          break;
        }
      }
      return boost::optional<std::string>();
    }
    return boost::optional<std::string>(std::string(charge, result_size));
  }

  void
  atom_group::append_atom_with_other_parent(hierarchy::atom const& atom)
  {
    data->atoms.push_back(atom);
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

  str4
  atom_group::confid_small_str() const
  {
    char blank = ' ';
    str4 result;
    data->altloc.copy_left_justified(result.elems, 1U, blank);
    data->resname.copy_right_justified(result.elems+1, 3U, blank);
    result.elems[4] = '\0';
    return result;
  }

namespace {

  std::string
  make_resid(str4 const& resseq, str1 const& icode)
  {
    char blank = ' ';
    char result[6];
    resseq.copy_right_justified(result, 4U, blank);
    icode.copy_left_justified(result+4, 1U, blank);
    result[5] = '\0';
    return std::string(result);
  }

} // namespace <anonymous>

  std::string
  residue_group::resid() const
  {
    return make_resid(data->resseq, data->icode);
  }

  str5
  residue_group::resid_small_str() const
  {
    char blank = ' ';
    str5 result(small_str_no_init);
    data->resseq.copy_right_justified(result.elems, 4U, blank);
    data->icode.copy_left_justified(result.elems+4, 1U, blank);
    result.elems[5] = '\0';
    return result;
  }

  std::string
  residue::resid() const
  {
    return make_resid(data->resseq, data->icode);
  }

  boost::optional<atom>
  residue::find_atom_by(char const* name) const
  {
    boost::optional<atom> result;
    if (name != 0) {
      std::vector<atom> const& ats = atoms();
      unsigned n_ats = atoms_size();
      for(unsigned i_at=0;i_at<n_ats;i_at++) {
        atom const& a = ats[i_at];
        if (std::strcmp(a.data->name.elems, name) == 0) {
          return boost::optional<atom>(a);
        }
      }
    }
    return boost::optional<atom>();
  }

  std::string
  atom_with_labels::resid() const
  {
    return make_resid(resseq, icode);
  }

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
    IOTBX_ASSERT(secondary.data->altloc == primary.data->altloc);
    IOTBX_ASSERT(secondary.data->resname == primary.data->resname);
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
    IOTBX_ASSERT(secondary.atoms_size() == 0);
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
          hierarchy::atom atom = ag.atoms()[i_atom];
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
    IOTBX_ASSERT(secondary.data->resseq == primary.data->resseq);
    IOTBX_ASSERT(secondary.data->icode == primary.data->icode);
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
    IOTBX_ASSERT(secondary.atom_groups_size() == 0);
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
            IOTBX_ASSERT(secondary.atom_groups_size() == 0);
            IOTBX_ASSERT(secondary.parent_ptr().get() == 0);
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
      IOTBX_ASSERT(result.size() == result_size);
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

  residue::residue(
    hierarchy::root const& root)
  :
    root_(root)
  {
    std::vector<model> const& models = root.models();
    IOTBX_ASSERT(models.size() == 1);
    std::vector<chain> const& chains = models[0].chains();
    IOTBX_ASSERT(chains.size() == 1);
    af::shared<conformer> conformers = chains[0].conformers();
    IOTBX_ASSERT(conformers.size() == 1);
    std::vector<residue> const& residues = conformers[0].residues();
    IOTBX_ASSERT(residues.size() == 1);
    data = residues[0].data;
  }

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
  conformer::build_from_residue_groups(
    const hierarchy::chain* chain,
    const residue_group* residue_groups,
    unsigned residue_groups_size)
  {
    const char nulc = '\0';
    std::vector<char> altlocs;
    typedef std::map<char, unsigned> mcu;
    mcu altloc_indices;
    bool have_at_least_one_atom_group = false;
    for(unsigned i_rg=0;i_rg<residue_groups_size;i_rg++) {
      residue_group const& rg = residue_groups[i_rg];
      unsigned n_ag = rg.atom_groups_size();
      if (n_ag != 0) have_at_least_one_atom_group = true;
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
    if (!have_at_least_one_atom_group) return af::shared<conformer>();
    unsigned n_cf = static_cast<unsigned>(altloc_indices.size());
    if (n_cf == 0) {
      altlocs.push_back(nulc);
      n_cf++;
    }
    af::shared<conformer> result((af::reserve(n_cf)));
    for(unsigned i_cf=0;i_cf<n_cf;i_cf++) {
      if (chain != 0) {
        result.push_back(conformer(*chain, str1(altlocs[i_cf]).elems));
      }
      else {
        result.push_back(conformer(str1(altlocs[i_cf]).elems));
      }
    }
    std::vector<str3> resnames; // allocate once
    resnames.reserve(32U); // not critical
    std::set<str3> resnames_with_altloc; // allocate once
    std::vector<std::vector<atom_group> > altloc_ags(n_cf); // allocate once
    for(unsigned i_rg=0;i_rg<residue_groups_size;i_rg++) {
      residue_group const& rg = residue_groups[i_rg];
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
          IOTBX_ASSERT(rai != resname_atoms.end());
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

  af::shared<conformer>
  chain::conformers() const
  {
    unsigned n_rg = residue_groups_size();
    if (n_rg == 0) return af::shared<conformer>();
    return conformer::build_from_residue_groups(
      this, &*residue_groups().begin(), n_rg);
  }

  af::shared<conformer>
  residue_group::conformers() const
  {
    chain ch;
    const chain* ch_ptr = 0;
    shared_ptr<chain_data> p = data->parent.lock();
    if (p.get() != 0) {
      ch = chain(p, true);
      ch_ptr = &ch;
    }
    return conformer::build_from_residue_groups(ch_ptr, this, 1U);
  }

  bool
  atom_group::is_identical_hierarchy(
    atom_group const& other) const
  {
    std::vector<atom> const& ats = data->atoms;
    std::vector<atom> const& oats = other.data->atoms;
    unsigned n = static_cast<unsigned>(ats.size());
    if (n != static_cast<unsigned>(oats.size())) return false;
    for(unsigned i=0;i<n;i++) {
      atom_data const& a = *ats[i].data;
      atom_data const& oa = *oats[i].data;
      if (a.name != oa.name) return false;
      if (a.element != oa.element) return false;
      if (a.charge != oa.charge) return false;
      if (a.serial != oa.serial) return false;
      if (a.hetero != oa.hetero) return false;
    }
    return true;
  }

  bool
  residue_group::is_identical_hierarchy(
    residue_group const& other) const
  {
    std::vector<atom_group> const& ags = data->atom_groups;
    std::vector<atom_group> const& oags = other.data->atom_groups;
    unsigned n = static_cast<unsigned>(ags.size());
    if (n != static_cast<unsigned>(oags.size())) return false;
    for(unsigned i=0;i<n;i++) {
      if (ags[i].data->altloc != oags[i].data->altloc) return false;
      if (ags[i].data->resname != oags[i].data->resname) return false;
      if (!ags[i].is_identical_hierarchy(oags[i])) return false;
    }
    return true;
  }

namespace {

  void
  get_confid_atom_names(
    std::map<str4, std::vector<str4> >& result,
    residue_group const& rg)
  {
    unsigned n_ag = rg.atom_groups_size();
    std::vector<atom_group> const& ags = rg.atom_groups();
    for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
      atom_group const& ag = ags[i_ag];
      std::vector<str4>& atom_names = result[ag.confid_small_str()];
      unsigned n_at = ag.atoms_size();
      std::vector<atom> const& ats = ag.atoms();
      for(unsigned i_at=0;i_at<n_at;i_at++) {
        atom_names.push_back(ats[i_at].data->name);
      }
    }
    typedef std::map<str4, std::vector<str4> >::iterator it;
    it i_end = result.end();
    for(it i=result.begin();i!=i_end;i++) {
      std::sort(i->second.begin(), i->second.end());
    }
  }

} // namespace <anonymous>

  bool
  residue_group::is_similar_hierarchy(
    residue_group const& other) const
  {
    std::map<str4, std::vector<str4> > confid_atom_names[2];
    get_confid_atom_names(confid_atom_names[0], *this);
    get_confid_atom_names(confid_atom_names[1], other);
    typedef std::map<str4, std::vector<str4> >::const_iterator it;
    if (   confid_atom_names[0].size()
        != confid_atom_names[1].size()) return false;
    it i_end = confid_atom_names[0].end();
    it j_end = confid_atom_names[1].end();
    for(it i=confid_atom_names[0].begin();i!=i_end;i++) {
      it j = confid_atom_names[1].find(i->first);
      if (j == j_end) return false;
      if (j->second != i->second) return false;
    }
    return true;
  }

namespace {

  bool
  chain_equivalence(
    chain const& self,
    chain const& other,
    int equivalence_type)
  {
    std::vector<residue_group> const& rgs = self.residue_groups();
    std::vector<residue_group> const& orgs = other.residue_groups();
    unsigned n = static_cast<unsigned>(rgs.size());
    if (n != static_cast<unsigned>(orgs.size())) return false;
    for(unsigned i=0;i<n;i++) {
      if (rgs[i].data->resseq != orgs[i].data->resseq) return false;
      if (rgs[i].data->icode != orgs[i].data->icode) return false;
      if (equivalence_type == 0) {
        if (!rgs[i].is_identical_hierarchy(orgs[i])) return false;
      }
      else {
        if (!rgs[i].is_similar_hierarchy(orgs[i])) return false;
      }
    }
    return true;
  }

  bool
  model_equivalence(
    model const& self,
    model const& other,
    int equivalence_type)
  {
    std::vector<chain> const& chs = self.chains();
    std::vector<chain> const& ochs = other.chains();
    unsigned n = static_cast<unsigned>(chs.size());
    if (n != static_cast<unsigned>(ochs.size())) return false;
    for(unsigned i=0;i<n;i++) {
      if (chs[i].data->id != ochs[i].data->id) return false;
      if (equivalence_type == 0) {
        if (!chs[i].is_identical_hierarchy(ochs[i])) return false;
      }
      else {
        if (!chs[i].is_similar_hierarchy(ochs[i])) return false;
      }
    }
    return true;
  }

} // namespace <anonymous>

  bool
  chain::is_identical_hierarchy(
    chain const& other) const { return chain_equivalence(*this, other, 0); }

  bool
  chain::is_similar_hierarchy(
    chain const& other) const { return chain_equivalence(*this, other, 1); }

  bool
  model::is_identical_hierarchy(
    model const& other) const { return model_equivalence(*this, other, 0); }

  bool
  model::is_similar_hierarchy(
    model const& other) const { return model_equivalence(*this, other, 1); }

  bool
  root::is_similar_hierarchy(
    root const& other) const
  {
    std::vector<model> const& mds = models();
    std::vector<model> const& omds = other.models();
    unsigned n = static_cast<unsigned>(mds.size());
    if (n != static_cast<unsigned>(omds.size())) return false;
    for(unsigned i=0;i<n;i++) {
      if (mds[i].data->id != omds[i].data->id) return false;
      if (!mds[i].is_similar_hierarchy(omds[i])) return false;
    }
    return true;
  }

  void
  model::transfer_chains_from_other(
    model& other)
  {
    std::vector<chain>& other_chains = other.data->chains;
    unsigned n = other.chains_size();
    pre_allocate_chains(n);
    for(unsigned i=0;i<n;i++) {
      other_chains[i].clear_parent();
      append_chain(other_chains[i]);
    }
    std::vector<chain> empty;
    other_chains.swap(empty);
  }

namespace {

  struct interleaved_conf_helper
  {
    const atom* atom_ptr;
    unsigned resname_index;
    unsigned atom_name_index;
    unsigned sequential_index;

    interleaved_conf_helper(
      const atom* atom_ptr_,
      unsigned resname_index_,
      unsigned atom_name_index_,
      unsigned sequential_index_)
    :
      atom_ptr(atom_ptr_),
      resname_index(resname_index_),
      atom_name_index(atom_name_index_),
      sequential_index(sequential_index_)
    {}

    bool
    operator<(interleaved_conf_helper const& other) const
    {
      if (resname_index < other.resname_index) return true;
      if (resname_index > other.resname_index) return false;
      if (atom_name_index < other.atom_name_index) return true;
      if (atom_name_index > other.atom_name_index) return false;
      return (sequential_index < other.sequential_index);
    }
  };

} // namespace <anonymous>

  void
  residue_group::atoms_interleaved_conf_impl(
    bool group_residue_names,
    af::shared<atom>& result) const
  {
    std::map<str3, unsigned> resname_indices;
    unsigned resname_indices_size = 0;
    unsigned resname_index = 0;
    std::vector<std::map<str4, unsigned> > resname_atom_name_indices;
    std::map<str4, unsigned>* atom_name_indices = 0;
    unsigned atom_name_indices_size = 0;
    if (!group_residue_names) {
      resname_atom_name_indices.resize(1U);
      atom_name_indices = &resname_atom_name_indices[0];
    }
    typedef std::vector<interleaved_conf_helper> tab_t;
    tab_t tab;
    tab.reserve(100U); // not critical
    unsigned sequential_index = 0;
    unsigned n_ag = atom_groups_size();
    std::vector<atom_group> const& ags = atom_groups();
    for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
      atom_group const& ag = ags[i_ag];
      if (group_residue_names) {
        resname_index = resname_indices[ag.data->resname];
        if (resname_index == 0) {
          resname_indices[ag.data->resname]
            = resname_index
            = ++resname_indices_size;
          resname_atom_name_indices.resize(resname_indices_size);
        }
        atom_name_indices = &resname_atom_name_indices[resname_index-1U];
        atom_name_indices_size = static_cast<unsigned>(
          atom_name_indices->size());
      }
      unsigned n_at = ag.atoms_size();
      std::vector<atom> const& ats = ag.atoms();
      for(unsigned i_at=0;i_at<n_at;i_at++) {
        atom const& a = ats[i_at];
        unsigned atom_name_index = (*atom_name_indices)[a.data->name];
        if (atom_name_index == 0) {
          (*atom_name_indices)[a.data->name]
            = atom_name_index
            = ++atom_name_indices_size;
        }
        tab.push_back(interleaved_conf_helper(
          &a,
          resname_index,
          atom_name_index,
          sequential_index++));
      }
    }
    std::sort(tab.begin(), tab.end());
    typedef tab_t::const_iterator it;
    it i_end = tab.end();
    for(it i=tab.begin();i!=i_end;i++) {
      result.push_back(*(i->atom_ptr));
    }
  }

  af::shared<atom>
  residue_group::atoms_interleaved_conf(
    bool group_residue_names) const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    atoms_interleaved_conf_impl(group_residue_names, result);
    return result;
  }

  af::shared<std::string>
  residue_group::unique_resnames() const
  {
    unsigned n_ag = atom_groups_size();
    std::vector<atom_group> const& ags = atom_groups();
    typedef std::set<str3> ss3;
    ss3 resname_set;
    for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
      atom_group const& ag = ags[i_ag];
      resname_set.insert(ag.data->resname);
    }
    af::shared<std::string> result((af::reserve(resname_set.size())));
    typedef ss3::const_iterator it;
    for(it i=resname_set.begin();i!=resname_set.end();i++) {
      result.push_back(std::string(i->elems));
    }
    return result;
  }

  af::shared<atom>
  chain::atoms_interleaved_conf(
    bool group_residue_names) const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_CHAIN_RESIDUE_GROUPS_LOOPS
      rgs[i_rg].atoms_interleaved_conf_impl(group_residue_names, result);
    }
    return result;
  }

  af::shared<atom>
  model::atoms_interleaved_conf(
    bool group_residue_names) const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_MODEL_RESIDUE_GROUPS_LOOPS
      rgs[i_rg].atoms_interleaved_conf_impl(group_residue_names, result);
    }}
    return result;
  }

  af::shared<atom>
  root::atoms_interleaved_conf(
    bool group_residue_names) const
  {
    af::shared<atom> result((af::reserve(atoms_size())));
    IOTBX_PDB_HIERARCHY_CPP_ROOT_RESIDUE_GROUPS_LOOPS
      rgs[i_rg].atoms_interleaved_conf_impl(group_residue_names, result);
    }}}
    return result;
  }

#define IOTBX_LOC(T) \
  af::shared<atom> \
  T::atoms( \
    int interleaved_conf) const \
  { \
    if (interleaved_conf <= 0) return atoms_sequential_conf(); \
    return atoms_interleaved_conf((interleaved_conf < 2)); \
  }

  IOTBX_LOC(residue_group)
  IOTBX_LOC(chain)
  IOTBX_LOC(model)
  IOTBX_LOC(root)

#undef IOTBX_LOC

  af::shared<atom>
  root::atoms_with_i_seq_mismatch() const
  {
    af::shared<atom> result;
    unsigned i_seq = 0;
    IOTBX_PDB_HIERARCHY_CPP_ROOT_ATOM_GROUPS_LOOPS
      std::vector<atom> const& ats = ags[i_ag].atoms();
      for(unsigned i_at=0;i_at<n_ats;i_at++) {
        if (ats[i_at].data->i_seq != i_seq++) {
          result.push_back(ats[i_at]);
        }
      }
    }}}}
    return result;
  }

  atom_with_labels::atom_with_labels()
  :
    is_first_in_chain(false),
    is_first_after_break(false)
  {}

  atom_with_labels::atom_with_labels(
    atom const& atom_,
    const char* model_id_,
    const char* chain_id_,
    const char* resseq_,
    const char* icode_,
    const char* altloc_,
    const char* resname_,
    bool is_first_in_chain_,
    bool is_first_after_break_)
  :
    atom(atom_),
    model_id(model_id_),
    chain_id(chain_id_),
    resseq(resseq_),
    icode(icode_),
    altloc(altloc_),
    resname(resname_),
    is_first_in_chain(is_first_in_chain_),
    is_first_after_break(is_first_after_break_)
  {}

  atom_with_labels
  atom_with_labels::detached_copy() const
  {
    return atom_with_labels(
      atom::detached_copy(),
      model_id.elems,
      chain_id.elems,
      resseq.elems,
      icode.elems,
      altloc.elems,
      resname.elems,
      is_first_in_chain,
      is_first_after_break);
  }

  int
  atom::serial_as_int() const
  {
    str5 const& s = data->serial;
    int result = -1;
    const char* errmsg = hy36decode(5U, s.elems, s.size(), &result);
    if (errmsg) {
      throw std::invalid_argument(
        "invalid atom serial number: \""
        + std::string(s.elems, s.size()) + "\"");
    }
    return result;
  }

  int
  atom_with_labels::serial_as_int() const
  {
    str5 const& s = data->serial;
    int result = -1;
    const char* errmsg = hy36decode(5U, s.elems, s.size(), &result);
    if (errmsg) {
      throw std::invalid_argument(
        "invalid atom serial number:\n"
        "  " + format_atom_record() + "\n"
        "        ^^^^^");
    }
    return result;
  }

namespace {

  std::string
  invalid_resseq_message(atom_with_labels const& awl)
  {
    return
      "invalid residue sequence number:\n"
      "  " + awl.format_atom_record() + "\n"
      "                        ^^^^";
  }

  std::string
  invalid_resseq_message(
    str4 const& resseq,
    std::size_t atoms_size,
    const atom* atoms)
  {
    if (atoms_size != 0) {
      try {
        atom_with_labels awl = atoms[0].fetch_labels();
        awl.resseq = resseq;
        return invalid_resseq_message(awl);
      }
      catch (std::exception const&) {}
    }
    return
      "invalid residue sequence number: \""
      + std::string(resseq.elems, resseq.size()) + "\"";
  }

} // namespace <anonymous>

  int
  atom_with_labels::resseq_as_int() const
  {
    str4 const& s = resseq;
    int result = -1;
    const char* errmsg = hy36decode(4U, s.elems, s.size(), &result);
    if (errmsg) {
      throw std::invalid_argument(invalid_resseq_message(*this));
    }
    return result;
  }

  int
  residue_group::resseq_as_int() const
  {
    str4 const& s = data->resseq;
    int result = -1;
    const char* errmsg = hy36decode(4U, s.elems, s.size(), &result);
    if (errmsg) {
      af::shared<atom> ats = atoms_sequential_conf();
      throw std::invalid_argument(
        invalid_resseq_message(s, ats.size(), ats.begin()));
    }
    return result;
  }

  int
  residue::resseq_as_int() const
  {
    str4 const& s = data->resseq;
    int result = -1;
    const char* errmsg = hy36decode(4U, s.elems, s.size(), &result);
    if (errmsg) {
      std::size_t n = data->atoms.size();
      throw std::invalid_argument(
        invalid_resseq_message(s, n, (n == 0 ? 0 : &*data->atoms.begin())));
    }
    return result;
  }

}}} // namespace iotbx::pdb::hierarchy
