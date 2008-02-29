#include <iotbx/pdb/hierarchy_v2.h>
#include <cctbx/eltbx/chemical_elements.h>

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

#define IOTBX_PDB_HIERARCHY_V2_SET_PARENT_ETC(P, T) \
  T& \
  T::set_parent(P const& parent) \
  { \
    if (data->parent.lock().get() != 0) { \
      throw std::runtime_error(#T " has another parent " #P " already."); \
    } \
    data->parent = parent.data; \
    return *this; \
  } \
\
  boost::optional<P> \
  T::parent() const \
  { \
    shared_ptr<P##_data> parent = data->parent.lock(); \
    if (parent.get() == 0) return boost::optional<P>(); \
    return boost::optional<P>(P(parent, true)); \
  }

#define IOTBX_PDB_HIERARCHY_V2_APPEND_ETC(T, C) \
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
  T::new_##C##s(unsigned number_of_additional_##C##s) \
  { \
    pre_allocate_##C##s(number_of_additional_##C##s); \
    for(unsigned i=0;i<number_of_additional_##C##s;i++) { \
      data->C##s.push_back(C(*this)); \
    } \
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
  } \
\
  unsigned \
  T::reset_atom_tmp(int new_value) const \
  { \
    unsigned n = C##s_size(); \
    if (n == 0) return 0; \
    unsigned result = 0; \
    const C* c = &*data->C##s.begin(); \
    for(unsigned i=0;i<n;i++) { \
      result += (c++)->reset_atom_tmp(new_value); \
    } \
    return result; \
  }

#define IOTBX_PDB_HIERARCHY_V2_DETACHED_COPY_ETC(P, T, C) \
  IOTBX_PDB_HIERARCHY_V2_SET_PARENT_ETC(P, T) \
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
    T result; \
    detach_copy_children(result, result.data->C##s, data->C##s); \
    return result; \
  } \
  IOTBX_PDB_HIERARCHY_V2_APPEND_ETC(T, C)

  IOTBX_PDB_HIERARCHY_V2_APPEND_ETC(root, model)
  IOTBX_PDB_HIERARCHY_V2_DETACHED_COPY_ETC(root, model, chain)
  IOTBX_PDB_HIERARCHY_V2_DETACHED_COPY_ETC(model, chain, residue_group)
  IOTBX_PDB_HIERARCHY_V2_DETACHED_COPY_ETC(chain, residue_group, atom_group)
  IOTBX_PDB_HIERARCHY_V2_DETACHED_COPY_ETC(residue_group, atom_group, atom)
  IOTBX_PDB_HIERARCHY_V2_SET_PARENT_ETC(atom_group, atom)

  root
  root::deep_copy() const
  {
    root result;
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
  atom::format_atom_record(
    char* result,
    const char* altloc,
    const char* resname,
    const char* resseq,
    const char* icode,
    const char* chain_id) const
  {
    char blank = ' ';
    atom_data const& d = *data;
    std::memcpy(result, (d.hetero ? "HETATM" : "ATOM  "), 6U);
    d.serial.copy_padded(result+6, 5U, blank);
    result[11] = blank;
    d.name.copy_padded(result+12, 4U, blank);
    copy_padded(result+16, 1U, altloc, 1U, blank);
    copy_padded(result+17, 3U, resname, 3U, blank);
    if (chain_id == 0 || chain_id[0] == '\0') {
      result[20] = blank;
      result[21] = blank;
    }
    else if (chain_id[1] == '\0') {
      result[20] = blank;
      result[21] = chain_id[0];
    }
    else {
      result[20] = chain_id[0];
      result[21] = chain_id[1];
    }
    copy_padded(result+22, 4U, resseq, 4U, blank);
    copy_padded(result+26, 4U, icode, 1U, blank);
    char *r = result + 30;
    for(unsigned i=0;i<3;i++) {
      std::sprintf(r, "%8.3f", d.xyz[i]);
      r += 8;
      if (*r != '\0') {
        throw std::runtime_error(
          std::string("atom ") + "XYZ"[i] + " coordinate value"
          " does not fit into F8.3 format.");
      }
    }
    std::sprintf(r, "%6.2f", d.occ);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom ") + "occupancy factor"
        " does not fit into F6.2 format.");
    }
    std::sprintf(r, "%6.2f", d.b);
    r += 6;
    if (*r != '\0') {
      throw std::runtime_error(
        std::string("atom ") + "B-factor"
        " does not fit into F6.2 format.");
    }
    copy_padded(r, 6U, 0, 0U, blank);
    d.segid.copy_padded(result+72, 4U, blank);
    d.element.copy_padded(result+76, 2U, blank);
    d.charge.copy_padded(result+78, 2U, blank);
    for(unsigned i=79;i!=71;i--) {
      if (result[i] != blank) {
        result[i+1] = '\0';
        return i+1;
      }
    }
    result[66] = '\0';
    return 66U;
  }

  unsigned
  atom::format_atom_record_using_parents(
    char* result) const
  {
    shared_ptr<atom_group_data> ag_lock = data->parent.lock();
    const atom_group_data* ag = ag_lock.get();
    if (ag == 0) {
      return format_atom_record(result,
        0, 0,
        0, 0,
        0);
    }
    shared_ptr<residue_group_data> rg_lock = ag->parent.lock();
    const residue_group_data* rg = rg_lock.get();
    if (rg == 0) {
      return format_atom_record(result,
        ag->altloc.elems, ag->resname.elems,
        0, 0,
        0);
    }
    shared_ptr<chain_data> c_lock = rg->parent.lock();
    chain_data const* c = c_lock.get();
    return format_atom_record(result,
      ag->altloc.elems, ag->resname.elems,
      rg->resseq.elems, rg->icode.elems,
      (c == 0 ? 0 : c->id.c_str()));
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

}}} // namespace iotbx::pdb::hierarchy_v2
