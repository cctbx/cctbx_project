#include <iotbx/pdb/hierarchy_v2_atoms.h>
#include <iotbx/pdb/hybrid_36_c.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 { namespace atoms {

  af::shared<std::string>
  extract_serial(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy_v2::atom* atoms_end = atoms.end(); \
    for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++) {
      result.push_back(a->data->serial.elems);
    }
    return result;
  }

  af::shared<std::string>
  extract_name(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy_v2::atom* atoms_end = atoms.end(); \
    for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++) {
      result.push_back(a->data->name.elems);
    }
    return result;
  }

#define IOTBX_LOC(attr, attr_type) \
  af::shared<attr_type > \
  extract_##attr( \
    af::const_ref<atom> const& atoms) \
  { \
    af::shared<attr_type > result( \
      atoms.size(), af::init_functor_null<attr_type >()); \
    attr_type* r = result.begin(); \
    const hierarchy_v2::atom* atoms_end = atoms.end(); \
    for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++) { \
      *r++ = a->data->attr; \
    } \
    return result; \
  }

  IOTBX_LOC(xyz, vec3)
  IOTBX_LOC(sigxyz, vec3)
  IOTBX_LOC(occ, double)
  IOTBX_LOC(sigocc, double)
  IOTBX_LOC(b, double)
  IOTBX_LOC(sigb, double)
  IOTBX_LOC(uij, sym_mat3)
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
  IOTBX_LOC(siguij, sym_mat3)
#endif
  IOTBX_LOC(i_seq, std::size_t)

#undef IOTBX_LOC

  af::shared<std::size_t>
  extract_hetero(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::size_t> result;
    const hierarchy_v2::atom* atoms_end = atoms.end();
    std::size_t i_seq = 0;
    for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++,i_seq++) {
      if (a->data->hetero) result.push_back(i_seq);
    }
    return result;
  }

  af::shared<std::string>
  extract_element(
    af::const_ref<atom> const& atoms,
    bool strip)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy_v2::atom* atoms_end = atoms.end();
    if (strip) {
      for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++) {
        result.push_back(a->data->element.strip().elems);
      }
    }
    else {
      for(const hierarchy_v2::atom* a=atoms.begin();a!=atoms_end;a++) {
        result.push_back(a->data->element.elems);
      }
    }
    return result;
  }

#define IOTBX_LOC(attr, attr_type) \
  void \
  set_##attr( \
    af::ref<atom> const& atoms, \
    af::const_ref<attr_type > const& new_##attr) \
  { \
    SCITBX_ASSERT(new_##attr.size() == atoms.size()); \
    unsigned n = static_cast<unsigned>(atoms.size()); \
    for(unsigned i=0;i<n;i++) { \
      atoms[i].data->attr = new_##attr[i]; \
    } \
  }

  IOTBX_LOC(xyz, vec3)
  IOTBX_LOC(sigxyz, vec3)
  IOTBX_LOC(occ, double)
  IOTBX_LOC(sigocc, double)
  IOTBX_LOC(b, double)
  IOTBX_LOC(sigb, double)
  IOTBX_LOC(uij, sym_mat3)
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
  IOTBX_LOC(siguij, sym_mat3)
#endif

#undef IOTBX_LOC

  void
  reset_serial(
    af::const_ref<atom> const& atoms,
    int first_value)
  {
    int value = first_value;
    for(const atom* a=atoms.begin();a!=atoms.end();a++) {
      const char* errmsg = hy36encode(5U, value++, a->data->serial.elems);
      if (errmsg != 0) {
        if (std::strcmp(errmsg, "value out of range.")) {
          errmsg = "PDB atom serial number out of range.";
        }
        throw std::runtime_error(errmsg);
      }
    }
  }

  void
  reset_i_seq(
    af::const_ref<atom> const& atoms)
  {
    unsigned value = 0;
    for(const atom* a=atoms.begin();a!=atoms.end();a++) {
      a->data->i_seq = value++;
    }
  }

  void
  reset_tmp(
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
  reset_tmp_for_occupancy_groups_simple(
    af::const_ref<atom> const& atoms)
  {
    int value = 0;
    for(const atom* a=atoms.begin();a!=atoms.end();a++,value++) {
      a->data->tmp = (a->element_is_hydrogen() ? -1 : value);
    }
  }

}}}} // namespace iotbx::pdb::hierarchy_v2::atoms
