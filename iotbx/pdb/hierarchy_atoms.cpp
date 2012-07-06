#include <iotbx/pdb/hierarchy_atoms.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <iotbx/error.h>
#include <memory>

namespace iotbx { namespace pdb { namespace hierarchy { namespace atoms {

  af::shared<std::string>
  extract_serial(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy::atom* atoms_end = atoms.end(); \
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
      result.push_back(a->data->serial.elems);
    }
    return result;
  }

  af::shared<std::string>
  extract_name(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy::atom* atoms_end = atoms.end(); \
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
      result.push_back(a->data->name.elems);
    }
    return result;
  }

  af::shared<std::string>
  extract_segid(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::string> result((af::reserve(atoms.size())));
    const hierarchy::atom* atoms_end = atoms.end(); \
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
      result.push_back(a->data->segid.elems);
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
    const hierarchy::atom* atoms_end = atoms.end(); \
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) { \
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
  extract_tmp_as_size_t(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::size_t> result(
      atoms.size(), af::init_functor_null<std::size_t>());
    std::size_t* r = result.begin();
    const hierarchy::atom* atoms_end = atoms.end();
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
      int tmp = a->data->tmp;
      if (tmp < 0) {
        throw std::runtime_error(
          "atom.tmp less than zero: cannot convert to unsigned value.");
      }
      *r++ = static_cast<std::size_t>(tmp);
    }
    return result;
  }

  af::shared<std::size_t>
  extract_hetero(
    af::const_ref<atom> const& atoms)
  {
    af::shared<std::size_t> result;
    const hierarchy::atom* atoms_end = atoms.end();
    std::size_t i_seq = 0;
    for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++,i_seq++) {
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
    const hierarchy::atom* atoms_end = atoms.end();
    if (strip) {
      for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
        result.push_back(a->data->element.strip().elems);
      }
    }
    else {
      for(const hierarchy::atom* a=atoms.begin();a!=atoms_end;a++) {
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
    IOTBX_ASSERT(new_##attr.size() == atoms.size()); \
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

} // namespace atoms

  void
  root::atoms_reset_serial(
    int interleaved_conf,
    int first_value) const
  {
    std::vector<model> const& models = this->models();
    unsigned n_mds = models_size();
    for(unsigned i_md=0;i_md<n_mds;i_md++) {
      atoms::reset_serial(
        models[i_md].atoms(interleaved_conf).const_ref(),
        first_value);
    }
  }

namespace atoms {

  void
  reset_i_seq(
    af::const_ref<atom> const& atoms)
  {
    unsigned value = 0;
    for(const atom* a=atoms.begin();a!=atoms.end();a++) {
      a->data->i_seq = value++;
    }
  }

  std::size_t
  set_chemical_element_simple_if_necessary(
    af::ref<atom> const& atoms,
    bool tidy_existing)
  {
    std::size_t result = 0;
    for(atom* a=atoms.begin();a!=atoms.end();a++) {
      if (a->set_chemical_element_simple_if_necessary(tidy_existing)) {
        result++;
      }
    }
    return result;
  }

  atom_tmp_sentinel::atom_tmp_sentinel(
    af::const_ref<atom> const& atoms)
  :
    atoms_(atoms.begin(), atoms.end()) // copy
  {
    typedef std::vector<atom>::iterator it;
    it i_end = atoms_.end();
    for(it i=atoms_.begin();i!=i_end;i++) {
      atom_data* d = i->data.get();
      if (d->have_sentinel) {
        throw std::runtime_error(
          "Another associated atom_tmp_sentinel instance still exists.");
      }
      d->have_sentinel = true;
    }
  }

  atom_tmp_sentinel::~atom_tmp_sentinel()
  {
    typedef std::vector<atom>::iterator it;
    it i_end = atoms_.end();
    for(it i=atoms_.begin();i!=i_end;i++) {
      atom_data* d = i->data.get();
      d->tmp = 0;
      d->have_sentinel = false;
    }
  }

  std::auto_ptr<atom_tmp_sentinel>
  reset_tmp(
    af::const_ref<atom> const& atoms,
    int first_value,
    int increment)
  {
    std::auto_ptr<atom_tmp_sentinel> result(new atom_tmp_sentinel(atoms));
    int value = first_value;
    for(const atom* a=atoms.begin();a!=atoms.end();a++) {
      a->data->tmp = value;
      value += increment;
    }
    return result;
  }

  std::auto_ptr<atom_tmp_sentinel>
  reset_tmp_for_occupancy_groups_simple(
    af::const_ref<atom> const& atoms)
  {
    std::auto_ptr<atom_tmp_sentinel> result(new atom_tmp_sentinel(atoms));
    int value = 0;
    for(const atom* a=atoms.begin();a!=atoms.end();a++,value++) {
      a->data->tmp = (a->element_is_hydrogen() ? -1 : value);
    }
    return result;
  }

}}}} // namespace iotbx::pdb::hierarchy::atoms
