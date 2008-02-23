#include <iotbx/pdb/hierarchy.h>
#include <iotbx/pdb/common_residue_names.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <cctbx/eltbx/chemical_elements.h>
#include <boost/scoped_array.hpp>
#include <ctype.h>

namespace iotbx { namespace pdb {

  void
  hierarchy::new_models(unsigned number_of_additional_models)
  {
    pre_allocate_models(number_of_additional_models);
    for(unsigned i=0;i<number_of_additional_models;i++) {
      data->models.push_back(model(*this));
    }
  }

  unsigned
  hierarchy::reset_atom_tmp(int new_value) const
  {
    unsigned n = models_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const model* m = &*data->models.begin();
    for(unsigned i=0;i<n;i++) {
      result += (m++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  boost::shared_ptr<hierarchy::overall_counts_holder>
  hierarchy::overall_counts() const
  {
    reset_atom_tmp(0);
    boost::shared_ptr<overall_counts_holder>
      result(new overall_counts_holder);
    std::map<str3, unsigned> residue_name_counts_str3;
    unsigned ms_sz = models_size();
    result->n_models = ms_sz;
    for(unsigned jm=0;jm<ms_sz;jm++) {
      model const& m = models()[jm];
      unsigned cs_sz = m.chains_size();
      result->n_chains += cs_sz;
      for(unsigned jc=0;jc<cs_sz;jc++) {
        chain const& c = m.chains()[jc];
        result->chain_ids[c.data->id]++;
        unsigned fs_sz = c.conformers_size();
        result->n_conformers += fs_sz;
        for(unsigned jf=0;jf<fs_sz;jf++) {
          conformer const& f = c.conformers()[jf];
          result->conformer_ids[f.data->id]++;
          for(unsigned jr=0;jr<f.residues().size();jr++) {
            residue const& r = f.residues()[jr];
            unsigned n_atoms = 0;
            unsigned as_sz = r.atoms_size();
            const atom* a = (as_sz == 0 ? 0 : &*r.atoms().begin());
            for(unsigned ja=0;ja<as_sz;ja++,a++) {
              if (a->data->tmp == 0) {
                a->data->tmp = 1;
                n_atoms++;
              }
            }
            if (n_atoms != 0) {
              result->n_atoms += n_atoms;
              residue_name_counts_str3[r.data->name]++;
            }
          }
        }
      }
    }
    for(std::map<str3, unsigned>::const_iterator
          i=residue_name_counts_str3.begin();
          i!=residue_name_counts_str3.end();i++) {
      result->n_residues += i->second;
      result->residue_names[i->first.elems] = i->second;
      result->residue_name_classes[
        common_residue_names::get_class(i->first)] += i->second;
    }
    return result;
  }

  boost::optional<hierarchy>
  model::parent() const
  {
    shared_ptr<hierarchy_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<hierarchy>();
    return boost::optional<hierarchy>(hierarchy(parent, true));
  }

  void
  model::new_chains(unsigned number_of_additional_chains)
  {
    pre_allocate_chains(number_of_additional_chains);
    for(unsigned i=0;i<number_of_additional_chains;i++) {
      data->chains.push_back(chain(*this));
    }
  }

  unsigned
  model::reset_atom_tmp(int new_value) const
  {
    unsigned n = chains_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const chain* c = &*data->chains.begin();
    for(unsigned i=0;i<n;i++) {
      result += (c++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  boost::optional<model>
  chain::parent() const
  {
    shared_ptr<model_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<model>();
    return boost::optional<model>(model(parent, true));
  }

  void
  chain::new_conformers(unsigned number_of_additional_conformers)
  {
    pre_allocate_conformers(number_of_additional_conformers);
    for(unsigned i=0;i<number_of_additional_conformers;i++) {
      data->conformers.push_back(conformer(*this));
    }
  }

  unsigned
  chain::number_of_atoms() const
  {
    unsigned n = conformers_size();
    if (n == 0) return 0;
    const conformer* f = &*data->conformers.begin();
    unsigned result = (f++)->number_of_atoms();
    for(unsigned i=1;i<n;i++) {
      result += (f++)->number_of_alternative_atoms();
    }
    return result;
  }

  unsigned
  chain::reset_atom_tmp(int new_value) const
  {
    unsigned n = conformers_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const conformer* f = &*data->conformers.begin();
    for(unsigned i=0;i<n;i++) {
      result += (f++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  af::shared<vec3>
  chain::extract_sites_cart() const
  {
    reset_atom_tmp(1);
    af::shared<vec3> result((af::reserve(reset_atom_tmp(0))));
    unsigned fs_sz = conformers_size();
    for(unsigned jf=0;jf<fs_sz;jf++) {
      conformer const& f = data->conformers[jf];
      unsigned rs_sz = f.residues_size();
      for(unsigned jr=0;jr<rs_sz;jr++) {
        residue const& r = f.residues()[jr];
        unsigned as_sz = r.atoms_size();
        const atom* a = (as_sz == 0 ? 0 : &*r.atoms().begin());
        for(unsigned ja=0;ja<as_sz;ja++,a++) {
          if (a->data->tmp == 0) {
            a->data->tmp = 1;
            result.push_back(a->data->xyz);
          }
        }
      }
    }
    return result;
  }

  void
  combine_conformers::process_single_conformer(
    pdb::conformer const& conformer)
  {
    std::string const& conformer_id = conformer.data->id;
    std::vector<residue> const& residues = conformer.residues();
    unsigned n = static_cast<unsigned>(residues.size());
    const residue* r = (n == 0 ? 0 : &*residues.begin());
    for(unsigned i=0;i<n;i++) {
      process_next_residue(conformer_id, true, *r++, (i==0));
    }
  }

  void
  combine_conformers::process_conformers(
    std::vector<conformer> const& conformers)
  {
    unsigned n_conf = static_cast<unsigned>(conformers.size());
    if (n_conf == 1) { // shortcut to increase efficiency
      process_single_conformer(conformers[0]);
      return;
    }
    boost::scoped_array<unsigned>
      i_residue_to_do_owner(new unsigned[n_conf]);
    unsigned* i_residue_to_do = i_residue_to_do_owner.get();
    std::fill_n(i_residue_to_do, n_conf, 0U);
    std::map<const residue_data*, unsigned> residue_todo_data_ptrs;
    for(unsigned i_conf=0;i_conf<n_conf;i_conf++) {
      std::vector<residue> const& residues = conformers[i_conf].residues();
      if (residues.size() == 0) continue;
      residue_todo_data_ptrs[residues[0].data.get()] = i_conf;
    }
    while (residue_todo_data_ptrs.size() != 0) {
      std::map<const residue_data*, unsigned>::const_iterator
        residue_todo_data_ptrs_end = residue_todo_data_ptrs.end();
      // find a conformer that shares atoms only with pivot residues
      typedef std::set<unsigned> i_confs_t;
      i_confs_t i_confs;
      for(unsigned i_conf=0;i_conf<n_conf;i_conf++) {
        std::vector<residue> const&
          residues = conformers[i_conf].residues();
        unsigned i_res = i_residue_to_do[i_conf];
        if (i_res == residues.size()) continue;
        // test if this conformer shares atoms only with pivot residues
        i_confs.clear();
        i_confs.insert(i_conf);
        typedef std::vector<atom> atoms_t;
        atoms_t const& atoms = residues[i_res].atoms();
        atoms_t::const_iterator atoms_end = atoms.end();
        for(atoms_t::const_iterator atom = atoms.begin();
            atom != atoms_end;
            atom++) {
          if (atom->is_alternative()) continue;
          typedef std::vector<weak_ptr<residue_data> > parents_t;
          parents_t const& parents = atom->data->parents;
          parents_t::const_iterator parents_end = parents.end();
          for(parents_t::const_iterator parent = parents.begin();
              parent != parents_end;
              parent++) {
            const residue_data* p = parent->lock().get();
            std::map<const residue_data*, unsigned>::const_iterator
              m = residue_todo_data_ptrs.find(p);
            if (m == residue_todo_data_ptrs_end) {
              goto try_next_i_conf;
            }
            i_confs.insert(m->second);
          }
        }
        goto process_i_confs;
        try_next_i_conf:;
      }
      throw std::runtime_error("conformers form a knot.");
      process_i_confs:
      bool is_first_of_group = true;
      i_confs_t::const_iterator i_confs_end = i_confs.end();
      for(i_confs_t::const_iterator i_conf_iter=i_confs.begin();
          i_conf_iter != i_confs_end;
          i_conf_iter++) {
        unsigned i_conf = *i_conf_iter;
        conformer const& c = conformers[i_conf];
        unsigned i_res = i_residue_to_do[i_conf]++;
        std::vector<residue>::const_iterator
          r = c.residues().begin() + i_res;
        process_next_residue(
          c.data->id, is_first_of_group, *r, (i_res == 0));
        is_first_of_group = false;
        residue_todo_data_ptrs.erase(r->data.get());
        if (++i_res != c.residues_size()) {
          residue_todo_data_ptrs[(++r)->data.get()] = i_conf;
        }
      }
    }
  }

  bool
  conformer_data::parent_chain_has_multiple_conformers() const
  {
    shared_ptr<chain_data> p = parent.lock();
    if (p.get() == 0) return false;
    return p->has_multiple_conformers();
  }

  boost::optional<chain>
  conformer::parent() const
  {
    shared_ptr<chain_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<chain>();
    return boost::optional<chain>(chain(parent, true));
  }

  void
  conformer::new_residues(unsigned number_of_additional_residues)
  {
    pre_allocate_residues(number_of_additional_residues);
    for(unsigned i=0;i<number_of_additional_residues;i++) {
      data->residues.push_back(residue(*this));
    }
  }

  unsigned
  conformer::number_of_atoms() const
  {
    unsigned n = residues_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const residue* r = &*data->residues.begin();
    for(unsigned i=0;i<n;i++) {
      result += (r++)->atoms().size();
    }
    return result;
  }

  unsigned
  conformer::number_of_alternative_atoms() const
  {
    unsigned n = residues_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const residue* r = &*data->residues.begin();
    for(unsigned i=0;i<n;i++) {
      result += (r++)->number_of_alternative_atoms();
    }
    return result;
  }

  unsigned
  conformer::reset_atom_tmp(int new_value) const
  {
    unsigned n = residues_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const residue* r = &*data->residues.begin();
    for(unsigned i=0;i<n;i++) {
      result += (r++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  af::shared<vec3>
  conformer::residue_centers_of_geometry() const
  {
    unsigned n = residues_size();
    af::shared<vec3> result((af::reserve(n)));
    const residue* r = (n == 0 ? 0 : &*data->residues.begin());
    for(unsigned i=0;i<n;i++) {
      result.push_back((r++)->center_of_geometry());
    }
    return result;
  }

  void
  conformer::select_residues_in_place(
    af::const_ref<std::size_t> const& selection)
  {
    unsigned ns = static_cast<unsigned>(selection.size());
    unsigned nr = residues_size();
    std::vector<residue> remaining_residues;
    remaining_residues.reserve(ns);
    const residue* r = (nr == 0 ? 0 : &*data->residues.begin());
    for(unsigned i=0;i<ns;i++) {
      std::size_t j = selection[i];
      if (j >= nr) {
        throw std::runtime_error("selection index out of range.");
      }
      remaining_residues.push_back(r[j]);
    }
    data->residues.swap(remaining_residues);
  }

  conformer
  conformer::select_residues(
    af::const_ref<std::size_t> const& selection)
  {
    unsigned ns = static_cast<unsigned>(selection.size());
    unsigned nr = residues_size();
    conformer result(data->id);
    result.data->residues.reserve(ns);
    const residue* r = (nr == 0 ? 0 : &*data->residues.begin());
    for(unsigned i=0;i<ns;i++) {
      std::size_t j = selection[i];
      if (j >= nr) {
        throw std::runtime_error("selection index out of range.");
      }
      result.data->residues.push_back(residue(result, r[j]));
    }
    return result;
  }

  af::shared<std::size_t>
  conformer::residue_class_selection(
    std::set<str3> const& class_set,
    bool negate) const
  {
    unsigned n = residues_size();
    af::shared<std::size_t> result((af::reserve(n)));
    std::set<str3>::const_iterator class_set_end = class_set.end();
    const residue* r = (n == 0 ? 0 : &*data->residues.begin());
    if (!negate) {
      for(unsigned i=0;i<n;i++) {
        if (class_set.find((r++)->data->name) != class_set_end) {
          result.push_back(i);
        }
      }
    }
    else {
      for(unsigned i=0;i<n;i++) {
        if (class_set.find((r++)->data->name) == class_set_end) {
          result.push_back(i);
        }
      }
    }
    return result;
  }

  af::shared<std::size_t>
  conformer::residue_class_selection(
    std::string const& class_name,
    bool negate) const
  {
    if (class_name == "common_amino_acid")
      return residue_class_selection(
        common_residue_names::amino_acid_set(), negate);
    if (class_name == "common_rna_dna")
      return residue_class_selection(
        common_residue_names::rna_dna_set(), negate);
    if (class_name == "ccp4_mon_lib_rna_dna")
      return residue_class_selection(
        common_residue_names::ccp4_mon_lib_rna_dna_set(), negate);
    if (class_name == "common_water")
      return residue_class_selection(
        common_residue_names::water_set(), negate);
    if (class_name == "common_small_molecule")
      return residue_class_selection(
        common_residue_names::small_molecule_set(), negate);
    if (class_name == "common_element")
      return residue_class_selection(
        common_residue_names::element_set(), negate);
    throw std::runtime_error("unknown class_name=\""+class_name+"\"");
  }

  af::shared<std::size_t>
  conformer::residue_class_selection(
    af::const_ref<std::string> const& residue_names,
    bool negate) const
  {
    std::set<str3> class_set;
    common_residue_names::initialize_set(class_set, residue_names);
    return residue_class_selection(class_set, negate);
  }

  void
  conformer::select_residue_class_in_place(
    std::string const& class_name,
    bool negate)
  {
    af::shared<std::size_t>
      sel = residue_class_selection(class_name, negate);
    if (sel.size() != residues().size()) {
      select_residues_in_place(sel.const_ref());
    }
  }

  void
  conformer::select_residue_class_in_place(
    af::const_ref<std::string> const& residue_names,
    bool negate)
  {
    af::shared<std::size_t>
      sel = residue_class_selection(residue_names, negate);
    if (sel.size() != residues().size()) {
      select_residues_in_place(sel.const_ref());
    }
  }

  std::string
  residue_data::id() const
  {
    char buf[32];
    std::sprintf(buf, "%-3s%4s%s", name.elems, seq.elems, icode.elems);
    return std::string(buf);
  }

  bool
  residue_data::parent_chain_has_multiple_conformers() const
  {
    shared_ptr<conformer_data> p = parent.lock();
    if (p.get() == 0) return false;
    return p->parent_chain_has_multiple_conformers();
  }

  residue::residue(
    conformer const& parent,
    residue const& other)
  :
    data(new residue_data(parent.data, *other.data))
  {
    std::vector<atom> const& other_atoms = other.data->atoms;
    unsigned n = other.atoms_size();
    data->atoms.reserve(n);
    const atom* oa = (n == 0 ? 0 : &*other.atoms().begin());
    for(unsigned i=0;i<n;i++) {
      data->atoms.push_back(atom(*this, *oa++));
    }
  }

  boost::optional<conformer>
  residue::parent() const
  {
    shared_ptr<conformer_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<conformer>();
    return boost::optional<conformer>(conformer(parent, true));
  }

  void
  residue::new_atoms(unsigned number_of_additional_atoms)
  {
    pre_allocate_atoms(number_of_additional_atoms);
    for(unsigned i=0;i<number_of_additional_atoms;i++) {
      atom new_atom;
      new_atom.pre_allocate_parents(1U);
      new_atom.add_parent(*this);
      data->atoms.push_back(new_atom);
    }
  }

  af::shared<std::string>
  residue::atom_names() const
  {
    af::shared<std::string> result;
    unsigned n = atoms_size();
    if (n == 0) return result;
    const atom* a = &*data->atoms.begin();
    for(unsigned i=0;i<n;i++) {
      result.push_back(std::string((a++)->data->name.elems));
    }
    return result;
  }

  unsigned
  residue::number_of_alternative_atoms() const
  {
    unsigned n = atoms_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const atom* a = &*data->atoms.begin();
    for(unsigned i=0;i<n;i++) {
      if ((a++)->is_alternative()) result++;
    }
    return result;
  }

  unsigned
  residue::reset_atom_tmp(int new_value) const
  {
    unsigned n = atoms_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const atom* a = &*data->atoms.begin();
    for(unsigned i=0;i<n;i++,a++) {
      if (a->data->tmp != new_value) {
        a->data->tmp = new_value;
        result++;
      }
    }
    return result;
  }

  vec3
  residue::center_of_geometry() const
  {
    vec3 result(0,0,0);
    unsigned n = atoms_size();
    const atom* a = (n == 0 ? 0 : &*data->atoms.begin());
    for(unsigned i=0;i<n;i++) {
      result += (a++)->data->xyz;
    }
    if (n != 0) result /= static_cast<double>(n);
    return result;
  }

  atom_data::atom_data(
    weak_ptr<residue_data> const& parent_,
    atom_data const& other)
  :
    parents(1, parent_),
    name(other.name),
    segid(other.segid),
    element(other.element),
    charge(other.charge),
    xyz(other.xyz), sigxyz(other.sigxyz),
    occ(other.occ), sigocc(other.sigocc),
    b(other.b), sigb(other.sigb),
    uij(other.uij), siguij(other.siguij),
    hetero(other.hetero),
    flag_altloc(other.flag_altloc),
    tmp(0)
  {}

  bool
  atom_data::is_alternative() const
  {
    if (parents.size() != 1) return false;
    shared_ptr<residue_data> p0 = parents[0].lock();
    if (p0.get() == 0) return false;
    return p0->parent_chain_has_multiple_conformers();
  }

  unsigned
  atom::reset_parents()
  {
    if (data->parents.capacity() == 0) return 0U;
    unsigned result = parents_size();
    std::vector<weak_ptr<residue_data> > buffer;
    data->parents.swap(buffer);
    return result;
  }

  unsigned
  atom::parents_size() const
  {
    unsigned result = 0;
    unsigned n = static_cast<unsigned>(data->parents.size());
    for(unsigned i=0;i<n;i++) {
      shared_ptr<residue_data> parent = data->parents[i].lock();
      if (parent.get() != 0) result++;
    }
    return result;
  }

  af::shared<residue>
  atom::parents() const
  {
    af::shared<residue> result;
    unsigned n = static_cast<unsigned>(data->parents.size());
    result.reserve(n);
    for(unsigned i=0;i<n;i++) {
      shared_ptr<residue_data> parent = data->parents[i].lock();
      if (parent.get() != 0) result.push_back(residue(parent, true));
    }
    return result;
  }

  void
  atom::add_parent(residue const& new_parent)
  {
    unsigned n = static_cast<unsigned>(data->parents.size());
    for(unsigned i=0;i<n;i++) {
      shared_ptr<residue_data> parent = data->parents[i].lock();
      if (parent.get() == 0) {
        data->parents[i] = new_parent.data;
        return;
      }
    }
    data->parents.push_back(new_parent.data);
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

  unsigned
  atom::format_atom_record(
    char* result,
    int serial,
    const char* resname,
    const char* resseq,
    const char* icode,
    const char* altloc,
    const char* chain_id) const
  {
    char blank = ' ';
    atom_data const& d = *data;
    std::memcpy(result, (d.hetero ? "HETATM" : "ATOM  "), 6U);
    const char* errmsg = hy36encode(5U, serial, result+6);
    if (errmsg != 0) throw std::runtime_error(
      std::string("hy36encode: ") + errmsg);
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
    char* result,
    int serial) const
  {
    atom_data const& d = *data;
    boost::shared_ptr<conformer_data> f_lock;
    boost::shared_ptr<chain_data> c_lock;
    unsigned n = static_cast<unsigned>(d.parents.size());
    for(unsigned i=0;i<n;i++) {
      shared_ptr<residue_data> res1 = data->parents[i].lock();
      const residue_data* r1 = res1.get();
      if (r1 != 0) {
        f_lock = r1->parent.lock();
        const conformer_data* f = f_lock.get();
        const chain_data* c = 0;
        if (f != 0) {
          c_lock = f->parent.lock();
          c = c_lock.get();
        }
        for(unsigned j=i+1;j<n;j++) {
          shared_ptr<residue_data> res2 = data->parents[j].lock();
          if (res2.get() != 0) {
            return format_atom_record(
              result,
              serial,
              r1->name.elems, r1->seq.elems, r1->icode.elems,
              0,
              (c != 0 ? c->id.c_str() : 0));
          }
        }
        return format_atom_record(
          result,
          serial,
          r1->name.elems, r1->seq.elems, r1->icode.elems,
          (f != 0 ? f->id.c_str() : 0),
          (c != 0 ? c->id.c_str() : 0));
      }
    }
    return format_atom_record(result, serial, 0, 0, 0, 0, 0);
  }

}} // namespace iotbx::pdb
