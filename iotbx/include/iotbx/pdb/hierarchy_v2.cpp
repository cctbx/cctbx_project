#include <iotbx/pdb/hierarchy_v2.h>
#include <iotbx/pdb/common_residue_names.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <cctbx/eltbx/chemical_elements.h>
#include <boost/scoped_array.hpp>
#include <ctype.h>

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

  root
  root::deep_copy() const
  {
    root result;
    detach_copy_children(result, result.data->models, data->models);
    return result;
  }

  void
  root::new_models(unsigned number_of_additional_models)
  {
    pre_allocate_models(number_of_additional_models);
    for(unsigned i=0;i<number_of_additional_models;i++) {
      data->models.push_back(model(*this));
    }
  }

  unsigned
  root::reset_atom_tmp(int new_value) const
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

  model::model(
    root const& parent,
    model const& other)
  :
    data(new model_data(parent.data, *other.data))
  {
    detach_copy_children(*this, data->chains, other.data->chains);
  }

  model
  model::detached_copy() const
  {
    model result;
    detach_copy_children(result, result.data->chains, data->chains);
    return result;
  }

  model&
  model::set_parent(root const& parent)
  {
    bool root_has_parent = (data->parent.lock().get() != 0);
    SCITBX_ASSERT(!root_has_parent);
    data->parent = parent.data;
    return *this;
  }

  boost::optional<root>
  model::parent() const
  {
    shared_ptr<root_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<root>();
    return boost::optional<root>(root(parent, true));
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

  chain::chain(
    model const& parent,
    chain const& other)
  :
    data(new chain_data(parent.data, *other.data))
  {
    detach_copy_children(
      *this, data->residue_groups, other.data->residue_groups);
  }

  chain
  chain::detached_copy() const
  {
    chain result;
    detach_copy_children(
      result, result.data->residue_groups, data->residue_groups);
    return result;
  }

  chain&
  chain::set_parent(model const& parent)
  {
    bool chain_has_parent = (data->parent.lock().get() != 0);
    SCITBX_ASSERT(!chain_has_parent);
    data->parent = parent.data;
    return *this;
  }

  boost::optional<model>
  chain::parent() const
  {
    shared_ptr<model_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<model>();
    return boost::optional<model>(model(parent, true));
  }

  void
  chain::new_residue_groups(unsigned number_of_additional_residue_groups)
  {
    pre_allocate_residue_groups(number_of_additional_residue_groups);
    for(unsigned i=0;i<number_of_additional_residue_groups;i++) {
      data->residue_groups.push_back(residue_group(*this));
    }
  }

  unsigned
  chain::reset_atom_tmp(int new_value) const
  {
    unsigned n = residue_groups_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const residue_group* rg = &*data->residue_groups.begin();
    for(unsigned i=0;i<n;i++) {
      result += (rg++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  residue_group::residue_group(
    chain const& parent,
    residue_group const& other)
  :
    data(new residue_group_data(parent.data, *other.data))
  {
    detach_copy_children(*this, data->atom_groups, other.data->atom_groups);
  }

  residue_group
  residue_group::detached_copy() const
  {
    residue_group result;
    detach_copy_children(result, result.data->atom_groups, data->atom_groups);
    return result;
  }

  residue_group&
  residue_group::set_parent(chain const& parent)
  {
    bool residue_group_has_parent = (data->parent.lock().get() != 0);
    SCITBX_ASSERT(!residue_group_has_parent);
    data->parent = parent.data;
    return *this;
  }

  boost::optional<chain>
  residue_group::parent() const
  {
    shared_ptr<chain_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<chain>();
    return boost::optional<chain>(chain(parent, true));
  }

  void
  residue_group::new_atom_groups(unsigned number_of_additional_atom_groups)
  {
    pre_allocate_atom_groups(number_of_additional_atom_groups);
    for(unsigned i=0;i<number_of_additional_atom_groups;i++) {
      data->atom_groups.push_back(atom_group(*this));
    }
  }

  unsigned
  residue_group::reset_atom_tmp(int new_value) const
  {
    unsigned n = atom_groups_size();
    if (n == 0) return 0;
    unsigned result = 0;
    const atom_group* ag = &*data->atom_groups.begin();
    for(unsigned i=0;i<n;i++) {
      result += (ag++)->reset_atom_tmp(new_value);
    }
    return result;
  }

  atom_group::atom_group(
    residue_group const& parent,
    atom_group const& other)
  :
    data(new atom_group_data(parent.data, *other.data))
  {
    detach_copy_children(*this, data->atoms, other.data->atoms);
  }

  atom_group
  atom_group::detached_copy() const
  {
    atom_group result;
    detach_copy_children(result, result.data->atoms, data->atoms);
    return result;
  }

  atom_group&
  atom_group::set_parent(residue_group const& parent)
  {
    bool atom_group_has_parent = (data->parent.lock().get() != 0);
    SCITBX_ASSERT(!atom_group_has_parent);
    data->parent = parent.data;
    return *this;
  }

  boost::optional<residue_group>
  atom_group::parent() const
  {
    shared_ptr<residue_group_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<residue_group>();
    return boost::optional<residue_group>(residue_group(parent, true));
  }

  void
  atom_group::new_atoms(unsigned number_of_additional_atoms)
  {
    pre_allocate_atoms(number_of_additional_atoms);
    for(unsigned i=0;i<number_of_additional_atoms;i++) {
      data->atoms.push_back(atom(*this));
    }
  }

  long
  atom_group::find_atom_index(
    hierarchy_v2::atom const& atom,
    bool must_be_present) const
  {
    long n = static_cast<long>(data->atoms.size());
    for(long i=0;i<n;i++) {
      if (data->atoms[i].data.get() == atom.data.get()) return i;
    }
    if (must_be_present) {
      throw std::runtime_error("atom not in atom_group.");
    }
    return -1;
  }

  unsigned
  atom_group::reset_atom_tmp(int new_value) const
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

  atom&
  atom::set_parent(atom_group const& parent)
  {
    bool atom_has_parent = (data->parent.lock().get() != 0);
    SCITBX_ASSERT(!atom_has_parent);
    data->parent = parent.data;
    return *this;
  }

  boost::optional<atom_group>
  atom::parent() const
  {
    shared_ptr<atom_group_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<atom_group>();
    return boost::optional<atom_group>(atom_group(parent, true));
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
