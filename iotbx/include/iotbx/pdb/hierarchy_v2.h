#ifndef IOTBX_PDB_HIERARCHY_V2_H
#define IOTBX_PDB_HIERARCHY_V2_H

#include <iotbx/pdb/namespace.h>
#include <iotbx/pdb/small_str.h>
#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <map>
#include <set>
#include <vector>
#include <string>

namespace iotbx { namespace pdb {

//! Hierarchy organizing PDB atom objects.
/*! Data hierarchy:
<pre>
  model
    parent
    id
    chain(s)
      parent
      id
      residue_group(s)
        parent
        resseq
        icode
        link_to_previous
        atom_group(s)
          parent
          altloc
          resname
          atom(s)
            parent
            name
            segid
            element
            charge
            serial
            xyz sigxyz
            occ sigocc
            b sigb
            uij siguij
            hetero
</pre>
Convenience objects:
<pre>
  residue
    altloc
    name
    resseq
    icode
    atoms
</pre>
A residue object is NOT a parent of the atoms.
<p>
<pre>
  conformer
    altloc
    residues
</pre>
 */

namespace hierarchy_v2 {

  class root_data;
  class root;
  class model_data;
  class model;
  class chain_data;
  class chain;
  class residue_group_data;
  class residue_group;
  class atom_group_data;
  class atom_group;
  class atom_data;
  class atom;

  class conformer;
  class residue;

  //! Holder for root attributes (to be held by a shared_ptr).
  class root_data : boost::noncopyable
  {
    protected:
      friend class root;
    public:
      af::shared<std::string> info;
    protected:
      std::vector<model> models;
  };

  //! Holder for model attributes (to be held by a shared_ptr).
  class model_data : boost::noncopyable
  {
    protected:
      friend class model;
      weak_ptr<root_data> parent;
    public:
      std::string id;
    protected:
      std::vector<chain> chains;

      inline
      model_data(
        weak_ptr<root_data> const& parent,
        std::string const& id);

      inline
      model_data(
        std::string const& id);

      inline
      model_data(
        weak_ptr<root_data> const& parent_,
        model_data const& other);
  };

  //! Holder for chain attributes (to be held by a shared_ptr).
  class chain_data : boost::noncopyable
  {
    protected:
      friend class chain;
      weak_ptr<model_data> parent;
    public:
      std::string id;
    protected:
      std::vector<residue_group> residue_groups;

      inline
      chain_data(
        weak_ptr<model_data> const& parent,
        std::string const& id);

      inline
      chain_data(
        std::string const& id);

      inline
      chain_data(
        weak_ptr<model_data> const& parent_,
        chain_data const& other);
  };

  //! Holder for residue_group attributes (to be held by a shared_ptr).
  class residue_group_data : boost::noncopyable
  {
    protected:
      friend class residue_group;
      weak_ptr<chain_data> parent;
    public:
      str4 resseq;
      str1 icode;
      bool link_to_previous;
    protected:
      std::vector<atom_group> atom_groups;

      inline
      residue_group_data(
        weak_ptr<chain_data> const& parent,
        const char* resseq_,
        const char* icode_,
        bool link_to_previous_);

      inline
      residue_group_data(
        const char* resseq_,
        const char* icode_,
        bool link_to_previous_);

      inline
      residue_group_data(
        weak_ptr<chain_data> const& parent_,
        residue_group_data const& other);
  };

  //! Holder for atom_group attributes (to be held by a shared_ptr).
  class atom_group_data : boost::noncopyable
  {
    protected:
      friend class atom_group;
      weak_ptr<residue_group_data> parent;
    public:
      str1 altloc;
      str3 resname;
    protected:
      std::vector<atom> atoms;

      inline
      atom_group_data(
        weak_ptr<residue_group_data> const& parent_,
        const char* altloc_,
        const char* resname_);

      inline
      atom_group_data(
        const char* altloc_,
        const char* resname_);

      inline
      atom_group_data(
        weak_ptr<residue_group_data> const& parent_,
        atom_group_data const& other);
  };

  //! Holder for atom attributes (to be held by a shared_ptr).
  class atom_data : boost::noncopyable
  {
    protected:
      friend class atom;
      friend class combine_conformers;
      weak_ptr<atom_group_data> parent;
    public:
      str4 name;
      str4 segid;
      str2 element;
      str2 charge;
      str5 serial;
      vec3 xyz;
      vec3 sigxyz;
      double occ;
      double sigocc;
      double b;
      double sigb;
      sym_mat3 uij;
      sym_mat3 siguij;
      bool hetero;
      mutable int tmp;

    protected:
      atom_data(
        weak_ptr<atom_group_data> const& parent_,
        str4 name_, str4 segid_,
        str2 element_, str2 charge_, str5 serial_,
        vec3 const& xyz_, vec3 const& sigxyz_,
        double occ_, double sigocc_,
        double b_, double sigb_,
        sym_mat3 const& uij_, sym_mat3 const& siguij_,
        bool hetero_)
      :
        parent(parent_),
        name(name_), segid(segid_),
        element(element_), charge(charge_), serial(serial_),
        xyz(xyz_), sigxyz(sigxyz_),
        occ(occ_), sigocc(sigocc_),
        b(b_), sigb(sigb_),
        uij(uij_), siguij(siguij_),
        hetero(hetero_),
        tmp(0)
      {}

      atom_data(
        str4 name_, str4 segid_,
        str2 element_, str2 charge_, str5 serial_,
        vec3 const& xyz_, vec3 const& sigxyz_,
        double occ_, double sigocc_,
        double b_, double sigb_,
        sym_mat3 const& uij_, sym_mat3 const& siguij_,
        bool hetero_)
      :
        name(name_), segid(segid_),
        element(element_), charge(charge_), serial(serial_),
        xyz(xyz_), sigxyz(sigxyz_),
        occ(occ_), sigocc(sigocc_),
        b(b_), sigb(sigb_),
        uij(uij_), siguij(siguij_),
        hetero(hetero_),
        tmp(0)
      {}

      atom_data(
        weak_ptr<atom_group_data> const& parent_,
        atom_data const& other);
  };

  //! Atom attributes.
  class atom
  {
    public:
      shared_ptr<atom_data> data;

    protected:
      friend class atom_group; // to support atom_group::append_atom()

      atom&
      set_parent(atom_group const& parent);

      void
      clear_parent() { data->parent.reset(); }

      unsigned
      reset_atom_tmp(int new_value) const
      {
        data->tmp = new_value;
        return 1U;
      }

    public:
      atom(shared_ptr<atom_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      atom(
        atom_group const& parent,
        const char* name="", const char* segid="",
        const char* element="", const char* charge="", const char* serial="",
        vec3 const& xyz=vec3(0,0,0), vec3 const& sigxyz=vec3(0,0,0),
        double occ=0, double sigocc=0,
        double b=0, double sigb=0,
        sym_mat3 const& uij=sym_mat3(-1,-1,-1,-1,-1,-1),
        sym_mat3 const& siguij=sym_mat3(-1,-1,-1,-1,-1,-1),
        bool hetero=false);

      atom(
        const char* name="", const char* segid="",
        const char* element="", const char* charge="", const char* serial="",
        vec3 const& xyz=vec3(0,0,0), vec3 const& sigxyz=vec3(0,0,0),
        double occ=0, double sigocc=0,
        double b=0, double sigb=0,
        sym_mat3 const& uij=sym_mat3(-1,-1,-1,-1,-1,-1),
        sym_mat3 const& siguij=sym_mat3(-1,-1,-1,-1,-1,-1),
        bool hetero=false)
      :
        data(new atom_data(
          name, segid,
          element, charge, serial,
          xyz, sigxyz,
          occ, sigocc,
          b, sigb,
          uij, siguij,
          hetero))
      {}

      inline
      atom(
        atom_group const& parent,
        atom const& other);

      atom
      detached_copy() const;

      atom&
      set_name(const char* new_name)
      {
        data->name = new_name;
        return *this;
      }

      atom&
      set_segid(const char* new_segid)
      {
        data->segid = new_segid;
        return *this;
      }

      atom&
      set_element(const char* new_element)
      {
        data->element = new_element;
        return *this;
      }

      atom&
      set_charge(const char* new_charge)
      {
        data->charge = new_charge;
        return *this;
      }

      atom&
      set_serial(const char* new_serial)
      {
        data->serial = new_serial;
        return *this;
      }

      atom&
      set_xyz(vec3 const& new_xyz)
      {
        data->xyz = new_xyz;
        return *this;
      }

      atom&
      set_sigxyz(vec3 const& new_sigxyz)
      {
        data->sigxyz = new_sigxyz;
        return *this;
      }

      atom&
      set_occ(double new_occ)
      {
        data->occ = new_occ;
        return *this;
      }

      atom&
      set_sigocc(double new_sigocc)
      {
        data->sigocc = new_sigocc;
        return *this;
      }

      atom&
      set_b(double new_b)
      {
        data->b = new_b;
        return *this;
      }

      atom&
      set_sigb(double new_sigb)
      {
        data->sigb = new_sigb;
        return *this;
      }

      atom&
      set_uij(sym_mat3 const& new_uij)
      {
        data->uij = new_uij;
        return *this;
      }

      atom&
      set_siguij(sym_mat3 const& new_siguij)
      {
        data->siguij = new_siguij;
        return *this;
      }

      atom&
      set_hetero(double new_hetero)
      {
        data->hetero = new_hetero;
        return *this;
      }

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      static
      std::size_t
      sizeof_data() { return sizeof(atom_data); }

      boost::optional<atom_group>
      parent() const;

      bool
      uij_is_defined() const
      {
        return !data->uij.const_ref().all_eq(-1);
      }

      bool
      siguij_is_defined() const
      {
        return !data->siguij.const_ref().all_eq(-1);
      }

      boost::optional<std::string>
      determine_chemical_element_simple() const;
  };

  //! atom_group attributes.
  class atom_group
  {
    public:
      shared_ptr<atom_group_data> data;

    protected:
      friend class atom; // to support atom::parent()

      atom_group(
        shared_ptr<atom_group_data> const& data_, bool) : data(data_) {}

      friend class residue_group; // to support atom_group::append_atom_group()

      atom_group&
      set_parent(residue_group const& parent);

      void
      clear_parent() { data->parent.reset(); }

    public:
      atom_group(shared_ptr<atom_group_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      atom_group(
        residue_group const& parent,
        const char* altloc="",
        const char* resname="");

      atom_group(
        const char* altloc="",
        const char* resname="")
      :
        data(new atom_group_data(altloc, resname))
      {}

      atom_group(
        residue_group const& parent,
        atom_group const& other);

      atom_group
      detached_copy() const;

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<residue_group>
      parent() const;

      unsigned
      atoms_size() const;

      std::vector<atom> const&
      atoms() const;

      std::vector<atom>::iterator
      atoms_begin() const { return data->atoms.begin(); }

      long
      find_atom_index(
        hierarchy_v2::atom const& atom,
        bool must_be_present=false) const;

      void
      pre_allocate_atoms(unsigned number_of_additional_atoms);

      void
      new_atoms(unsigned number_of_additional_atoms);

      void
      insert_atom(long i, hierarchy_v2::atom& atom);

      void
      append_atom(hierarchy_v2::atom& atom);

      void
      remove_atom(long i);

      void
      remove_atom(hierarchy_v2::atom& atom);

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  //! residue_group attributes.
  class residue_group
  {
    public:
      shared_ptr<residue_group_data> data;

    protected:
      friend class atom_group; // to support atom_group::parent()

      residue_group(
        shared_ptr<residue_group_data> const& data_, bool) : data(data_) {}

      friend class chain; // to support chain::append_residue_group()

      residue_group&
      set_parent(chain const& parent);

      void
      clear_parent() { data->parent.reset(); }

    public:
      residue_group(shared_ptr<residue_group_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      residue_group(
        chain const& parent,
        const char* resseq="",
        const char* icode="",
        bool link_to_previous=true);

      inline
      residue_group(
        const char* resseq="",
        const char* icode="",
        bool link_to_previous=true)
      :
        data(new residue_group_data(resseq, icode, link_to_previous))
      {}

      residue_group(
        chain const& parent,
        residue_group const& other);

      residue_group
      detached_copy() const;

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<chain>
      parent() const;

      unsigned
      atom_groups_size() const;

      std::vector<atom_group> const&
      atom_groups() const;

      std::vector<atom_group>::iterator
      atom_groups_begin() const { return data->atom_groups.begin(); }

      long
      find_atom_group_index(
        hierarchy_v2::atom_group const& atom_group,
        bool must_be_present=false) const;

      void
      pre_allocate_atom_groups(unsigned number_of_additional_atom_groups);

      void
      new_atom_groups(unsigned number_of_additional_atom_groups);

      void
      insert_atom_group(long i, hierarchy_v2::atom_group& atom_group);

      void
      append_atom_group(hierarchy_v2::atom_group& atom_group);

      void
      remove_atom_group(long i);

      void
      remove_atom_group(hierarchy_v2::atom_group& atom_group);

      atom_group
      new_atom_group(
        const char* altloc="",
        const char* resname="")
      {
        data->atom_groups.push_back(
          atom_group(*this, altloc, resname));
        return data->atom_groups.back();
      }

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  //! Chain attributes.
  class chain
  {
    public:
      shared_ptr<chain_data> data;

    protected:
      friend class residue_group; // to support residue_group::parent()

      chain(shared_ptr<chain_data> const& data_, bool) : data(data_) {}

      friend class model; // to support model::append_chain()

      chain&
      set_parent(model const& parent);

      void
      clear_parent() { data->parent.reset(); }

    public:
      chain(shared_ptr<chain_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      chain(
        model const& parent,
        std::string const& id="");

      chain(
        std::string const& id="")
      :
        data(new chain_data(id))
      {}

      chain(
        model const& parent,
        chain const& other);

      chain
      detached_copy() const;

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<model>
      parent() const;

      unsigned
      residue_groups_size() const;

      std::vector<residue_group> const&
      residue_groups() const;

      std::vector<residue_group>::iterator
      residue_groups_begin() const { return data->residue_groups.begin(); }

      long
      find_residue_group_index(
        hierarchy_v2::residue_group const& residue_group,
        bool must_be_present=false) const;

      void
      pre_allocate_residue_groups(unsigned number_of_additional_residue_groups);

      void
      new_residue_groups(unsigned number_of_additional_residue_groups);

      void
      insert_residue_group(long i, hierarchy_v2::residue_group& residue_group);

      void
      append_residue_group(hierarchy_v2::residue_group& residue_group);

      void
      remove_residue_group(long i);

      void
      remove_residue_group(hierarchy_v2::residue_group& residue_group);

      residue_group
      new_residue_group(
        const char* resseq="",
        const char* icode="",
        bool link_to_previous=true)
      {
        data->residue_groups.push_back(
          residue_group(*this, resseq, icode, link_to_previous));
        return data->residue_groups.back();
      }

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  //! Model attributes.
  class model
  {
    public:
      shared_ptr<model_data> data;

    protected:
      friend class chain; // to support chain::parent()
      model(shared_ptr<model_data> const& data_, bool) : data(data_) {}

      friend class root; // to support root::append_model()

      model&
      set_parent(root const& parent);

      void
      clear_parent() { data->parent.reset(); }

    public:
      model(
        shared_ptr<model_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      model(
        root const& parent,
        std::string const& id="");

      model(std::string const& id="") : data(new model_data(id)) {}

      model(
        root const& parent,
        model const& other);

      model
      detached_copy() const;

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<root>
      parent() const;

      unsigned
      chains_size() const;

      std::vector<chain> const&
      chains() const;

      std::vector<chain>::iterator
      chains_begin() const { return data->chains.begin(); }

      long
      find_chain_index(
        hierarchy_v2::chain const& chain,
        bool must_be_present=false) const;

      void
      pre_allocate_chains(unsigned number_of_additional_chains);

      void
      new_chains(unsigned number_of_additional_chains);

      void
      insert_chain(long i, hierarchy_v2::chain& chain);

      void
      append_chain(hierarchy_v2::chain& chain);

      void
      remove_chain(long i);

      void
      remove_chain(hierarchy_v2::chain& chain);

      chain
      new_chain(std::string const& id="")
      {
        data->chains.push_back(chain(*this, id));
        return data->chains.back();
      }

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  //! Root attributes.
  class root
  {
    public:
      shared_ptr<root_data> data;

    protected:
      friend class model; // to support model::parent()

      root(shared_ptr<root_data> const& data_, bool) : data(data_) {}

    public:
      root(
        shared_ptr<root_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      root() : data(new root_data) {}

      root
      deep_copy() const;

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      unsigned
      models_size() const;

      std::vector<model> const&
      models() const;

      std::vector<model>::iterator
      models_begin() const { return data->models.begin(); }

      long
      find_model_index(
        hierarchy_v2::model const& model,
        bool must_be_present=false) const;

      void
      pre_allocate_models(unsigned number_of_additional_models);

      void
      new_models(unsigned number_of_additional_models);

      void
      insert_model(long i, hierarchy_v2::model& model);

      void
      append_model(hierarchy_v2::model& model);

      void
      remove_model(long i);

      void
      remove_model(hierarchy_v2::model& model);

      model
      new_model(std::string const& id="")
      {
        data->models.push_back(model(*this, id));
        return data->models.back();
      }

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  inline
  model_data::model_data(
    weak_ptr<root_data> const& parent_,
    std::string const& id_)
  :
    parent(parent_),
    id(id_)
  {}

  inline
  model_data::model_data(std::string const& id_) : id(id_) {}

  inline
  model_data::model_data(
    weak_ptr<root_data> const& parent_,
    model_data const& other)
  :
    parent(parent_),
    id(other.id)
  {}

  inline
  model::model(
    root const& parent,
    std::string const& id)
  :
    data(new model_data(parent.data, id))
  {}

  inline
  chain_data::chain_data(
    weak_ptr<model_data> const& parent_,
    std::string const& id_)
  :
    parent(parent_),
    id(id_)
  {}

  inline
  chain_data::chain_data(
    std::string const& id_)
  :
    id(id_)
  {}

  inline
  chain_data::chain_data(
    weak_ptr<model_data> const& parent_,
    chain_data const& other)
  :
    parent(parent_),
    id(other.id)
  {}

  inline
  chain::chain(
    model const& parent,
    std::string const& id)
  :
    data(new chain_data(parent.data, id))
  {}

  inline
  residue_group_data::residue_group_data(
    weak_ptr<chain_data> const& parent_,
    const char* resseq_,
    const char* icode_,
    bool link_to_previous_)
  :
    parent(parent_),
    resseq(resseq_),
    icode(icode_),
    link_to_previous(link_to_previous_)
  {}

  inline
  residue_group_data::residue_group_data(
    const char* resseq_,
    const char* icode_,
    bool link_to_previous_)
  :
    resseq(resseq_),
    icode(icode_),
    link_to_previous(link_to_previous_)
  {}

  inline
  residue_group_data::residue_group_data(
    weak_ptr<chain_data> const& parent_,
    residue_group_data const& other)
  :
    parent(parent_),
    resseq(other.resseq),
    icode(other.icode),
    link_to_previous(other.link_to_previous)
  {}

  inline
  residue_group::residue_group(
    chain const& parent,
    const char* resseq,
    const char* icode,
    bool link_to_previous)
  :
    data(new residue_group_data(parent.data, resseq, icode, link_to_previous))
  {}

  inline
  atom_group_data::atom_group_data(
    weak_ptr<residue_group_data> const& parent_,
    const char* altloc_,
    const char* resname_)
  :
    parent(parent_),
    altloc(altloc_),
    resname(resname_)
  {}

  inline
  atom_group_data::atom_group_data(
    const char* altloc_,
    const char* resname_)
  :
    altloc(altloc_),
    resname(resname_)
  {}

  inline
  atom_group_data::atom_group_data(
    weak_ptr<residue_group_data> const& parent_,
    atom_group_data const& other)
  :
    parent(parent_),
    altloc(other.altloc),
    resname(other.resname)
  {}

  inline
  atom_group::atom_group(
    residue_group const& parent,
    const char* altloc,
    const char* resname)
  :
    data(new atom_group_data(parent.data, altloc, resname))
  {}

  inline
  atom::atom(
    atom_group const& parent,
    const char* name, const char* segid,
    const char* element, const char* charge, const char* serial,
    vec3 const& xyz, vec3 const& sigxyz,
    double occ, double sigocc,
    double b, double sigb,
    sym_mat3 const& uij,
    sym_mat3 const& siguij,
    bool hetero)
  :
    data(new atom_data(
      parent.data,
      name, segid,
      element, charge, serial,
      xyz, sigxyz,
      occ, sigocc,
      b, sigb,
      uij, siguij,
      hetero))
  {}

  inline
  atom::atom(
    atom_group const& parent,
    atom const& other)
  :
    data(new atom_data(parent.data, *other.data))
  {}

}}} // namespace iotbx::pdb::hierarchy_v2

#endif // IOTBX_PDB_HIERARCHY_V2_H
