#ifndef IOTBX_PDB_HIERARCHY_H
#define IOTBX_PDB_HIERARCHY_H

#include <iotbx/pdb/small_str.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/sym_mat3.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <map>
#include <set>
#include <vector>
#include <string>

namespace iotbx {

//! Handling of files in PDB format.
namespace pdb {

  namespace af = scitbx::af;

  typedef scitbx::vec3<double> vec3;
  typedef scitbx::sym_mat3<double> sym_mat3;

  using boost::shared_ptr;
  using boost::weak_ptr;

  class hierarchy_data;
  class hierarchy;
  class model_data;
  class model;
  class chain_data;
  class chain;
  class conformer_data;
  class conformer;
  class residue_data;
  class residue;
  class atom;

  class combine_conformers;

  //! Holder for hierarchy attributes (to be held by a shared_ptr).
  class hierarchy_data : boost::noncopyable
  {
    protected:
      friend class hierarchy;
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
      weak_ptr<hierarchy_data> parent;
    public:
      int id;
    protected:
      std::vector<chain> chains;

      inline
      model_data(
        weak_ptr<hierarchy_data> const& parent,
        int id);

      inline
      model_data(
        int id);
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
      std::vector<conformer> conformers;

      inline
      chain_data(
        weak_ptr<model_data> const& parent,
        std::string const& id);

      inline
      chain_data(
        std::string const& id);

    public:
      inline
      bool
      has_multiple_conformers() const;
  };

  //! Holder for conformer attributes (to be held by a shared_ptr).
  class conformer_data : boost::noncopyable
  {
    protected:
      friend class atom;
      friend class conformer;
      weak_ptr<chain_data> parent;
    public:
      std::string id;
    protected:
      std::vector<residue> residues;

      inline
      conformer_data(
        weak_ptr<chain_data> const& parent,
        std::string const& id);

      inline
      conformer_data(
        std::string const& id);

    public:
      bool
      parent_chain_has_multiple_conformers() const;
  };

  //! Holder for residue attributes (to be held by a shared_ptr).
  class residue_data : boost::noncopyable
  {
    protected:
      friend class atom;
      friend class residue;
      weak_ptr<conformer_data> parent;
    public:
      str3 name;
      str4 seq;
      str1 icode;
      bool link_to_previous;
    protected:
      std::vector<atom> atoms;

      inline
      residue_data(
        weak_ptr<conformer_data> const& parent_,
        const char* name_,
        const char* seq_,
        const char* const& icode_,
        bool link_to_previous_);

      inline
      residue_data(
        const char* name_,
        const char* seq_,
        const char* icode_,
        bool link_to_previous_);

      inline
      residue_data(
        weak_ptr<conformer_data> const& parent_,
        residue_data const& other);

    public:
      std::string
      id() const;

      bool
      parent_chain_has_multiple_conformers() const;
  };

  //! Holder for atom attributes (to be held by a shared_ptr).
  class atom_data : boost::noncopyable
  {
    protected:
      friend class atom;
      friend class combine_conformers;
      std::vector<weak_ptr<residue_data> > parents;
    public:
      str4 name;
      str4 segid;
      str2 element;
      str2 charge;
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
        str4 name_, str4 segid_, str2 element_, str2 charge_,
        vec3 const& xyz_, vec3 const& sigxyz_,
        double occ_, double sigocc_,
        double b_, double sigb_,
        sym_mat3 const& uij_, sym_mat3 const& siguij_,
        bool hetero_)
      :
        name(name_), segid(segid_), element(element_), charge(charge_),
        xyz(xyz_), sigxyz(sigxyz_),
        occ(occ_), sigocc(sigocc_),
        b(b_), sigb(sigb_),
        uij(uij_), siguij(siguij_),
        hetero(hetero_),
        tmp(0)
      {}

      atom_data(
        weak_ptr<residue_data> const& parent_,
        atom_data const& other);

      bool
      is_alternative() const;
  };

  //! Atom attributes.
  class atom
  {
    public:
      shared_ptr<atom_data> data;

      atom(shared_ptr<atom_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      atom(
        str4 name, str4 segid,
        str2 element, str2 charge,
        vec3 const& xyz, vec3 const& sigxyz,
        double occ, double sigocc,
        double b, double sigb,
        sym_mat3 const& uij, sym_mat3 const& siguij,
        bool hetero)
      :
        data(new atom_data(
          name, segid, element, charge,
          xyz, sigxyz,
          occ, sigocc,
          b, sigb,
          uij, siguij,
          hetero))
      {}

      atom(
        const char* name="", const char* segid="",
        const char* element="", const char* charge="",
        vec3 const& xyz=vec3(0,0,0), vec3 const& sigxyz=vec3(0,0,0),
        double occ=0, double sigocc=0,
        double b=0, double sigb=0,
        sym_mat3 const& uij=sym_mat3(-1,-1,-1,-1,-1,-1),
        sym_mat3 const& siguij=sym_mat3(-1,-1,-1,-1,-1,-1),
        bool hetero=false)
      :
        data(new atom_data(
          name, segid, element, charge,
          xyz, sigxyz,
          occ, sigocc,
          b, sigb,
          uij, siguij,
          hetero))
      {}

      inline
      atom(
        residue const& parent,
        atom const& other);

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

      //! Not available in Python.
      unsigned
      reset_parents();

      void
      pre_allocate_parents(unsigned number_of_additional_parents)
      {
        data->parents.reserve(parents_size()+number_of_additional_parents);
      }

      unsigned
      parents_size() const;

      af::shared<residue>
      parents() const;

      void
      add_parent(residue const& new_parent);

      bool
      is_alternative() const { return data->is_alternative(); }

      boost::optional<std::string>
      determine_chemical_element_simple() const;

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

      //! Not available in Python.
      /*! result must point to an array of size 81 (or greater).
          On return, result is null-terminated.
       */
      unsigned
      format_atom_record(
        char* result,
        int serial,
        const char* resname,
        const char* resseq,
        const char* icode,
        const char* altloc,
        const char* chain_id) const;

      //! Not available in Python.
      /*! result must point to an array of size 81 (or greater).
          On return, result is null-terminated.
       */
      unsigned
      format_atom_record_using_parents(
        char* result,
        int serial) const;
  };

  //! Residue attributes.
  class residue
  {
    public:
      shared_ptr<residue_data> data;

    protected:
      friend class atom;
      residue(shared_ptr<residue_data> const& data_, bool) : data(data_) {}

    public:
      residue(shared_ptr<residue_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      residue(
        conformer const& parent,
        const char* name="",
        const char* seq="",
        const char* icode="",
        bool link_to_previous=true);

      residue(
        const char* name="",
        const char* seq="",
        const char* icode="",
        bool link_to_previous=true)
      :
        data(new residue_data(name, seq, icode, link_to_previous))
      {}

      residue(
        conformer const& parent,
        residue const& other);

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<conformer>
      parent() const;

      inline
      void
      set_parent(conformer const& new_parent);

      inline
      bool
      parent_chain_has_multiple_conformers() const
      {
        return data->parent_chain_has_multiple_conformers();
      }

      void
      pre_allocate_atoms(unsigned number_of_additional_atoms)
      {
        data->atoms.reserve(data->atoms.size()+number_of_additional_atoms);
      }

      void
      new_atoms(unsigned number_of_additional_atoms);

      void
      add_atom(atom& new_atom)
      {
        new_atom.add_parent(*this);
        data->atoms.push_back(new_atom);
      }

      unsigned
      atoms_size() const
      {
        return static_cast<unsigned>(data->atoms.size());
      }

      std::vector<atom> const&
      atoms() const { return data->atoms; }

      af::shared<std::string>
      atom_names() const;

      unsigned
      number_of_alternative_atoms() const;

      unsigned
      reset_atom_tmp(int new_value) const;

      vec3
      center_of_geometry() const;
  };

  //! Conformer attributes.
  class conformer
  {
    public:
      shared_ptr<conformer_data> data;

    protected:
      friend class residue;
      conformer(shared_ptr<conformer_data> const& data_, bool) : data(data_) {}

    public:
      conformer(shared_ptr<conformer_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      inline
      conformer(
        chain const& parent,
        std::string const& id="");

      inline
      conformer(
        std::string const& id="")
      :
        data(new conformer_data(id))
      {}

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<chain>
      parent() const;

      inline
      void
      set_parent(chain const& new_parent);

      bool
      parent_chain_has_multiple_conformers() const
      {
        return data->parent_chain_has_multiple_conformers();
      }

      void
      pre_allocate_residues(unsigned number_of_additional_residues)
      {
        data->residues.reserve(
          data->residues.size()+number_of_additional_residues);
      }

      void
      new_residues(unsigned number_of_additional_residues);

      residue
      new_residue(
        const char* name="",
        const char* seq="",
        const char* icode="",
        bool link_to_previous=true)
      {
        data->residues.push_back(
          residue(*this, name, seq, icode, link_to_previous));
        return data->residues.back();
      }

      unsigned
      residues_size() const
      {
        return static_cast<unsigned>(data->residues.size());
      }

      std::vector<residue> const&
      residues() const { return data->residues; }

      unsigned
      number_of_atoms() const;

      unsigned
      number_of_alternative_atoms() const;

      unsigned
      reset_atom_tmp(int new_value) const;

      af::shared<vec3>
      residue_centers_of_geometry() const;

      void
      select_residues_in_place(
        af::const_ref<std::size_t> const& selection);

      conformer
      select_residues(
        af::const_ref<std::size_t> const& selection);

      //! Not available in Python.
      af::shared<std::size_t>
      residue_class_selection(
        std::set<str3> const& class_set,
        bool negate) const;

      af::shared<std::size_t>
      residue_class_selection(
        std::string const& class_name,
        bool negate=false) const;

      af::shared<std::size_t>
      residue_class_selection(
        af::const_ref<std::string> const& residue_names,
        bool negate=false) const;

      void
      select_residue_class_in_place(
        std::string const& class_name,
        bool negate=false);

      void
      select_residue_class_in_place(
        af::const_ref<std::string> const& residue_names,
        bool negate=false);
  };

  //! Not available from Python.
  struct combine_conformers
  {
    virtual
    void
    process_next_residue(
      std::string const& conformer_id,
      bool is_first_of_group,
      const pdb::residue& residue,
      bool suppress_chain_break) = 0;

    void
    process_single_conformer(pdb::conformer const& conformer);

    void
    process_conformers(std::vector<conformer> const& conformers);
  };

  //! Chain attributes.
  class chain
  {
    public:
      shared_ptr<chain_data> data;

    protected:
      friend class conformer;
      chain(shared_ptr<chain_data> const& data_, bool) : data(data_) {}

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

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<model>
      parent() const;

      inline
      void
      set_parent(model const& new_parent);

      void
      pre_allocate_conformers(unsigned number_of_additional_conformers)
      {
        data->conformers.reserve(
          data->conformers.size()+number_of_additional_conformers);
      }

      void
      new_conformers(unsigned number_of_additional_conformers);

      unsigned
      conformers_size() const
      {
        return static_cast<unsigned>(data->conformers.size());
      }

      std::vector<conformer> const&
      conformers() const { return data->conformers; }

      //! Not available from Python.
      conformer
      get_conformer(unsigned i_conformer) const
      {
        return data->conformers[i_conformer];
      }

      bool
      has_multiple_conformers() const
      {
        return data->has_multiple_conformers();
      }

      unsigned
      number_of_atoms() const;

      unsigned
      reset_atom_tmp(int new_value) const;

      af::shared<vec3>
      extract_sites_cart() const;
  };

  //! Model attributes.
  class model
  {
    public:
      shared_ptr<model_data> data;

    protected:
      friend class chain;
      model(shared_ptr<model_data> const& data_, bool) : data(data_) {}

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
        hierarchy const& parent,
        int id=0);

      model(int id=0) : data(new model_data(id)) {}

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      boost::optional<hierarchy>
      parent() const;

      inline
      void
      set_parent(hierarchy const& new_parent);

      void
      pre_allocate_chains(unsigned number_of_additional_chains)
      {
        data->chains.reserve(
          data->chains.size()+number_of_additional_chains);
      }

      void
      new_chains(unsigned number_of_additional_chains);

      unsigned
      chains_size() const
      {
        return static_cast<unsigned>(data->chains.size());
      }

      std::vector<chain> const&
      chains() const { return data->chains; }

      //! Not available from Python.
      chain
      get_chain(unsigned i_chain) const
      {
        return data->chains[i_chain];
      }

      chain
      new_chain(std::string const& chain_id)
      {
        data->chains.push_back(chain(*this, chain_id));
        return data->chains.back();
      }

      void
      adopt_chain(chain& new_chain)
      {
        new_chain.set_parent(*this);
        data->chains.push_back(new_chain);
      }

      unsigned
      reset_atom_tmp(int new_value) const;
  };

  //! Hierarchy attributes.
  class hierarchy
  {
    public:
      shared_ptr<hierarchy_data> data;

    protected:
      friend class model;
      hierarchy(shared_ptr<hierarchy_data> const& data_, bool) : data(data_) {}

    public:
      hierarchy(
        shared_ptr<hierarchy_data> const& data_)
      :
        data(data_)
      {
        SCITBX_ASSERT(data.get() != 0);
      }

      hierarchy() : data(new hierarchy_data) {}

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      void
      pre_allocate_models(unsigned number_of_additional_models)
      {
        data->models.reserve(
          data->models.size()+number_of_additional_models);
      }

      void
      new_models(unsigned number_of_additional_models);

      unsigned
      models_size() const
      {
        return static_cast<unsigned>(data->models.size());
      }

      std::vector<model> const&
      models() const { return data->models; }

      model
      new_model(int model_id)
      {
        data->models.push_back(model(*this, model_id));
        return data->models.back();
      }

      void
      adopt_model(model& new_model)
      {
        new_model.set_parent(*this);
        data->models.push_back(new_model);
      }

      unsigned
      reset_atom_tmp(int new_value) const;

      struct overall_counts_holder
      {
        unsigned n_models;
        unsigned n_chains;
        std::map<std::string, unsigned> chain_ids;
        unsigned n_conformers;
        std::map<std::string, unsigned> conformer_ids;
        unsigned n_residues;
        std::map<std::string, unsigned> residue_names;
        unsigned n_atoms;
        std::map<std::string, unsigned> residue_name_classes;

        overall_counts_holder()
        :
          n_models(0),
          n_chains(0),
          n_conformers(0),
          n_residues(0),
          n_atoms(0)
        {}
      };

      boost::shared_ptr<overall_counts_holder>
      overall_counts() const;
  };

  inline
  model_data::model_data(
    weak_ptr<hierarchy_data> const& parent_,
    int id_)
  :
    parent(parent_),
    id(id_)
  {}

  inline
  model_data::model_data(int id_) : id(id_) {}

  inline
  model::model(
    hierarchy const& parent,
    int id)
  :
    data(new model_data(parent.data, id))
  {}

  inline
  void
  model::set_parent(
    hierarchy const& new_parent)
  {
    data->parent = new_parent.data;
  }

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

  bool
  chain_data::has_multiple_conformers() const
  {
    return (conformers.size() > 1);
  }

  inline
  chain::chain(
    model const& parent,
    std::string const& id)
  :
    data(new chain_data(parent.data, id))
  {}

  inline
  void
  chain::set_parent(
    model const& new_parent)
  {
    data->parent = new_parent.data;
  }

  inline
  conformer_data::conformer_data(
    weak_ptr<chain_data> const& parent_,
    std::string const& id_)
  :
    parent(parent_),
    id(id_)
  {}

  inline
  conformer_data::conformer_data(
    std::string const& id_)
  :
    id(id_)
  {}

  inline
  conformer::conformer(
    chain const& parent,
    std::string const& id)
  :
    data(new conformer_data(parent.data, id))
  {}

  inline
  void
  conformer::set_parent(
    chain const& new_parent)
  {
    data->parent = new_parent.data;
  }

  inline
  residue_data::residue_data(
    weak_ptr<conformer_data> const& parent_,
    const char* name_,
    const char* seq_,
    const char* const& icode_,
    bool link_to_previous_)
  :
    parent(parent_),
    name(name_),
    seq(seq_),
    icode(icode_),
    link_to_previous(link_to_previous_)
  {}

  inline
  residue_data::residue_data(
    const char* name_,
    const char* seq_,
    const char* icode_,
    bool link_to_previous_)
  :
    name(name_),
    seq(seq_),
    icode(icode_),
    link_to_previous(link_to_previous_)
  {}

  inline
  residue_data::residue_data(
    weak_ptr<conformer_data> const& parent_,
    residue_data const& other)
  :
    parent(parent_),
    name(other.name),
    seq(other.seq),
    icode(other.icode),
    link_to_previous(other.link_to_previous)
  {}

  inline
  residue::residue(
    conformer const& parent,
    const char* name,
    const char* seq,
    const char* icode,
    bool link_to_previous)
  :
    data(new residue_data(parent.data, name, seq, icode, link_to_previous))
  {}

  inline
  void
  residue::set_parent(conformer const& new_parent)
  {
    data->parent = new_parent.data;
  }

  inline
  atom::atom(
    residue const& parent,
    atom const& other)
  :
    data(new atom_data(parent.data, *other.data))
  {}

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_HIERARCHY_H
