#ifndef IOTBX_PDB_HIERARCHY_H
#define IOTBX_PDB_HIERARCHY_H

#include <iotbx/pdb/small_str.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/sym_mat3.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <boost/cstdint.hpp>
#include <vector>
#include <string>

namespace iotbx {

//! Handling of files in PDB format.
namespace pdb {

  namespace af = scitbx::af;

  typedef scitbx::vec3<double> vec3;
  typedef scitbx::sym_mat3<double> sym_mat3;
  typedef boost::int32_t int32_t;

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
  };

  //! Holder for conformer attributes (to be held by a shared_ptr).
  class conformer_data : boost::noncopyable
  {
    protected:
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
  };

  //! Holder for residue attributes (to be held by a shared_ptr).
  class residue_data : boost::noncopyable
  {
    protected:
      friend class residue;
      weak_ptr<conformer_data> parent;
    public:
      str4 name;
      int32_t seq;
      str4 icode;
    protected:
      std::vector<atom> atoms;

      inline
      residue_data(
        weak_ptr<conformer_data> const& parent,
        const char* name,
        int32_t seq,
        const char* const& icode);

      inline
      residue_data(
        const char* name_,
        int32_t seq_,
        const char* icode_);

    public:
      std::string
      id() const
      {
        char buf[32];
        std::sprintf(buf, "%-4s%4d%s", name.elems, seq, icode.elems);
        return std::string(buf);
      }
  };

  //! Holder for atom attributes (to be held by a shared_ptr).
  class atom_data : boost::noncopyable
  {
    protected:
      friend class atom;
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
        hetero(hetero_)
      {}
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
        sym_mat3 const& uij=sym_mat3(0,0,0,0,0,0),
        sym_mat3 const& siguij=sym_mat3(0,0,0,0,0,0),
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

      void
      pre_allocate_parents(unsigned number_of_additional_parents)
      {
        data->parents.reserve(
          data->parents.size()+number_of_additional_parents);
      }

      unsigned
      parents_size() const
      {
        return static_cast<unsigned>(data->parents.size());
      }

      inline
      af::shared<residue>
      parents() const;

      inline
      void
      add_parent(residue const& new_parent);
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
        int32_t seq=0,
        const char* icode="");

      residue(
        const char* name="",
        int32_t seq=0,
        const char* icode="")
      :
        data(new residue_data(name, seq, icode))
      {}

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      inline
      boost::optional<conformer>
      parent() const;

      inline
      void
      set_parent(conformer const& new_parent);

      void
      pre_allocate_atoms(unsigned number_of_additional_atoms)
      {
        data->atoms.reserve(data->atoms.size()+number_of_additional_atoms);
      }

      void
      new_atoms(unsigned number_of_additional_atoms)
      {
        pre_allocate_atoms(number_of_additional_atoms);
        for(unsigned i=0;i<number_of_additional_atoms;i++) {
          atom new_atom;
          new_atom.pre_allocate_parents(1);
          new_atom.add_parent(*this);
          data->atoms.push_back(new_atom);
        }
      }

      void
      add_atom(atom& new_atom)
      {
        new_atom.add_parent(*this);
        data->atoms.push_back(new_atom);
      }

      std::vector<atom> const&
      atoms() const { return data->atoms; }
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

      inline
      boost::optional<chain>
      parent() const;

      inline
      void
      set_parent(chain const& new_parent);

      void
      pre_allocate_residues(unsigned number_of_additional_residues)
      {
        data->residues.reserve(
          data->residues.size()+number_of_additional_residues);
      }

      void
      new_residues(unsigned number_of_additional_residues)
      {
        pre_allocate_residues(number_of_additional_residues);
        for(unsigned i=0;i<number_of_additional_residues;i++) {
          data->residues.push_back(residue(*this));
        }
      }

      residue
      new_residue(
        const char* name="",
        int32_t seq=0,
        const char* icode="")
      {
        data->residues.push_back(residue(*this, name, seq, icode));
        return data->residues.back();
      }

      std::vector<residue> const&
      residues() const { return data->residues; }
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

      inline
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
      new_conformers(unsigned number_of_additional_conformers)
      {
        pre_allocate_conformers(number_of_additional_conformers);
        for(unsigned i=0;i<number_of_additional_conformers;i++) {
          data->conformers.push_back(conformer(*this));
        }
      }

      std::vector<conformer> const&
      conformers() const { return data->conformers; }

      //! Not available from Python.
      conformer
      get_conformer(unsigned i_conformer) const
      {
        return data->conformers[i_conformer];
      }
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

      inline
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
      new_chains(unsigned number_of_additional_chains)
      {
        pre_allocate_chains(number_of_additional_chains);
        for(unsigned i=0;i<number_of_additional_chains;i++) {
          data->chains.push_back(chain(*this));
        }
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
      new_models(unsigned number_of_additional_models)
      {
        pre_allocate_models(number_of_additional_models);
        for(unsigned i=0;i<number_of_additional_models;i++) {
          data->models.push_back(model(*this));
        }
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

      std::map<std::string, unsigned>
      residue_name_counts() const
      {
        std::map<std::string, unsigned> result;
        std::map<str4, unsigned> result_str4;
        for(unsigned jm=0;jm<models().size();jm++) {
          model const& m = models()[jm];
          for(unsigned jc=0;jc<m.chains().size();jc++) {
            chain const& c = m.chains()[jc];
            for(unsigned jf=0;jf<c.conformers().size();jf++) {
              conformer const& f = c.conformers()[jf];
              for(unsigned jr=0;jr<f.residues().size();jr++) {
                residue const& r = f.residues()[jr];
                result_str4[r.data->name]++;
              }
            }
          }
        }
        for(std::map<str4, unsigned>::const_iterator
              i=result_str4.begin();i!=result_str4.end();i++) {
          result[i->first.elems] = i->second;
        }
        return result;
      }
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
  boost::optional<hierarchy>
  model::parent() const
  {
    shared_ptr<hierarchy_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<hierarchy>();
    return boost::optional<hierarchy>(hierarchy(parent, true));
  }

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

  inline
  chain::chain(
    model const& parent,
    std::string const& id)
  :
    data(new chain_data(parent.data, id))
  {}

  inline
  boost::optional<model>
  chain::parent() const
  {
    shared_ptr<model_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<model>();
    return boost::optional<model>(model(parent, true));
  }

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
  boost::optional<chain>
  conformer::parent() const
  {
    shared_ptr<chain_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<chain>();
    return boost::optional<chain>(chain(parent, true));
  }

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
    int32_t seq_,
    const char* const& icode_)
  :
    parent(parent_),
    name(name_),
    seq(seq_),
    icode(icode_)
  {}

  inline
  residue_data::residue_data(
    const char* name_,
    int32_t seq_,
    const char* icode_)
  :
    name(name_),
    seq(seq_),
    icode(icode_)
  {}

  inline
  residue::residue(
    conformer const& parent,
    const char* name,
    int32_t seq,
    const char* icode)
  :
    data(new residue_data(parent.data, name, seq, icode))
  {}

  inline
  boost::optional<conformer>
  residue::parent() const
  {
    shared_ptr<conformer_data> parent = data->parent.lock();
    if (parent.get() == 0) return boost::optional<conformer>();
    return boost::optional<conformer>(conformer(parent, true));
  }

  inline
  void
  residue::set_parent(conformer const& new_parent)
  {
    data->parent = new_parent.data;
  }

  inline
  af::shared<residue>
  atom::parents() const
  {
    af::shared<residue> result;
    unsigned n = parents_size();
    result.reserve(n);
    for(unsigned i=0;i<n;i++) {
      shared_ptr<residue_data> parent = data->parents[i].lock();
      if (parent.get() != 0) result.push_back(residue(parent, true));
    }
    return result;
  }

  inline
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

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_HIERARCHY_H
