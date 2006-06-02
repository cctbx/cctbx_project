#ifndef IOTBX_PDB_HIERARCHY_H
#define IOTBX_PDB_HIERARCHY_H

#include <iotbx/pdb/common_residue_names.h>
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

    public:
      inline
      bool
      has_multiple_conformers() const;
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

    public:
      bool
      parent_chain_has_multiple_conformers() const
      {
        shared_ptr<chain_data> p = parent.lock();
        if (p.get() == 0) return false;
        return p->has_multiple_conformers();
      }
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
      bool link_to_previous;
    protected:
      std::vector<atom> atoms;

      inline
      residue_data(
        weak_ptr<conformer_data> const& parent_,
        const char* name_,
        int32_t seq_,
        const char* const& icode_,
        bool link_to_previous_);

      inline
      residue_data(
        const char* name_,
        int32_t seq_,
        const char* icode_,
        bool link_to_previous_);

      inline
      residue_data(
        weak_ptr<conformer_data> const& parent_,
        residue_data const& other);

    public:
      std::string
      id() const
      {
        char buf[32];
        std::sprintf(buf, "%-4s%4d%s", name.elems, seq, icode.elems);
        return std::string(buf);
      }

      bool
      parent_chain_has_multiple_conformers() const
      {
        shared_ptr<conformer_data> p = parent.lock();
        if (p.get() == 0) return false;
        return p->parent_chain_has_multiple_conformers();
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

      inline
      atom_data(
        weak_ptr<residue_data> const& parent_,
        atom_data const& other);

      bool
      is_alternative() const
      {
        if (parents.size() != 1) return false;
        shared_ptr<residue_data> p0 = parents[0].lock();
        if (p0.get() == 0) return false;
        return p0->parent_chain_has_multiple_conformers();
      }
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

      void
      pre_allocate_parents(unsigned number_of_additional_parents)
      {
        data->parents.reserve(parents_size()+number_of_additional_parents);
      }

      unsigned
      parents_size() const
      {
        unsigned result = 0;
        unsigned n = static_cast<unsigned>(data->parents.size());
        for(unsigned i=0;i<n;i++) {
          shared_ptr<residue_data> parent = data->parents[i].lock();
          if (parent.get() != 0) result++;
        }
        return result;
      }

      inline
      af::shared<residue>
      parents() const;

      inline
      void
      add_parent(residue const& new_parent);

      bool
      is_alternative() const { return data->is_alternative(); }
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
        const char* icode="",
        bool link_to_previous=true);

      residue(
        const char* name="",
        int32_t seq=0,
        const char* icode="",
        bool link_to_previous=true)
      :
        data(new residue_data(name, seq, icode, link_to_previous))
      {}

      inline
      residue(
        conformer const& parent,
        residue const& other);

      std::size_t
      memory_id() const { return reinterpret_cast<std::size_t>(data.get()); }

      inline
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

      unsigned
      atoms_size() const
      {
        return static_cast<unsigned>(data->atoms.size());
      }

      std::vector<atom> const&
      atoms() const { return data->atoms; }

      unsigned
      number_of_alternative_atoms() const
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
      reset_atom_tmp(int new_value) const
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
      center_of_geometry() const
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
      number_of_atoms() const
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
      number_of_alternative_atoms() const
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
      reset_atom_tmp(int new_value) const
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
      residue_centers_of_geometry() const
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
      select_residues_in_place(af::const_ref<std::size_t> const& selection)
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
      select_residues(af::const_ref<std::size_t> const& selection)
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

      //! Not available in Python.
      af::shared<std::size_t>
      residue_class_selection(
        std::set<str4> const& class_set,
        bool negate) const
      {
        unsigned n = residues_size();
        af::shared<std::size_t> result((af::reserve(n)));
        std::set<str4>::const_iterator class_set_end = class_set.end();
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
      residue_class_selection(
        std::string const& class_name,
        bool negate=false) const
      {
        if (class_name == "common_amino_acid")
          return residue_class_selection(
            common_residue_names::amino_acid_set(), negate);
        if (class_name == "common_rna_dna")
          return residue_class_selection(
            common_residue_names::rna_dna_set(), negate);
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
      residue_class_selection(
        af::const_ref<std::string> const& residue_names,
        bool negate=false) const
      {
        std::set<str4> class_set;
        common_residue_names::initialize_set(class_set, residue_names);
        return residue_class_selection(class_set, negate);
      }

      void
      select_residue_class_in_place(
        std::string const& class_name,
        bool negate=false)
      {
        af::shared<std::size_t>
          sel = residue_class_selection(class_name, negate);
        if (sel.size() != residues().size()) {
          select_residues_in_place(sel.const_ref());
        }
      }

      void
      select_residue_class_in_place(
        af::const_ref<std::string> const& residue_names,
        bool negate=false)
      {
        af::shared<std::size_t>
          sel = residue_class_selection(residue_names, negate);
        if (sel.size() != residues().size()) {
          select_residues_in_place(sel.const_ref());
        }
      }
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
      reset_atom_tmp(int new_value) const
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
      extract_sites_cart() const
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
      reset_atom_tmp(int new_value) const
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
      reset_atom_tmp(int new_value) const
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
      overall_counts() const
      {
        reset_atom_tmp(0);
        boost::shared_ptr<overall_counts_holder>
          result(new overall_counts_holder);
        std::map<str4, unsigned> residue_name_counts_str4;
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
                  residue_name_counts_str4[r.data->name]++;
                }
              }
            }
          }
        }
        for(std::map<str4, unsigned>::const_iterator
              i=residue_name_counts_str4.begin();
              i!=residue_name_counts_str4.end();i++) {
          result->n_residues += i->second;
          result->residue_names[i->first.elems] = i->second;
          result->residue_name_classes[
            common_residue_names::get_class(i->first)] += i->second;
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
    int32_t seq_,
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
    int32_t seq,
    const char* icode,
    bool link_to_previous)
  :
    data(new residue_data(parent.data, name, seq, icode, link_to_previous))
  {}

  inline
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
    tmp(0)
  {}

  inline
  atom::atom(
    residue const& parent,
    atom const& other)
  :
    data(new atom_data(parent.data, *other.data))
  {}

  inline
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
