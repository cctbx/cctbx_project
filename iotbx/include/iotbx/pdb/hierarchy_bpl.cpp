#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/str.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_arg.hpp>
#include <iotbx/pdb/hierarchy.h>

namespace iotbx { namespace pdb {
namespace {

#define IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(attr) \
    static \
    boost::python::str \
    get_##attr(w_t const& self) \
    { \
      return boost::python::str(self.data->attr.elems); \
    } \
    static \
    void \
    set_##attr(w_t& self, const char* value) \
    { \
      if (!self.data->attr.replace_with(value)) { \
        char buf[128]; \
        std::sprintf(buf, \
          "string is too long for " #attr " attribute" \
          " (maximum length is %u characters, %lu given).", \
            self.data->attr.capacity(), \
            static_cast<unsigned long>(std::strlen(value))); \
        PyErr_SetString(PyExc_ValueError, buf); \
        boost::python::throw_error_already_set(); \
      } \
    }

  struct atom_wrappers
  {
    typedef atom w_t;

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(name)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(segid)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(element)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(charge)

    static vec3
    get_xyz(w_t const& self) { return self.data->xyz; }

    static void
    set_xyz(w_t const& self, vec3 const& new_xyz)
    {
      self.data->xyz = new_xyz;
    }

    static vec3
    get_sigxyz(w_t const& self) { return self.data->sigxyz; }

    static void
    set_sigxyz(w_t const& self, vec3 const& new_sigxyz)
    {
      self.data->sigxyz = new_sigxyz;
    }

    static double
    get_occ(w_t const& self) { return self.data->occ; }

    static void
    set_occ(w_t const& self, double new_occ)
    {
      self.data->occ = new_occ;
    }

    static double
    get_sigocc(w_t const& self) { return self.data->sigocc; }

    static void
    set_sigocc(w_t const& self, double new_sigocc)
    {
      self.data->sigocc = new_sigocc;
    }

    static double
    get_b(w_t const& self) { return self.data->b; }

    static void
    set_b(w_t const& self, double new_b)
    {
      self.data->b = new_b;
    }

    static double
    get_sigb(w_t const& self) { return self.data->sigb; }

    static void
    set_sigb(w_t const& self, double new_sigb)
    {
      self.data->sigb = new_sigb;
    }

    static sym_mat3
    get_uij(w_t const& self) { return self.data->uij; }

    static void
    set_uij(w_t const& self, sym_mat3 const& new_uij)
    {
      self.data->uij = new_uij;
    }

    static sym_mat3
    get_siguij(w_t const& self) { return self.data->siguij; }

    static void
    set_siguij(w_t const& self, sym_mat3 const& new_siguij)
    {
      self.data->siguij = new_siguij;
    }

    static bool
    get_hetero(w_t const& self) { return self.data->hetero; }

    static void
    set_hetero(w_t const& self, bool new_hetero)
    {
      self.data->hetero = new_hetero;
    }

    static
    boost::python::list
    parents(w_t const& self)
    {
      boost::python::list result;
      af::shared<residue> parents = self.parents();
      unsigned n = static_cast<unsigned>(parents.size());
      const residue* p = parents.begin();
      for(unsigned i=0;i<n;i++) result.append(*p++);
      return result;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("atom", no_init)
        .def(init<>())
        .def("set_name", set_name, (arg_("new_name")), return_self<>())
        .def("set_segid", set_segid, (arg_("new_segid")), return_self<>())
        .def("set_element", set_element, (arg_("new_element")), return_self<>())
        .def("set_charge", set_charge, (arg_("new_charge")), return_self<>())
        .def("set_xyz", set_xyz, (arg_("new_xyz")), return_self<>())
        .def("set_sigxyz", set_sigxyz, (arg_("new_sigxyz")), return_self<>())
        .def("set_occ", set_occ, (arg_("new_occ")), return_self<>())
        .def("set_sigocc", set_sigocc, (arg_("new_sigocc")), return_self<>())
        .def("set_b", set_b, (arg_("new_b")), return_self<>())
        .def("set_sigb", set_sigb, (arg_("new_sigb")), return_self<>())
        .def("set_uij", set_uij, (arg_("new_uij")), return_self<>())
        .def("set_siguij", set_siguij, (arg_("new_siguij")), return_self<>())
        .def("set_hetero", set_hetero, (arg_("new_hetero")), return_self<>())
        .add_property("name",
          make_function(get_name), make_function(set_name))
        .add_property("segid",
          make_function(get_segid), make_function(set_segid))
        .add_property("element",
          make_function(get_element), make_function(set_element))
        .add_property("charge",
          make_function(get_charge), make_function(set_charge))
        .add_property("xyz",
          make_function(get_xyz), make_function(set_xyz))
        .add_property("sigxyz",
          make_function(get_sigxyz), make_function(set_sigxyz))
        .add_property("occ",
          make_function(get_occ), make_function(set_occ))
        .add_property("sigocc",
          make_function(get_sigocc), make_function(set_sigocc))
        .add_property("b",
          make_function(get_b), make_function(set_b))
        .add_property("sigb",
          make_function(get_sigb), make_function(set_sigb))
        .add_property("uij",
          make_function(get_uij), make_function(set_uij))
        .add_property("siguij",
          make_function(get_siguij), make_function(set_siguij))
        .add_property("hetero",
          make_function(get_hetero), make_function(set_hetero))
        .def("memory_id", &w_t::memory_id)
        .def("pre_allocate_parents", &w_t::pre_allocate_parents, (
          arg_("number_of_additional_parents")))
        .def("parents_size", &w_t::parents_size)
        .def("parents", parents)
        .def("add_parent", &w_t::add_parent, (arg_("new_parent")))
      ;
    }
  };

  template <typename ChildType, typename ParentType>
  struct get_parent
  {
    static
    boost::python::object
    wrapper(ChildType const& child)
    {
      boost::optional<ParentType> parent = child.parent();
      if (!parent) return boost::python::object();
      return boost::python::object(*parent);
    }
  };

#define IOTBX_PDB_HIERARCHY_GET_CHILDREN(parent_t, child_t, method) \
  static \
  boost::python::list \
  get_##method(parent_t const& parent) \
  { \
    boost::python::list result; \
    std::vector<child_t> const& children = parent.method(); \
    unsigned n = static_cast<unsigned>(children.size()); \
    for(unsigned i=0;i<n;i++) result.append(children[i]); \
    return result; \
  }

  struct residue_wrappers
  {
    typedef residue w_t;

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(name)

    static int
    get_seq(w_t const& self) { return self.data->seq; }

    static void
    set_seq(w_t const& self, int new_seq)
    {
      self.data->seq = new_seq;
    }

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(icode)

    static std::string
    id(w_t const& self) { return self.data->id(); }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(residue, atom, atoms)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("residue", no_init)
        .def(init<
          conformer const&,
            optional<const char*, int32_t, const char*> >((
              arg_("parent"),
              arg_("name")="", arg_("seq")=0, arg_("icode")="")))
        .def(init<
          optional<const char*, int32_t, const char*> >((
            arg_("name")="", arg_("seq")=0, arg_("icode")="")))
        .add_property("name",
          make_function(get_name), make_function(set_name))
        .add_property("seq",
          make_function(get_seq), make_function(set_seq))
        .add_property("icode",
          make_function(get_icode), make_function(set_icode))
        .def("id", id)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue, conformer>::wrapper)
        .def("set_parent", &w_t::set_parent, (arg_("new_parent")))
        .def("pre_allocate_atoms", &w_t::pre_allocate_atoms,
          (arg_("number_of_additional_atoms")))
        .def("new_atoms", &w_t::new_atoms,
          (arg_("number_of_additional_atoms")))
        .def("add_atom", &w_t::add_atom, (arg_("new_atom")))
        .def("atoms", get_atoms)
      ;
    }
  };

  struct conformer_wrappers
  {
    typedef conformer w_t;

    static std::string
    get_id(w_t const& self) { return self.data->id; }

    static void
    set_id(w_t const& self, std::string const& new_id)
    {
      self.data->id = new_id;
    }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(conformer, residue, residues)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("conformer", no_init)
        .def(init<chain const&, optional<std::string const&> >((
          arg_("parent"), arg_("id")="")))
        .def(init<std::string const&>((
          arg_("id")="")))
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<conformer, chain>::wrapper)
        .def("set_parent", &w_t::set_parent, (arg_("new_parent")))
        .def("pre_allocate_residues", &w_t::pre_allocate_residues, (
          arg_("number_of_additional_residues")))
        .def("new_residues", &w_t::new_residues, (
          arg_("number_of_additional_residues")))
        .def("new_residue", &w_t::new_residue, (
          arg_("name")="", arg_("seq")=0, arg_("icode")=""))
        .def("residues", get_residues)
      ;
    }
  };

  struct chain_wrappers
  {
    typedef chain w_t;

    static std::string
    get_id(w_t const& self) { return self.data->id; }

    static void
    set_id(w_t const& self, std::string const& new_id)
    {
      self.data->id = new_id;
    }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(chain, conformer, conformers)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("chain", no_init)
        .def(init<model const&, optional<std::string const&> >((
          arg_("parent"), arg_("id")="")))
        .def(init<std::string const&>((
          arg_("id")="")))
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<chain, model>::wrapper)
        .def("set_parent", &w_t::set_parent, (arg_("new_parent")))
        .def("pre_allocate_conformers", &w_t::pre_allocate_conformers, (
          arg_("number_of_additional_conformers")))
        .def("new_conformers", &w_t::new_conformers, (
          arg_("number_of_additional_conformers")))
        .def("conformers", get_conformers)
      ;
    }
  };

  struct model_wrappers
  {
    typedef model w_t;

    static int
    get_id(w_t const& self) { return self.data->id; }

    static void
    set_id(w_t const& self, int new_id) { self.data->id = new_id; }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(model, chain, chains)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("model", no_init)
        .def(init<hierarchy const&, optional<int> >((
          arg_("parent"), arg_("id")=0)))
        .def(init<int>((arg_("id")=0)))
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<model, hierarchy>::wrapper)
        .def("set_parent", &w_t::set_parent, (arg_("new_parent")))
        .def("pre_allocate_chains", &w_t::pre_allocate_chains,
          (arg_("number_of_additional_chains")))
        .def("new_chains", &w_t::new_chains,
          (arg_("number_of_additional_chains")))
        .def("chains", get_chains)
        .def("new_chain", &w_t::new_chain, (arg_("chain_id")))
        .def("adopt_chain", &w_t::adopt_chain, (arg_("new_chain")))
      ;
    }
  };

  struct hierarchy_wrappers
  {
    typedef hierarchy w_t;

    static af::shared<std::string>
    get_info(w_t const& self) { return self.data->info; }

    static void
    set_info(w_t const& self, af::shared<std::string> const& new_info)
    {
      self.data->info = new_info;
    }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(hierarchy, model, models)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("hierarchy", no_init)
        .def(init<>())
        .add_property("info", make_function(get_info), make_function(set_info))
        .def("memory_id", &w_t::memory_id)
        .def("pre_allocate_models", &w_t::pre_allocate_models,
          (arg_("number_of_additional_models")))
        .def("new_models", &w_t::new_models,
          (arg_("number_of_additional_models")))
        .def("models", get_models)
        .def("new_model", &w_t::new_model, (arg_("model_id")))
        .def("adopt_model", &w_t::adopt_model, (arg_("new_model")))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    atom_wrappers::wrap();
    residue_wrappers::wrap();
    conformer_wrappers::wrap();
    chain_wrappers::wrap();
    model_wrappers::wrap();
    hierarchy_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_hierarchy() { wrap_all(); }

}}} // namespace iotbx::pdb::boost_python
