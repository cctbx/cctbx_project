#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/list.hpp>
#include <boost/python/str.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/boost_python/stl_map_as_dict.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <iotbx/pdb/hierarchy_v2.h>

namespace boost_python_meta_ext { struct holder {}; }

namespace iotbx { namespace pdb { namespace hierarchy_v2 {
namespace {

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

#define IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(attr) \
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

    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(name)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(segid)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(element)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(charge)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(serial)

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

    static int
    get_tmp(w_t const& self) { return self.data->tmp; }

    static void
    set_tmp(w_t const& self, int new_tmp) { self.data->tmp = new_tmp; }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("atom", no_init)
        .def(init<>())
        .def(init<atom_group const&, atom const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .def("set_name", set_name, (arg_("new_name")), return_self<>())
        .def("set_segid", set_segid, (arg_("new_segid")), return_self<>())
        .def("set_element", set_element, (arg_("new_element")), return_self<>())
        .def("set_charge", set_charge, (arg_("new_charge")), return_self<>())
        .def("set_serial", set_serial, (arg_("new_serial")), return_self<>())
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
        .add_property("serial",
          make_function(get_serial), make_function(set_serial))
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
        .add_property("tmp",
          make_function(get_tmp), make_function(set_tmp))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<atom, atom_group>::wrapper)
        .def("uij_is_defined", &w_t::uij_is_defined)
        .def("siguij_is_defined", &w_t::siguij_is_defined)
        .def("determine_chemical_element_simple",
          &w_t::determine_chemical_element_simple)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<
          atom, rir>::wrap("af_shared_atom");
      }
    }
  };

#define IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(parent_t, child_t, method) \
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

  struct atom_group_wrappers
  {
    typedef atom_group w_t;

    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(altloc)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(resname)

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(atom_group, atom, atoms)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      find_atom_index_overloads, find_atom_index, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("atom_group", no_init)
        .def(init<
          residue_group const&,
            optional<const char*, const char*> >((
              arg_("parent"), arg_("altloc")="", arg_("resname")="")))
        .def(init<
          optional<const char*, const char*> >((
            arg_("altloc")="", arg_("resname")="")))
        .def(init<residue_group const&, atom_group const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .add_property("altloc",
          make_function(get_altloc), make_function(set_altloc))
        .add_property("resname",
          make_function(get_resname), make_function(set_resname))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<atom_group, residue_group>::wrapper)
        .def("pre_allocate_atoms", &w_t::pre_allocate_atoms,
          (arg_("number_of_additional_atoms")))
        .def("new_atoms", &w_t::new_atoms,
          (arg_("number_of_additional_atoms")))
        .def("insert_atom", &w_t::insert_atom, (arg_("i"), arg_("atom")))
        .def("append_atom", &w_t::append_atom, (arg_("atom")))
        .def("remove_atom",
          (void(w_t::*)(long)) &w_t::remove_atom, (arg_("i")))
        .def("remove_atom",
          (void(w_t::*)(atom&)) &w_t::remove_atom, (arg_("atom")))
        .def("atoms", get_atoms)
        .def("atoms_size", &w_t::atoms_size)
        .def("find_atom_index", &w_t::find_atom_index,
          find_atom_index_overloads((
            arg_("atom"), arg_("must_be_present")=false)))
        .def("reset_atom_tmp", &w_t::reset_atom_tmp, (arg_("new_value")))
      ;
    }
  };

  struct residue_group_wrappers
  {
    typedef residue_group w_t;

    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(resseq)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(icode)

    static bool
    get_link_to_previous(w_t const& self)
    {
      return self.data->link_to_previous;
    }

    static void
    set_link_to_previous(w_t const& self, bool new_link_to_previous)
    {
      self.data->link_to_previous = new_link_to_previous;
    }

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(residue_group, atom_group, atom_groups)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      new_atom_group_overloads, new_atom_group, 0, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("residue_group", no_init)
        .def(init<chain const&, optional<const char*, const char*, bool> >((
          arg_("parent"),
          arg_("resseq")="", arg_("icode")="", arg_("link_to_previous")=true)))
        .def(init<optional<const char*, const char*, bool> >((
          arg_("resseq")="", arg_("icode")="", arg_("link_to_previous")=true)))
        .def(init<chain const&, residue_group const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .add_property("resseq",
          make_function(get_resseq), make_function(set_resseq))
        .add_property("icode",
          make_function(get_icode), make_function(set_icode))
        .add_property("link_to_previous",
          make_function(get_link_to_previous),
          make_function(set_link_to_previous))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue_group, chain>::wrapper)
        .def("pre_allocate_atom_groups", &w_t::pre_allocate_atom_groups, (
          arg_("number_of_additional_atom_groups")))
        .def("new_atom_groups", &w_t::new_atom_groups, (
          arg_("number_of_additional_atom_groups")))
        .def("new_atom_group",
          &w_t::new_atom_group, new_atom_group_overloads((
            arg_("altloc")="", arg_("resname")="")))
        .def("append_atom_group", &w_t::append_atom_group, (
          arg_("atom_group")))
        .def("atom_groups_size", &w_t::atom_groups_size)
        .def("atom_groups", get_atom_groups)
        .def("reset_atom_tmp", &w_t::reset_atom_tmp, (arg_("new_value")))
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

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(chain, residue_group, residue_groups)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      new_residue_group_overloads, new_residue_group, 0, 3)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("chain", no_init)
        .def(init<model const&, optional<std::string const&> >((
          arg_("parent"), arg_("id")="")))
        .def(init<std::string const&>((
          arg_("id")="")))
        .def(init<model const&, chain const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<chain, model>::wrapper)
        .def("pre_allocate_residue_groups", &w_t::pre_allocate_residue_groups,(
          arg_("number_of_additional_residue_groups")))
        .def("new_residue_groups", &w_t::new_residue_groups, (
          arg_("number_of_additional_residue_groups")))
        .def("new_residue_group", &w_t::new_residue_group,
          new_residue_group_overloads((
            arg_("resseq")="",
            arg_("icode")="",
            arg_("link_to_previous")=true)))
        .def("append_residue_group", &w_t::append_residue_group, (
          arg_("residue_group")))
        .def("residue_groups_size", &w_t::residue_groups_size)
        .def("residue_groups", get_residue_groups)
        .def("reset_atom_tmp", &w_t::reset_atom_tmp, (arg_("new_value")))
      ;
    }
  };

  struct model_wrappers
  {
    typedef model w_t;

    static std::string
    get_id(w_t const& self) { return self.data->id; }

    static void
    set_id(w_t const& self, std::string const& new_id)
    {
      self.data->id = new_id;
    }

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(model, chain, chains)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("model", no_init)
        .def(init<root const&, optional<std::string> >((
          arg_("parent"), arg_("id")="")))
        .def(init<std::string>((arg_("id")="")))
        .def(init<root const&, model const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<model, root>::wrapper)
        .def("pre_allocate_chains", &w_t::pre_allocate_chains,
          (arg_("number_of_additional_chains")))
        .def("new_chains", &w_t::new_chains,
          (arg_("number_of_additional_chains")))
        .def("new_chain", &w_t::new_chain, (arg_("id")))
        .def("append_chain", &w_t::append_chain, (arg_("chain")))
        .def("chains_size", &w_t::chains_size)
        .def("chains", get_chains)
        .def("reset_atom_tmp", &w_t::reset_atom_tmp, (arg_("new_value")))
      ;
    }
  };

  struct root_wrappers
  {
    typedef root w_t;

    static af::shared<std::string>
    get_info(w_t const& self) { return self.data->info; }

    static void
    set_info(w_t const& self, af::shared<std::string> const& new_info)
    {
      self.data->info = new_info;
    }

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(root, model, models)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("root", no_init)
        .def(init<>())
        .def("deep_copy", &w_t::deep_copy)
        .add_property("info", make_function(get_info), make_function(set_info))
        .def("memory_id", &w_t::memory_id)
        .def("pre_allocate_models", &w_t::pre_allocate_models,
          (arg_("number_of_additional_models")))
        .def("new_models", &w_t::new_models,
          (arg_("number_of_additional_models")))
        .def("new_model", &w_t::new_model, (arg_("id")))
        .def("append_model", &w_t::append_model, (arg_("model")))
        .def("reset_atom_tmp", &w_t::reset_atom_tmp, (arg_("new_value")))
        .def("models_size", &w_t::models_size)
        .def("models", get_models)
      ;
    }
  };

  void
  wrap_hierarchy_v2()
  {
    atom_wrappers::wrap();
    atom_group_wrappers::wrap();
    residue_group_wrappers::wrap();
    chain_wrappers::wrap();
    model_wrappers::wrap();
    root_wrappers::wrap();
  }

}}}} // namespace iotbx::pdb::hierarchy_v2::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_hierarchy_v2_ext)
{
  iotbx::pdb::hierarchy_v2::wrap_hierarchy_v2();
}
