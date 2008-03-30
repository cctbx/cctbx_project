#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/boost_python/stl_map_as_dict.h>
#include <scitbx/boost_python/array_as_list.h>
#include <iotbx/pdb/hierarchy_v2.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 {

namespace atoms { void bpl_wrap(); }

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

#define IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET(attr) \
    static \
    boost::python::str \
    get_##attr(w_t const& self) \
    { \
      return boost::python::str(self.data->attr.elems); \
    }

#define IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_SET(attr) \
    static \
    void \
    set_##attr(w_t& self, const char* value) \
    { \
      self.data->attr.replace_with(value); \
    }

#define IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(attr) \
  IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET(attr) \
  IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_SET(attr)

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

    static
    boost::python::object
    format_atom_record(
      w_t const& self,
      const char* replace_floats_with=0)
    {
      boost::python::handle<> str_hdl(PyString_FromStringAndSize(0, 81));
      PyObject* str_obj = str_hdl.get();
      char* str_begin = PyString_AS_STRING(str_obj);
      unsigned str_len = self.format_atom_record(
        str_begin, 0, replace_floats_with);
      str_hdl.release();
      if (_PyString_Resize(&str_obj, static_cast<int>(str_len)) != 0) {
        boost::python::throw_error_already_set();
      }
      return boost::python::object(boost::python::handle<>(str_obj));
    }

#define IOTBX_LOC(R) \
    static \
    boost::python::object \
    format_##R##_record( \
      w_t const& self) \
    { \
      boost::python::handle<> str_hdl(PyString_FromStringAndSize(0, 81)); \
      PyObject* str_obj = str_hdl.get(); \
      char* str_begin = PyString_AS_STRING(str_obj); \
      unsigned str_len = self.format_##R##_record(str_begin, 0); \
      str_hdl.release(); \
      if (_PyString_Resize(&str_obj, static_cast<int>(str_len)) != 0) { \
        boost::python::throw_error_already_set(); \
      } \
      return boost::python::object(boost::python::handle<>(str_obj)); \
    }

    IOTBX_LOC(sigatm)
    IOTBX_LOC(anisou)
    IOTBX_LOC(siguij)

#undef IOTBX_LOC

    BOOST_PYTHON_FUNCTION_OVERLOADS(
      format_atom_record_overloads, format_atom_record, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
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
        .def("sizeof_data", &w_t::sizeof_data)
        .staticmethod("sizeof_data")
        .def("parent", get_parent<atom, atom_group>::wrapper)
        .def("uij_is_defined", &w_t::uij_is_defined)
        .def("siguij_is_defined", &w_t::siguij_is_defined)
        .def("pdb_label_columns", &w_t::pdb_label_columns)
        .def("pdb_element_charge_columns", &w_t::pdb_element_charge_columns)
        .def("format_atom_record", format_atom_record,
          format_atom_record_overloads((
            arg_("self"),
            arg_("replace_floats_with")=0)))
        .def("format_sigatm_record", format_sigatm_record)
        .def("format_anisou_record", format_anisou_record)
        .def("format_siguij_record", format_siguij_record)
        .def("element_is_hydrogen", &w_t::element_is_hydrogen)
        .def("determine_chemical_element_simple",
          &w_t::determine_chemical_element_simple)
      ;
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

#define IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(C) \
        .def(#C "s", get_##C##s) \
        .def(#C "s_size", &w_t::C##s_size) \
        .def("find_" #C "_index", &w_t::find_##C##_index, \
          find_##C##_index_overloads(( \
            arg_(#C), arg_("must_be_present")=false))) \
        .def("pre_allocate_" #C "s", &w_t::pre_allocate_##C##s, \
          (arg_("number_of_additional_" #C "s"))) \
        .def("insert_" #C, &w_t::insert_##C, (arg_("i"), arg_(#C))) \
        .def("append_" #C, &w_t::append_##C, (arg_(#C))) \
        .def("remove_" #C, \
          (void(w_t::*)(long)) &w_t::remove_##C, (arg_("i"))) \
        .def("remove_" #C, \
          (void(w_t::*)(C&)) &w_t::remove_##C, (arg_(#C)))

  template <typename ElementType>
  af::shared<ElementType>
  std_vector_as_af_shared(
    std::vector<ElementType> const& v)
  {
    if (v.size() == 0) return af::shared<ElementType>();
    return af::shared<ElementType>(&*v.begin(), &*v.end());
  }

  struct atom_group_wrappers
  {
    typedef atom_group w_t;

    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(altloc)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET_SET(resname)

    static
    af::shared<atom>
    get_atoms(w_t const& self)
    {
      return std_vector_as_af_shared(self.atoms());
    }

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
        .add_property("altloc",
          make_function(get_altloc), make_function(set_altloc))
        .add_property("resname",
          make_function(get_resname), make_function(set_resname))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<atom_group, residue_group>::wrapper)
        IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(atom)
        .def("confid", &w_t::confid)
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
      find_atom_group_index_overloads, find_atom_group_index, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      atoms_interleaved_conf_overloads, atoms_interleaved_conf, 0, 1)

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
        .add_property("resseq",
          make_function(get_resseq), make_function(set_resseq))
        .add_property("icode",
          make_function(get_icode), make_function(set_icode))
        .add_property("link_to_previous",
          make_function(get_link_to_previous),
          make_function(set_link_to_previous))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue_group, chain>::wrapper)
        IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(atom_group)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
        .def("resid", &w_t::resid)
        .def("have_conformers", &w_t::have_conformers)
        .def("merge_atom_groups", &w_t::merge_atom_groups, (
          arg_("primary"), arg_("secondary")))
        .def("move_blank_altloc_atom_groups_to_front",
          &w_t::move_blank_altloc_atom_groups_to_front)
        .def("edit_blank_altloc", &w_t::edit_blank_altloc)
        .def("atoms_interleaved_conf", &w_t::atoms_interleaved_conf,
          atoms_interleaved_conf_overloads((
            arg_("group_residue_names")=true)))
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
      find_residue_group_index_overloads, find_residue_group_index, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      find_pure_altloc_ranges_overloads, find_pure_altloc_ranges, 0, 1)

    static
    boost::python::object
    conformers(
      w_t const& self)
    {
      af::shared<conformer> result = self.conformers();
      return scitbx::boost_python::array_as_list(
        result.begin(), result.size());
    }

    static void
    append_atom_record_groups(
      w_t const& self,
      boost::python::list pdb_records,
      int interleaved_conf=0,
      bool atom_hetatm=true,
      bool sigatm=true,
      bool anisou=true,
      bool siguij=true)
    {
      boost::python::ssize_t max_str_len = 0;
      if (atom_hetatm) max_str_len += 81;
      if (sigatm) max_str_len += 81;
      if (anisou) max_str_len += 81;
      if (siguij) max_str_len += 81;
      if (max_str_len == 0) max_str_len = 1;
      atom_label_columns_formatter label_formatter;
      label_formatter.chain_id = self.data->id.c_str();
      unsigned n_rg = self.residue_groups_size();
      for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
        residue_group const& rg = self.residue_groups()[i_rg];
        if (i_rg != 0 && !rg.data->link_to_previous) {
          pdb_records.append("BREAK");
        }
        label_formatter.resseq = rg.data->resseq.elems;
        label_formatter.icode = rg.data->icode.elems;
        if (interleaved_conf <= 1) {
          unsigned n_ag = rg.atom_groups_size();
          for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
            atom_group const& ag = rg.atom_groups()[i_ag];
            label_formatter.altloc = ag.data->altloc.elems;
            label_formatter.resname = ag.data->resname.elems;
            typedef std::vector<atom> va;
            va const& atoms = ag.atoms();
            va::const_iterator atoms_end = atoms.end();
            for(va::const_iterator atom=atoms.begin();atom!=atoms_end;atom++) {
#define IOTBX_LOC \
              boost::python::handle<> str_hdl(PyString_FromStringAndSize( \
                0, max_str_len)); \
              PyObject* str_obj = str_hdl.get(); \
              char* str_begin = PyString_AS_STRING(str_obj); \
              unsigned str_len = atom->format_atom_record_group( \
                str_begin, &label_formatter, \
                atom_hetatm, sigatm, anisou, siguij); \
              str_hdl.release(); \
              if (_PyString_Resize(&str_obj, static_cast<int>(str_len))!=0) { \
                boost::python::throw_error_already_set(); \
              } \
              pdb_records.append(boost::python::handle<>(str_obj));
              IOTBX_LOC
            }
          }
        }
        else { // interleaved_conf > 0
          af::shared<atom> atoms_ilc = rg.atoms_interleaved_conf(
            /* group_residue_names */ (interleaved_conf < 2));
          af::const_ref<atom> ats = atoms_ilc.const_ref();
          unsigned n_at = static_cast<unsigned>(ats.size());
          for(unsigned i_at=0;i_at<n_at;i_at++) {
            hierarchy_v2::atom const* atom = &ats[i_at];
            shared_ptr<atom_group_data> ag_data = atom->parent_ptr();
            label_formatter.altloc = ag_data->altloc.elems;
            label_formatter.resname = ag_data->resname.elems;
            IOTBX_LOC
          }
        }
      }
#undef IOTBX_LOC
    }

    BOOST_PYTHON_FUNCTION_OVERLOADS(
      append_atom_record_groups_overloads, append_atom_record_groups, 2, 7)

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
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<chain, model>::wrapper)
        IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(residue_group)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
        .def("merge_residue_groups", &w_t::merge_residue_groups, (
          arg_("primary"), arg_("secondary")))
        .def("merge_disconnected_residue_groups_with_pure_altloc",
          &w_t::merge_disconnected_residue_groups_with_pure_altloc)
        .def("find_pure_altloc_ranges", &w_t::find_pure_altloc_ranges,
          find_pure_altloc_ranges_overloads((
            arg_("common_residue_name_class_only")=0)))
        .def("conformers", conformers)
        .def("append_atom_record_groups", append_atom_record_groups,
          append_atom_record_groups_overloads((
          arg_("self"),
          arg_("pdb_records"),
          arg_("interleaved_conf")=0,
          arg_("atom_hetatm")=true,
          arg_("sigatm")=true,
          arg_("anisou")=true,
          arg_("siguij")=true)))
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

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      find_chain_index_overloads, find_chain_index, 1, 2)

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
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<model, root>::wrapper)
        IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(chain)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
        .def("is_identical_topology", &w_t::is_identical_topology, (
          arg_("other")))
        .def("transfer_chains_from_other", &w_t::transfer_chains_from_other, (
          arg_("other")))
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

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      find_model_index_overloads, find_model_index, 1, 2)

    static void
    get_overall_counts(
      w_t const& self,
      boost::python::object result)
    {
      using scitbx::boost_python::array_as_list;
      using scitbx::boost_python::stl_map_as_dict;
#define IOTBX_LOC_SA(N) \
      result.attr(#N) = oc.N;
      //
#define IOTBX_LOC_SAA(N) \
      result.attr(#N) = array_as_list(oc.N.begin(), oc.N.size());
      //
#define IOTBX_LOC_SAM(N) \
      result.attr(#N) = stl_map_as_dict(oc.N);
      //
#define IOTBX_LOC_SAO(N) \
      { \
        boost::python::object v; \
        if (oc.N) v = boost::python::object(*oc.N); \
        result.attr(#N) = v; \
      }
      //
      hierarchy_v2::overall_counts oc(self);
      IOTBX_LOC_SA(root)
      IOTBX_LOC_SA(n_empty_models)
      IOTBX_LOC_SA(n_empty_chains)
      IOTBX_LOC_SA(n_empty_residue_groups)
      IOTBX_LOC_SA(n_empty_atom_groups)
      IOTBX_LOC_SA(n_duplicate_model_ids)
      IOTBX_LOC_SA(n_duplicate_chain_ids)
      IOTBX_LOC_SA(n_duplicate_atom_labels)
      IOTBX_LOC_SAA(duplicate_atom_labels)
      IOTBX_LOC_SA(n_models)
      IOTBX_LOC_SA(n_chains)
      IOTBX_LOC_SA(n_alt_conf)
      IOTBX_LOC_SA(n_residues)
      IOTBX_LOC_SA(n_residue_groups)
      IOTBX_LOC_SA(n_explicit_chain_breaks)
      IOTBX_LOC_SA(n_atoms)
      IOTBX_LOC_SAM(model_ids)
      IOTBX_LOC_SAM(chain_ids)
      IOTBX_LOC_SAM(alt_conf_ids)
      IOTBX_LOC_SAM(resnames)
      IOTBX_LOC_SAM(resname_classes)
      IOTBX_LOC_SAM(element_charge_types)
      IOTBX_LOC_SA(n_alt_conf_none)
      IOTBX_LOC_SA(n_alt_conf_pure)
      IOTBX_LOC_SA(n_alt_conf_proper)
      IOTBX_LOC_SA(n_alt_conf_improper)
      IOTBX_LOC_SAO(alt_conf_proper)
      IOTBX_LOC_SAO(alt_conf_improper)
      {
        boost::python::list l;
        std::size_t n = oc.consecutive_residue_groups_with_same_resid.size();
        for(std::size_t i=0;i<n;i++) {
          af::tiny<residue_group, 2> const&
            rgs = oc.consecutive_residue_groups_with_same_resid[i];
          l.append(boost::python::make_tuple(rgs[0], rgs[1]));
        }
        result.attr("consecutive_residue_groups_with_same_resid") = l;
      }
      IOTBX_LOC_SA(n_chains_with_mix_of_proper_and_improper_alt_conf)
      IOTBX_LOC_SAA(residue_groups_with_multiple_resnames_using_same_altloc)
      //
#undef IOTBX_LOC_SA
#undef IOTBX_LOC_SAA
#undef IOTBX_LOC_SAM
#undef IOTBX_LOC_SAO
    }

    static void
    get_atom_selection_cache(
      w_t const& self,
      boost::python::object result)
    {
      atom_selection_cache asc(self);
#define IOTBX_LOC(A) \
      result.attr(#A) = asc.A;
      IOTBX_LOC(n_seq)
      IOTBX_LOC(name)
      IOTBX_LOC(altloc)
      IOTBX_LOC(resname)
      IOTBX_LOC(chain_id)
      IOTBX_LOC(resseq)
      IOTBX_LOC(icode)
      IOTBX_LOC(resid)
      IOTBX_LOC(segid)
      IOTBX_LOC(model_id)
      IOTBX_LOC(element)
      IOTBX_LOC(charge)
      IOTBX_LOC(anisou)
#undef IOTBX_LOC
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("root", no_init)
        .def(init<>())
        .add_property("info", make_function(get_info), make_function(set_info))
        .def("deep_copy", &w_t::deep_copy)
        .def("memory_id", &w_t::memory_id)
        IOTBX_PDB_HIERARCHY_V2_DEF_APPEND_ETC(model)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
        .def("get_overall_counts", get_overall_counts)
        .def("get_atom_selection_cache", get_atom_selection_cache)
      ;
    }
  };

  struct residue_wrappers
  {
    typedef residue w_t;

    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET(resname)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET(resseq)
    IOTBX_PDB_HIERARCHY_V2_DATA_WRAPPERS_SMALL_STR_GET(icode)

    static bool
    get_link_to_previous(w_t const& self)
    {
      return self.data->link_to_previous;
    }

    static bool
    get_is_pure_main_conf(w_t const& self)
    {
      return self.data->is_pure_main_conf;
    }

    static
    af::shared<atom>
    get_atoms(w_t const& self)
    {
      return std_vector_as_af_shared(self.atoms());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("residue", no_init)
        .add_property("resname", make_function(get_resname))
        .add_property("resseq", make_function(get_resseq))
        .add_property("icode", make_function(get_icode))
        .add_property("link_to_previous", make_function(get_link_to_previous))
        .add_property("is_pure_main_conf",make_function(get_is_pure_main_conf))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue, conformer>::wrapper)
        .def("atoms", get_atoms)
        .def("atoms_size", &w_t::atoms_size)
        .def("resid", &w_t::resid)
      ;
    }
  };

  struct conformer_wrappers
  {
    typedef conformer w_t;

    static std::string
    get_altloc(w_t const& self) { return self.data->altloc; }

    IOTBX_PDB_HIERARCHY_V2_GET_CHILDREN(conformer, residue, residues)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("conformer", no_init)
        .def(init<chain const&, std::string const&>((
          arg_("parent"), arg_("altloc"))))
        .add_property("altloc", make_function(get_altloc))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<conformer, chain>::wrapper)
        .def("residues_size", &w_t::residues_size)
        .def("residues", get_residues)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
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

    residue_wrappers::wrap();
    conformer_wrappers::wrap();

    atoms::bpl_wrap();
  }

}}}} // namespace iotbx::pdb::hierarchy_v2::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_hierarchy_v2_ext)
{
  iotbx::pdb::hierarchy_v2::wrap_hierarchy_v2();
}
