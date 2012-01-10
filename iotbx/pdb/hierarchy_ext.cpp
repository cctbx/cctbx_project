#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <scitbx/boost_python/stl_map_as_dict.h>
#include <scitbx/boost_python/array_as_list.h>
#include <iotbx/pdb/hierarchy_atoms.h>
#include <iotbx/pdb/hierarchy_bpl.h>
#include <iotbx/pdb/write_utils_bpl.h>

namespace iotbx { namespace pdb { namespace hierarchy {

void atom_bpl_wrap();
namespace atoms { void bpl_wrap(); }

namespace {

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

#define IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(C) \
        .def(#C "s", get_##C##s) \
        .def(#C "s_size", &w_t::C##s_size) \
        .def("find_" #C "_index", &w_t::find_##C##_index, ( \
          arg(#C), arg("must_be_present")=false)) \
        .def("pre_allocate_" #C "s", &w_t::pre_allocate_##C##s, ( \
          arg("number_of_additional_" #C "s"))) \
        .def("insert_" #C, &w_t::insert_##C, (arg("i"), arg(#C))) \
        .def("append_" #C, &w_t::append_##C, (arg(#C))) \
        .def("remove_" #C, \
          (void(w_t::*)(long)) &w_t::remove_##C, (arg("i"))) \
        .def("remove_" #C, \
          (void(w_t::*)(C&)) &w_t::remove_##C, (arg(#C)))

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

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(altloc)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(resname)

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
      class_<w_t>("atom_group", no_init)
        .def(init<
          residue_group const&,
            optional<const char*, const char*> >((
              arg("parent"), arg("altloc")="", arg("resname")="")))
        .def(init<
          optional<const char*, const char*> >((
            arg("altloc")="", arg("resname")="")))
        .def(init<residue_group const&, atom_group const&>((
          arg("parent"), arg("other"))))
        .add_property("altloc",
          make_function(get_altloc), make_function(set_altloc))
        .add_property("resname",
          make_function(get_resname), make_function(set_resname))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<atom_group, residue_group>::wrapper, (
          arg("optional")=true))
        IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(atom)
        .def("append_atom_with_other_parent",
          &w_t::append_atom_with_other_parent, (arg("atom")))
        .def("confid", &w_t::confid)
      ;
    }
  };

  template <typename ChainOrResidueGroup>
  struct conformers_as_list
  {
    static
    boost::python::object
    get(ChainOrResidueGroup const& self)
    {
      af::shared<conformer> result = self.conformers();
      return scitbx::boost_python::array_as_list(
        result.begin(), result.size());
    }
  };

  struct residue_group_wrappers
  {
    typedef residue_group w_t;

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(resseq)
    IOTBX_PDB_HIERARCHY_WRAPPERS_SET_HY36(resseq, data->resseq, 4U,
      /* HY36_WIDTH_4_MIN */ -999,
      /* HY36_WIDTH_4_MAX */ 2436111)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(icode)

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

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(residue_group, atom_group, atom_groups)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("residue_group", no_init)
        .def(init<chain const&, optional<const char*, const char*, bool> >((
          arg("parent"),
          arg("resseq")="", arg("icode")="", arg("link_to_previous")=true)))
        .def(init<optional<const char*, const char*, bool> >((
          arg("resseq")="", arg("icode")="", arg("link_to_previous")=true)))
        .def(init<chain const&, residue_group const&>((
          arg("parent"), arg("other"))))
        .add_property("resseq",
          make_function(get_resseq), make_function(set_resseq))
        .add_property("icode",
          make_function(get_icode), make_function(set_icode))
        .add_property("link_to_previous",
          make_function(get_link_to_previous),
          make_function(set_link_to_previous))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue_group, chain>::wrapper, (
          arg("optional")=true))
        IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(atom_group)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms, (arg("interleaved_conf")=0))
        .def("resseq_as_int", &w_t::resseq_as_int)
        .def("resid", &w_t::resid)
        .def("have_conformers", &w_t::have_conformers)
        .def("merge_atom_groups", &w_t::merge_atom_groups, (
          arg("primary"), arg("secondary")))
        .def("move_blank_altloc_atom_groups_to_front",
          &w_t::move_blank_altloc_atom_groups_to_front)
        .def("edit_blank_altloc", &w_t::edit_blank_altloc)
        .def("is_identical_hierarchy", &w_t::is_identical_hierarchy, (
          arg("other")))
        .def("is_similar_hierarchy", &w_t::is_similar_hierarchy, (
          arg("other")))
        .def("conformers", conformers_as_list<w_t>::get)
        .def("unique_resnames", &w_t::unique_resnames)
      ;
    }
  };

  struct chain_wrappers
  {
    typedef chain w_t;

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(id)

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(chain, residue_group, residue_groups)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("chain", no_init)
        .def(init<model const&, optional<const char*> >((
          arg("parent"), arg("id")="")))
        .def(init<const char*>((
          arg("id")="")))
        .def(init<model const&, chain const&>((
          arg("parent"), arg("other"))))
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<chain, model>::wrapper, (
          arg("optional")=true))
        IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(residue_group)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms, (arg("interleaved_conf")=0))
        .def("merge_residue_groups", &w_t::merge_residue_groups, (
          arg("primary"), arg("secondary")))
        .def("merge_disconnected_residue_groups_with_pure_altloc",
          &w_t::merge_disconnected_residue_groups_with_pure_altloc)
        .def("find_pure_altloc_ranges", &w_t::find_pure_altloc_ranges, (
          arg("common_residue_name_class_only")=object()))
        .def("conformers", conformers_as_list<w_t>::get)
        .def("is_identical_hierarchy", &w_t::is_identical_hierarchy, (
          arg("other")))
        .def("is_similar_hierarchy", &w_t::is_similar_hierarchy, (
          arg("other")))
      ;
    }
  };

  struct model_wrappers
  {
    typedef model w_t;

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(id)

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(model, chain, chains)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("model", no_init)
        .def(init<root const&, optional<const char*> >((
          arg("parent"), arg("id")="")))
        .def(init<const char*>((arg("id")="")))
        .def(init<root const&, model const&>((
          arg("parent"), arg("other"))))
        .add_property("id", make_function(get_id), make_function(set_id))
        .def("detached_copy", &w_t::detached_copy)
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<model, root>::wrapper, (
          arg("optional")=true))
        IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(chain)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms, (arg("interleaved_conf")=0))
        .def("is_identical_hierarchy", &w_t::is_identical_hierarchy, (
          arg("other")))
        .def("is_similar_hierarchy", &w_t::is_similar_hierarchy, (
          arg("other")))
        .def("transfer_chains_from_other", &w_t::transfer_chains_from_other, (
          arg("other")))
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

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(root, model, models)

    static void
    as_pdb_string_cstringio(
      w_t const& self,
      boost::python::object cstringio,
      bool append_end=false,
      int interleaved_conf=0,
      boost::optional<int>
        atoms_reset_serial_first_value=boost::optional<int>(),
      bool atom_hetatm=true,
      bool sigatm=true,
      bool anisou=true,
      bool siguij=true)
    {
      if (atoms_reset_serial_first_value) {
        self.atoms_reset_serial(
          interleaved_conf, *atoms_reset_serial_first_value);
      }
      write_utils::cstringio_write write(cstringio.ptr());
      models_as_pdb_string(
        write,
        self.models(),
        append_end,
        interleaved_conf,
        atom_hetatm,
        sigatm,
        anisou,
        siguij);
    }

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
      hierarchy::overall_counts oc(self);
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
      IOTBX_LOC_SA(n_anisou)
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
      IOTBX_LOC(resid_list)
      IOTBX_LOC(chain_break_list)
#undef IOTBX_LOC
    }

    static
    boost::python::object
    altloc_indices(
      w_t const& self)
    {
      return boost::python::object(
        atom_selection_cache(self, /* altloc_only */ true).altloc);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("root", no_init)
        .def(init<>())
        .enable_pickling()
        .add_property("info", make_function(get_info), make_function(set_info))
        .def("deep_copy", &w_t::deep_copy)
        .def("memory_id", &w_t::memory_id)
        IOTBX_PDB_HIERARCHY_DEF_APPEND_ETC(model)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms, (arg("interleaved_conf")=0))
        .def("atoms_with_i_seq_mismatch", &w_t::atoms_with_i_seq_mismatch)
        .def("atoms_reset_serial", &w_t::atoms_reset_serial, (
          arg("interleaved_conf")=0,
          arg("first_value")=1))
        .def("is_similar_hierarchy", &w_t::is_similar_hierarchy, (
          arg("other")))
        .def("_as_pdb_string_cstringio", as_pdb_string_cstringio, (
          arg("self"),
          arg("cstringio"),
          arg("append_end"),
          arg("interleaved_conf"),
          arg("atoms_reset_serial_first_value"),
          arg("atom_hetatm"),
          arg("sigatm"),
          arg("anisou"),
          arg("siguij")))
        .def("_write_pdb_file", &w_t::write_pdb_file, (
          arg("file_name"),
          arg("open_append"),
          arg("append_end"),
          arg("interleaved_conf"),
          arg("atoms_reset_serial_first_value"),
          arg("atom_hetatm"),
          arg("sigatm"),
          arg("anisou"),
          arg("siguij")))
        .def("select",
          (root(w_t::*)(af::const_ref<bool> const&, bool) const)
            &w_t::select, (
              arg("atom_selection"), arg("copy_atoms")=false))
        .def("select",
          (root(w_t::*)(af::const_ref<std::size_t> const&, bool) const)
            &w_t::select, (
              arg("atom_selection"), arg("copy_atoms")=false))
        .def("get_overall_counts", get_overall_counts)
        .def("get_atom_selection_cache", get_atom_selection_cache)
        .def("altloc_indices", altloc_indices)
      ;
      def("get_resid_sequence", get_resid_sequence, (
        arg("resid_list"),
        arg("chain_break_list"),
        arg("start"),
        arg("stop")));
    }
  };

  struct residue_wrappers
  {
    typedef residue w_t;

    static boost::python::object
    get_root(w_t const& self)
    {
      boost::optional<root> result = self.root();
      if (!result) return boost::python::object();
      return boost::python::object(*result);
    }

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(resname)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(resseq)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(icode)

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
        .enable_pickling()
        .def(init<root const&>((arg("root"))))
        .def("root", get_root)
        .add_property("resname", make_function(get_resname))
        .add_property("resseq", make_function(get_resseq))
        .add_property("icode", make_function(get_icode))
        .add_property("link_to_previous", make_function(get_link_to_previous))
        .add_property("is_pure_main_conf",make_function(get_is_pure_main_conf))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<residue, conformer>::wrapper, (
          arg("optional")=true))
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", get_atoms)
        .def("resseq_as_int", &w_t::resseq_as_int)
        .def("resid", &w_t::resid)
        .def("id_str", &w_t::id_str, (arg("suppress_segid")=0))
        .def("find_atom_by", &w_t::find_atom_by, (arg("name")))
      ;
    }
  };

  struct conformer_wrappers
  {
    typedef conformer w_t;

    static std::string
    get_altloc(w_t const& self) { return self.data->altloc; }

    IOTBX_PDB_HIERARCHY_GET_CHILDREN(conformer, residue, residues)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("conformer", no_init)
        .def(init<chain const&, std::string const&>((
          arg("parent"), arg("altloc"))))
        .add_property("altloc", make_function(get_altloc))
        .def("memory_id", &w_t::memory_id)
        .def("parent", get_parent<conformer, chain>::wrapper, (
          arg("optional")=true))
        .def("residues_size", &w_t::residues_size)
        .def("residues", get_residues)
        .def("atoms_size", &w_t::atoms_size)
        .def("atoms", &w_t::atoms)
      ;
    }
  };

  void
  wrap_hierarchy()
  {
    atom_bpl_wrap();
    atom_group_wrappers::wrap();
    residue_group_wrappers::wrap();
    chain_wrappers::wrap();
    model_wrappers::wrap();
    root_wrappers::wrap();

    residue_wrappers::wrap();
    conformer_wrappers::wrap();

    atoms::bpl_wrap();
  }

}}}} // namespace iotbx::pdb::hierarchy::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_hierarchy_ext)
{
  scitbx::boost_python::cstringio_import();
  iotbx::pdb::hierarchy::wrap_hierarchy();
}
