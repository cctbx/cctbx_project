#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/return_arg.hpp>
#include <iotbx/pdb/hierarchy.h>
#include <iotbx/pdb/hierarchy_bpl.h>

namespace iotbx { namespace pdb { namespace hierarchy {

namespace {

  struct atom_wrappers
  {
    typedef atom w_t;

    static boost::python::dict
    data_offsets()
    {
      boost::python::dict result;
      atom a;
      atom_data* d = a.data.get();
      char* p = reinterpret_cast<char*>(d);
#define IOTBX_LOC(attr) \
      result[static_cast<long>(reinterpret_cast<char*>(&(d->attr)) - p)] = \
        #attr;
      IOTBX_LOC(xyz)
      IOTBX_LOC(sigxyz)
      IOTBX_LOC(occ)
      IOTBX_LOC(sigocc)
      IOTBX_LOC(b)
      IOTBX_LOC(sigb)
      IOTBX_LOC(uij)
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
      IOTBX_LOC(siguij)
#endif
      IOTBX_LOC(i_seq)
      IOTBX_LOC(tmp)
      IOTBX_LOC(have_sentinel)
      IOTBX_LOC(hetero)
      IOTBX_LOC(serial)
      IOTBX_LOC(name)
      IOTBX_LOC(segid)
      IOTBX_LOC(element)
      IOTBX_LOC(charge)
#undef IOTBX_LOC
      result[atom::sizeof_data()] = "atom::sizeof_data()";
      return result;
    }

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

#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
    static sym_mat3
    get_siguij(w_t const& self) { return self.data->siguij; }

    static void
    set_siguij(w_t const& self, sym_mat3 const& new_siguij)
    {
      self.data->siguij = new_siguij;
    }
#endif

    static unsigned
    get_i_seq(w_t const& self) { return self.data->i_seq; }

    static int
    get_tmp(w_t const& self) { return self.data->tmp; }

    static void
    set_tmp(w_t const& self, int new_tmp) { self.data->tmp = new_tmp; }

    static bool
    get_hetero(w_t const& self) { return self.data->hetero; }

    static void
    set_hetero(w_t const& self, bool new_hetero)
    {
      self.data->hetero = new_hetero;
    }

    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(serial)
    IOTBX_PDB_HIERARCHY_WRAPPERS_SET_HY36(serial, data->serial, 5U,
      /* HY36_WIDTH_5_MIN */ -9999,
      /* HY36_WIDTH_5_MAX */ 87440031)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(name)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(segid)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(element)
    IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(charge)

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

    static
    boost::python::object
    format_atom_record_group(
      w_t const& self,
      bool atom_hetatm=true,
      bool sigatm=true,
      bool anisou=true,
      bool siguij=true)
    {
      boost::python::handle<> str_hdl(PyString_FromStringAndSize(0, 324));
      PyObject* str_obj = str_hdl.get();
      char* str_begin = PyString_AS_STRING(str_obj);
      unsigned str_len = self.format_atom_record_group(
        str_begin, 0, atom_hetatm, sigatm, anisou, siguij);
      str_hdl.release();
      if (_PyString_Resize(&str_obj, static_cast<int>(str_len)) != 0) {
        boost::python::throw_error_already_set();
      }
      return boost::python::object(boost::python::handle<>(str_obj));
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(id_str_overloads, id_str, 0, 2)

    BOOST_PYTHON_FUNCTION_OVERLOADS(
      format_atom_record_overloads, format_atom_record, 1, 2)

    BOOST_PYTHON_FUNCTION_OVERLOADS(
      format_atom_record_group_overloads, format_atom_record_group, 1, 5)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(quote_overloads, quote, 0, 1)

    // not inline to work around bug in
    // g++ (GCC) 3.2.3 20030502 (Red Hat Linux 3.2.3-34) x86_64
    static void
    wrap();
  };

    void
    atom_wrappers::wrap()
    {
      using namespace boost::python;
      class_<w_t>("atom", no_init)
        .def(init<>())
        .def(init<atom_group const&, atom const&>((
          arg_("parent"), arg_("other"))))
        .def("detached_copy", &w_t::detached_copy)
        .def("set_xyz", set_xyz, (arg_("new_xyz")), return_self<>())
        .def("set_sigxyz", set_sigxyz, (arg_("new_sigxyz")), return_self<>())
        .def("set_occ", set_occ, (arg_("new_occ")), return_self<>())
        .def("set_sigocc", set_sigocc, (arg_("new_sigocc")), return_self<>())
        .def("set_b", set_b, (arg_("new_b")), return_self<>())
        .def("set_sigb", set_sigb, (arg_("new_sigb")), return_self<>())
        .def("set_uij", set_uij, (arg_("new_uij")), return_self<>())
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
        .def("set_siguij", set_siguij, (arg_("new_siguij")), return_self<>())
#endif
        .def("set_hetero", set_hetero, (arg_("new_hetero")), return_self<>())
        .def("set_serial", set_serial, (arg_("new_serial")), return_self<>())
        .def("set_name", set_name, (arg_("new_name")), return_self<>())
        .def("set_segid", set_segid, (arg_("new_segid")), return_self<>())
        .def("set_element", set_element, (arg_("new_element")), return_self<>())
        .def("set_charge", set_charge, (arg_("new_charge")), return_self<>())
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
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
        .add_property("siguij",
          make_function(get_siguij), make_function(set_siguij))
#endif
        .add_property("i_seq", make_function(get_i_seq))
        .add_property("tmp",
          make_function(get_tmp), make_function(set_tmp))
        .add_property("hetero",
          make_function(get_hetero), make_function(set_hetero))
        .add_property("serial",
          make_function(get_serial), make_function(set_serial))
        .def("serial_as_int", &w_t::serial_as_int)
        .add_property("name",
          make_function(get_name), make_function(set_name))
        .add_property("segid",
          make_function(get_segid), make_function(set_segid))
        .add_property("element",
          make_function(get_element), make_function(set_element))
        .add_property("charge",
          make_function(get_charge), make_function(set_charge))
        .def("memory_id", &w_t::memory_id)
        .def("sizeof_data", &w_t::sizeof_data)
        .staticmethod("sizeof_data")
        .def("data_offsets", data_offsets)
        .staticmethod("data_offsets")
        .def("parent", get_parent<atom, atom_group>::wrapper)
        .def("uij_is_defined", &w_t::uij_is_defined)
        .def("uij_erase", &w_t::uij_erase)
        .def("has_siguij", &w_t::has_siguij)
        .staticmethod("has_siguij")
        .def("siguij_is_defined", &w_t::siguij_is_defined)
        .def("siguij_erase", &w_t::siguij_erase)
        .def("pdb_label_columns", &w_t::pdb_label_columns)
        .def("pdb_element_charge_columns", &w_t::pdb_element_charge_columns)
        .def("id_str", &w_t::id_str, id_str_overloads((
          arg_("pdbres")=false,
          arg_("suppress_segid")=false)))
        .def("format_atom_record", format_atom_record,
          format_atom_record_overloads((
            arg_("self"),
            arg_("replace_floats_with")=0)))
        .def("format_sigatm_record", format_sigatm_record)
        .def("format_anisou_record", format_anisou_record)
        .def("format_siguij_record", format_siguij_record)
        .def("format_atom_record_group", format_atom_record_group,
          format_atom_record_group_overloads((
            arg_("self"),
            arg_("atom_hetatm")=true,
            arg_("sigatm")=true,
            arg_("anisou")=true,
            arg_("siguij")=true)))
        .def("quote", &w_t::quote, quote_overloads((arg_("full")=false)))
        .def("fetch_labels", &w_t::fetch_labels)
        .def("element_is_hydrogen", &w_t::element_is_hydrogen)
        .def("determine_chemical_element_simple",
          &w_t::determine_chemical_element_simple)
      ;
    }

  struct atom_with_labels_wrappers
  {
    typedef atom_with_labels w_t;

#define IOTBX_LOC_GET(attr) \
    static \
    boost::python::str \
    get_##attr(w_t const& self) \
    { \
      return boost::python::str(self.attr.elems); \
    }

#define IOTBX_LOC_SET(attr) \
    static \
    void \
    set_##attr(w_t& self, const char* value) \
    { \
      self.attr.replace_with(value); \
    }

#define IOTBX_LOC_GET_SET(attr) \
  IOTBX_LOC_GET(attr) \
  IOTBX_LOC_SET(attr)

    IOTBX_LOC_GET_SET(model_id)
    IOTBX_LOC_GET_SET(chain_id)
    IOTBX_PDB_HIERARCHY_WRAPPERS_SET_HY36(resseq, resseq, 4U,
      /* HY36_WIDTH_4_MIN */ -999,
      /* HY36_WIDTH_4_MAX */ 2436111)
    IOTBX_LOC_GET(resseq)
    IOTBX_LOC_GET_SET(icode)
    IOTBX_LOC_GET_SET(altloc)
    IOTBX_LOC_GET_SET(resname)

#undef IOTBX_LOC_GET
#undef IOTBX_LOC_SET
#undef IOTBX_LOC_GET_SET

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(id_str_overloads, id_str, 0, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      format_atom_record_overloads, format_atom_record, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      format_atom_record_group_overloads, format_atom_record_group, 0, 4)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(quote_overloads, quote, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<atom> >("atom_with_labels", no_init)
        .def(init<>())
        .def(init<
          atom const&,
          const char*, const char*, const char*,
          const char*, const char*, const char*,
          bool, bool>((
            arg_("atom"),
            arg_("model_id"),
            arg_("chain_id"),
            arg_("resseq"),
            arg_("icode"),
            arg_("altloc"),
            arg_("resname"),
            arg_("is_first_in_chain"),
            arg_("is_first_after_break"))))
        .def("detached_copy", &w_t::detached_copy)
        .add_property("model_id",
          make_function(get_model_id), make_function(set_model_id))
        .add_property("chain_id",
          make_function(get_chain_id), make_function(set_chain_id))
        .add_property("resseq",
          make_function(get_resseq), make_function(set_resseq))
        .add_property("icode",
          make_function(get_icode), make_function(set_icode))
        .add_property("altloc",
          make_function(get_altloc), make_function(set_altloc))
        .add_property("resname",
          make_function(get_resname), make_function(set_resname))
        .def_readwrite("is_first_in_chain", &w_t::is_first_in_chain)
        .def_readwrite("is_first_after_break", &w_t::is_first_after_break)
        .def("resseq_as_int", &w_t::resseq_as_int)
        .def("resid", &w_t::resid)
        .def("id_str", &w_t::id_str, id_str_overloads((
          arg_("pdbres")=false,
          arg_("suppress_segid")=false)))
        .def("format_atom_record", &w_t::format_atom_record,
          format_atom_record_overloads((
            arg_("replace_floats_with")=0)))
        .def("format_sigatm_record", &w_t::format_sigatm_record)
        .def("format_anisou_record", &w_t::format_anisou_record)
        .def("format_siguij_record", &w_t::format_siguij_record)
        .def("format_atom_record_group", &w_t::format_atom_record_group,
          format_atom_record_group_overloads((
            arg_("atom_hetatm")=true,
            arg_("sigatm")=true,
            arg_("anisou")=true,
            arg_("siguij")=true)))
        .def("quote", &w_t::quote, quote_overloads((arg_("full")=false)))
      ;
    }
  };

} // namespace <anonymous>

  void
  atom_bpl_wrap()
  {
    atom_wrappers::wrap();
    atom_with_labels_wrappers::wrap();
  }

}}} // namespace iotbx::pdb::hierarchy
