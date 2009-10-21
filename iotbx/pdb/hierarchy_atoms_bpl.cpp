#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/str.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/boost_python/selections_wrapper.h>
#include <iotbx/pdb/hierarchy_atoms.h>
#include <boost/format.hpp>

namespace iotbx { namespace pdb { namespace hierarchy { namespace atoms {

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    extract_element_overloads, extract_element, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(reset_serial_overloads, reset_serial, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    set_chemical_element_simple_if_necessary_overloads,
    set_chemical_element_simple_if_necessary, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(reset_tmp_overloads, reset_tmp, 1, 3)

  boost::python::dict
  build_dict(
    af::const_ref<atom> const& atoms,
    bool strip_names,
    bool upper_names,
    bool throw_runtime_error_if_duplicate_keys)
  {
    namespace bp = boost::python;
    bp::dict result;
    bp::object none;
    for(unsigned i=0;i<atoms.size();i++) {
      str4 name = atoms[i].data->name;
      if (strip_names) name = name.strip();
      if (upper_names) name.upper_in_place();
      bp::str key = name.elems;
      bp::object prev_atom = result.get(key);
      if (prev_atom.ptr() == none.ptr()) {
        result[key] = atoms[i];
      }
      else if (throw_runtime_error_if_duplicate_keys) {
        throw std::runtime_error((boost::format(
          "Duplicate keys in build_dict(strip_names=%s, upper_names=%s):\n"
          "  %s\n"
          "  %s")
            % (strip_names ? "true" : "false")
            % (upper_names ? "true" : "false")
            % bp::extract<atom const&>(prev_atom)().id_str()
            % atoms[i].id_str()).str());
      }
    }
    return result;
  }

} // namespace <anonymous>

  void
  bpl_wrap()
  {
    using namespace boost::python;
    class_<atom_tmp_sentinel,
           std::auto_ptr<atom_tmp_sentinel>,
           boost::noncopyable>("atom_data_tmp_sentinel", no_init);
    typedef scitbx::af::boost_python::shared_wrapper<atom> wat;
    class_<wat::w_t> wa = wat::wrap("af_shared_atom");
    scitbx::af::boost_python::select_wrappers<
      atom, af::shared<atom> >::wrap(wa);
    wa.def("extract_serial", extract_serial)
      .def("extract_name", extract_name)
      .def("extract_xyz", extract_xyz)
      .def("extract_sigxyz", extract_sigxyz)
      .def("extract_occ", extract_occ)
      .def("extract_sigocc", extract_sigocc)
      .def("extract_b", extract_b)
      .def("extract_sigb", extract_sigb)
      .def("extract_uij", extract_uij)
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
      .def("extract_siguij", extract_siguij)
#endif
      .def("extract_hetero", extract_hetero)
      .def("extract_element", extract_element, extract_element_overloads((
        arg_("self"),
        arg_("strip")=false)))
      .def("extract_i_seq", extract_i_seq)
      .def("extract_tmp_as_size_t", extract_tmp_as_size_t)
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
      .def("reset_serial", reset_serial, reset_serial_overloads((
        arg_("self"),
        arg_("first_value")=1)))
      .def("set_chemical_element_simple_if_necessary",
        set_chemical_element_simple_if_necessary,
          set_chemical_element_simple_if_necessary_overloads((
            arg_("self"),
            arg_("tidy_existing")=true)))
      .def("reset_i_seq", reset_i_seq)
      .def("reset_tmp", reset_tmp, reset_tmp_overloads((
        arg_("self"),
        arg_("first_value")=0,
        arg_("increment")=1)))
      .def("reset_tmp_for_occupancy_groups_simple",
        reset_tmp_for_occupancy_groups_simple)
      .def("build_dict", build_dict, (
        arg("strip_names")=false,
        arg("upper_names")=false,
        arg("throw_runtime_error_if_duplicate_keys")=true));
    ;
  }

}}}} // namespace iotbx::pdb::hierarchy::atoms
