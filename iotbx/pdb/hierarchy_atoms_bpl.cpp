#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/adptbx.h>
#include <cctbx/uctbx.h>

#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/str.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/boost_python/selections_wrapper.h>
#include <iotbx/pdb/hierarchy_atoms.h>
#include <boost/format.hpp>

namespace iotbx { namespace pdb { namespace hierarchy { namespace atoms {

namespace {

  boost::python::dict
  build_dict(
    af::const_ref<atom> const& atoms,
    bool strip_names,
    bool upper_names,
    bool convert_stars_to_primes,
    bool throw_runtime_error_if_duplicate_keys)
  {
    namespace bp = boost::python;
    bp::dict result;
    bp::object none;
    for(unsigned i=0;i<atoms.size();i++) {
      str4 name = atoms[i].data->name;
      if (strip_names) name = name.strip();
      if (upper_names) name.upper_in_place();
      if (convert_stars_to_primes) name.replace_in_place('*', '\'');
      bp::str key = name.elems;
      bp::object prev_atom = result.get(key);
      if (prev_atom.ptr() == none.ptr()) {
        result[key] = atoms[i];
      }
      else if (throw_runtime_error_if_duplicate_keys) {
        throw std::runtime_error((boost::format(
          "Duplicate keys in build_dict("
          "strip_names=%s, upper_names=%s, convert_stars_to_primes=%s):\n"
          "  %s\n"
          "  %s")
            % (strip_names ? "true" : "false")
            % (upper_names ? "true" : "false")
            % (convert_stars_to_primes ? "true" : "false")
            % bp::extract<atom const&>(prev_atom)().id_str()
            % atoms[i].id_str()).str());
      }
    }
    return result;
  }

  void
  set_adps_from_scatterers(
    af::const_ref<atom> const& atoms,
    af::const_ref<cctbx::xray::scatterer<> > const& scatterers,
    cctbx::uctbx::unit_cell const& unit_cell)
  {
    namespace adptbx = cctbx::adptbx;
    for (unsigned i = 0; i < atoms.size(); i++) {
      if (scatterers[i].flags.use_u_iso()) {
        atoms[i].data->b = adptbx::u_as_b(scatterers[i].u_iso);
        atoms[i].uij_erase();
      } else if (scatterers[i].flags.use_u_aniso()) {
        atoms[i].data->uij = adptbx::u_star_as_u_cart(unit_cell,
          scatterers[i].u_star);
        atoms[i].data->b = adptbx::u_as_b(adptbx::u_cart_as_u_iso(
          atoms[i].data->uij));
      }
    }
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
      .def("extract_element", extract_element, (arg("strip")=false))
      .def("extract_i_seq", extract_i_seq)
      .def("extract_tmp_as_size_t", extract_tmp_as_size_t)
      .def("set_xyz", set_xyz, (arg("new_xyz")), return_self<>())
      .def("set_sigxyz", set_sigxyz, (arg("new_sigxyz")), return_self<>())
      .def("set_occ", set_occ, (arg("new_occ")), return_self<>())
      .def("set_sigocc", set_sigocc, (arg("new_sigocc")), return_self<>())
      .def("set_b", set_b, (arg("new_b")), return_self<>())
      .def("set_sigb", set_sigb, (arg("new_sigb")), return_self<>())
      .def("set_uij", set_uij, (arg("new_uij")), return_self<>())
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
      .def("set_siguij", set_siguij, (arg("new_siguij")), return_self<>())
#endif
      .def("reset_serial", reset_serial, (arg("first_value")=1))
      .def("set_chemical_element_simple_if_necessary",
        set_chemical_element_simple_if_necessary, (
          arg("tidy_existing")=true))
      .def("reset_i_seq", reset_i_seq)
      .def("reset_tmp", reset_tmp, (
        arg("first_value")=0,
        arg("increment")=1))
      .def("reset_tmp_for_occupancy_groups_simple",
        reset_tmp_for_occupancy_groups_simple)
      .def("build_dict", build_dict, (
        arg("strip_names")=false,
        arg("upper_names")=false,
        arg("convert_stars_to_primes")=false,
        arg("throw_runtime_error_if_duplicate_keys")=true))
      .def("set_adps_from_scatterers", set_adps_from_scatterers, (
        arg("scatterers"),
        arg("unit_cell")));
    ;
  }

}}}} // namespace iotbx::pdb::hierarchy::atoms
