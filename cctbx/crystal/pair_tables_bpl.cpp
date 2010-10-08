#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/stl/map_wrapper.h>
#include <scitbx/stl/vector_wrapper.h>
#include <cctbx/crystal/pair_tables.h>
#include <cctbx/crystal/workarounds_bpl.h>

namespace cctbx { namespace crystal {
namespace {

  struct pair_sym_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_dict");
      scitbx::af::boost_python::shared_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_table")
        .def("proxy_select",
          (pair_sym_table(*)(
            af::const_ref<pair_sym_dict> const&,
            af::const_ref<std::size_t> const&))
              scitbx::af::array_of_map_proxy_select)
      ;
    }
  };

  struct pair_asu_table_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_dict");
      scitbx::af::boost_python::shared_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_table_table");
    }
  };

  struct pair_asu_table_wrappers
  {
    typedef pair_asu_table<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t, boost::shared_ptr<w_t> >("pair_asu_table", no_init)
        .def(init<
          boost::shared_ptr<direct_space_asu::asu_mappings<> > >(
            (arg("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings)
        .def("table", &w_t::table, ccr())
        .def("__contains__",
          (bool(w_t::*)(direct_space_asu::asu_mapping_index_pair const&) const)
            &w_t::contains, (
          arg("pair")))
        .def("contains",
          (bool(w_t::*)(unsigned, unsigned, unsigned) const)
            &w_t::contains, (
          arg("i_seq"), arg("j_seq"), arg("j_sym")))
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("pair_counts", &w_t::pair_counts)
        .def("cluster_pivot_selection", &w_t::cluster_pivot_selection, (
          arg("general_positions_only")=false,
          arg("max_clusters")=0,
          arg("estimated_reduction_factor")=4))
        .def("add_covalent_pairs", &w_t::add_covalent_pairs, (
          arg("scattering_types"),
          arg("exclude_scattering_types")=boost::python::object(),
          arg("distance_cutoff")=3.5,
          arg("min_cubicle_edge")=5,
          arg("tolerance")=0.5,
          arg("epsilon")=1e-6), return_self<>())
        .def("add_all_pairs", &w_t::add_all_pairs, (
          arg("distance_cutoff"),
          arg("min_cubicle_edge")=5,
          arg("epsilon")=1e-6), return_self<>())
        .def("add_pair_sym_table", &w_t::add_pair_sym_table, (
          arg("sym_table")), return_self<>())
        .def("add_pair", (pair_asu_table<>&(w_t::*)(
            direct_space_asu::asu_mapping_index_pair const&)) &w_t::add_pair,
          (arg("pair")), return_self<>())
        .def("add_pair", (pair_asu_table<>&(w_t::*)(
            unsigned, unsigned, sgtbx::rt_mx const&)) &w_t::add_pair,
          (arg("i_seq"), arg("j_seq"), arg("rt_mx_ji")), return_self<>())
        .def("add_pair", (pair_asu_table<>&(w_t::*)(
            af::tiny<unsigned, 2> const&)) &w_t::add_pair,
          (arg("i_seqs")), return_self<>())
        .def("extract_pair_sym_table", &w_t::extract_pair_sym_table, (
          arg("skip_j_seq_less_than_i_seq")=true,
          arg("all_interactions_from_inside_asu")=false))
        .def("angle_pair_asu_table", &w_t::angle_pair_asu_table)
      ;
    }
  };

  struct adp_iso_local_sphere_restraints_energies_wrappers
  {
    typedef adp_iso_local_sphere_restraints_energies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("adp_iso_local_sphere_restraints_energies", no_init)
        .def(init<
          af::const_ref<pair_sym_dict> const&,
          scitbx::mat3<double> const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<bool> const&,
          af::const_ref<bool> const&,
          double,
          double,
          double,
          double,
          bool,
          bool>((
            arg("pair_sym_table"),
            arg("orthogonalization_matrix"),
            arg("sites_frac"),
            arg("u_isos"),
            arg("selection"),
            arg("use_u_iso"),
            arg("grad_u_iso"),
            arg("sphere_radius"),
            arg("distance_power"),
            arg("average_power"),
            arg("min_u_sum"),
            arg("compute_gradients"),
            arg("collect"))))
        .def_readonly("number_of_restraints", &w_t::number_of_restraints)
        .def_readonly("residual_sum", &w_t::residual_sum)
        .add_property("gradients", make_getter(&w_t::gradients, rbv()))
        .add_property("u_i", make_getter(&w_t::u_i, rbv()))
        .add_property("u_j", make_getter(&w_t::u_j, rbv()))
        .add_property("r_ij", make_getter(&w_t::r_ij, rbv()))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;

    def("get_distances",
      (af::shared<double>(*)(
        af::const_ref<crystal::pair_sym_dict> const&,
        scitbx::mat3<double> const&,
        af::const_ref<scitbx::vec3<double> > const&)) get_distances, (
      arg("pair_sym_table"),
      arg("orthogonalization_matrix"),
      arg("sites_frac")));
    def("get_distances",
      (af::shared<double>(*)(
        af::const_ref<crystal::pair_sym_dict> const&,
        af::const_ref<scitbx::vec3<double> > const&)) get_distances, (
      arg("pair_sym_table"),
      arg("sites_cart")));
    pair_sym_table_wrappers::wrap();

    pair_asu_table_table_wrappers::wrap();
    pair_asu_table_wrappers::wrap();
    adp_iso_local_sphere_restraints_energies_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_pair_tables() { wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
