#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_arg.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/crystal/direct_space_asu.h>
#include <cctbx/crystal/workarounds_bpl.h>

namespace cctbx { namespace crystal { namespace direct_space_asu {

namespace {

  struct float_cut_plane_wrappers
  {
    typedef float_cut_plane<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("direct_space_asu_float_cut_plane", no_init)
        .def(init<fractional<double> const&, double>((arg("n"), arg("c"))))
        .add_property("n",
          make_getter(&w_t::n, rbv()),
          make_setter(&w_t::n, dcp()))
        .def_readwrite("c", &w_t::c)
        .def("evaluate", &w_t::evaluate, (arg("point")))
        .def("is_inside", &w_t::is_inside, (
          arg("point"), arg("epsilon")=0))
        .def("get_point_in_plane", &w_t::get_point_in_plane)
        .def("add_buffer", &w_t::add_buffer,
          (arg("unit_cell"), arg("thickness")))
      ;
    }
  };

  struct float_asu_wrappers
  {
    typedef float_asu<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("direct_space_asu_float_asu", no_init)
        .def(init<
          uctbx::unit_cell const&,
          w_t::cuts_t const&,
          double const&>((
            arg("unit_cell"),
            arg("cuts"),
            arg("is_inside_epsilon")=1e-6)))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("cuts", &w_t::cuts, ccr())
        .def("is_inside_epsilon", &w_t::is_inside_epsilon)
        .def("is_inside", &w_t::is_inside, (arg("point")))
        .def("is_inside_frac", &w_t::is_inside_frac, (arg("sites_frac")))
        .def("is_inside_cart", &w_t::is_inside_cart, (arg("sites_cart")))
        .def("_add_buffer", &w_t::add_buffer)
        .def("volume_vertices", &w_t::volume_vertices, (
          arg("cartesian")=false,
          arg("epsilon")=1e-6))
        .def("box_min", &w_t::box_min, (arg("cartesian")=false), ccr())
        .def("box_max", &w_t::box_max, (arg("cartesian")=false), ccr())
      ;
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping<w_t::cuts_t, fixed_capacity_policy>();
      }
    }
  };

  struct asu_mapping_wrappers
  {
    typedef asu_mapping<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("direct_space_asu_asu_mapping", no_init)
        .def("i_sym_op", &w_t::i_sym_op)
        .def("unit_shifts", &w_t::unit_shifts, ccr())
        .def("mapped_site", &w_t::mapped_site, ccr())
      ;
    }
  };

  struct asu_mapping_index_pair_wrappers
  {
    typedef asu_mapping_index_pair w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("direct_space_asu_asu_mapping_index_pair", no_init)
        .def_readonly("i_seq", &w_t::i_seq)
        .def_readonly("j_seq", &w_t::j_seq)
        .def_readonly("j_sym", &w_t::j_sym)
        .def("is_active", &w_t::is_active, (arg("minimal")=false))
      ;
    }
  };

  struct asu_mapping_index_pair_and_diff_wrappers
  {
    typedef asu_mapping_index_pair_and_diff<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<asu_mapping_index_pair> >(
        "direct_space_asu_asu_mapping_index_pair_and_diff", no_init)
        .add_property("diff_vec", make_getter(&w_t::diff_vec, rbv()))
        .def_readonly("dist_sq", &w_t::dist_sq)
      ;
    }
  };

  struct asu_mappings_wrappers
  {
    typedef asu_mappings<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t, boost::shared_ptr<w_t> >(
          "direct_space_asu_asu_mappings", no_init)
        .def(init<sgtbx::space_group const&,
                  float_asu<> const&,
                  double const&>((
          arg("space_group"),
          arg("asu"),
          arg("buffer_thickness"))))
        .def("reserve", &w_t::reserve, (arg("n_sites_final")))
        .def("space_group", &w_t::space_group, rir())
        .def("asu", &w_t::asu, rir())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("buffer_thickness", &w_t::buffer_thickness)
        .def("asu_buffer", &w_t::asu_buffer, rir())
        .def("buffer_covering_sphere", &w_t::buffer_covering_sphere, rir())
        .def("process",
          (w_t&(w_t::*)(fractional<> const&, double const&)) &w_t::process, (
            arg("original_site"), arg("min_distance_sym_equiv")=0.5),
              return_self<>())
        .def("process",
          (w_t&(w_t::*)(fractional<> const&, sgtbx::site_symmetry_ops const&))
              &w_t::process, (
          arg("original_site"), arg("site_symmetry_ops")),
          return_self<>())
        .def("process_sites_frac",
          (w_t&(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            double const&))
              &w_t::process_sites_frac, (
            arg("original_sites"),
            arg("min_distance_sym_equiv")=0.5),
              return_self<>())
        .def("process_sites_frac",
          (w_t&(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_frac, (
          arg("original_sites"), arg("site_symmetry_table")),
          return_self<>())
        .def("process_sites_cart",
          (w_t&(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            double const&))
              &w_t::process_sites_cart, (
            arg("original_sites"),
            arg("min_distance_sym_equiv")=0.5),
              return_self<>())
        .def("process_sites_cart",
          (w_t&(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_cart, (
          arg("original_sites"), arg("site_symmetry_table")),
          return_self<>())
        .def("n_sites_in_asu_and_buffer", &w_t::n_sites_in_asu_and_buffer)
        .def("n_sites_in_asu_and_buffer", &w_t::n_sites_in_asu_and_buffer)
        .def("mappings", &w_t::mappings, ccr())
        .def("mapped_sites_min", &w_t::mapped_sites_min, ccr())
        .def("mapped_sites_max", &w_t::mapped_sites_max, ccr())
        .def("mapped_sites_span", &w_t::mapped_sites_span)
        .def("special_op", &w_t::special_op, (arg("i_seq")), ccr())
        .def("site_symmetry_table", &w_t::site_symmetry_table, rir())
        .def("get_rt_mx",
          (sgtbx::rt_mx(w_t::*)(asu_mapping<> const&) const)
            &w_t::get_rt_mx, (
          arg("mapping")))
        .def("get_rt_mx",
          (sgtbx::rt_mx(w_t::*)(std::size_t, std::size_t) const)
            &w_t::get_rt_mx, (
          arg("i_seq"), arg("i_sym")))
        .def("get_rt_mx_i", &w_t::get_rt_mx_i, (arg("pair")))
        .def("get_rt_mx_j", &w_t::get_rt_mx_j, (arg("pair")))
        .def("diff_vec", &w_t::diff_vec, (arg("pair")))
        .def("map_moved_site_to_asu", &w_t::map_moved_site_to_asu, (
          arg("moved_original_site"), arg("i_seq"), arg("i_sym")))
        .def("r_inv_cart", &w_t::r_inv_cart, (arg("i_seq"), arg("i_sym")))
        .def("is_simple_interaction", &w_t::is_simple_interaction, (
          arg("pair")))
        .def("make_trial_pair", &w_t::make_trial_pair, (
          arg("i_seq"), arg("j_seq"), arg("j_sym")))
        .def("make_pair", &w_t::make_pair, (
          arg("i_seq"), arg("j_seq"), arg("j_sym")))
        .def("find_i_sym", &w_t::find_i_sym, (arg("i_seq"), arg("rt_mx")))
      ;
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping_variable_capacity<
          w_t::array_of_mappings_for_one_site>();
      }
      {
        scitbx::af::boost_python::shared_wrapper<
          w_t::array_of_mappings_for_one_site>::wrap(
            "direct_space_asu_array_of_array_of_mappings_for_one_site");
      }
    }
  };

}} // namespace direct_space_asu::<anoymous>

namespace boost_python {

  void wrap_direct_space_asu()
  {
    direct_space_asu::float_cut_plane_wrappers::wrap();
    direct_space_asu::float_asu_wrappers::wrap();
    direct_space_asu::asu_mapping_wrappers::wrap();
    direct_space_asu::asu_mapping_index_pair_wrappers::wrap();
    direct_space_asu::asu_mapping_index_pair_and_diff_wrappers::wrap();
    direct_space_asu::asu_mappings_wrappers::wrap();
  }

}}} // namespace cctbx::crystal::boost_python
