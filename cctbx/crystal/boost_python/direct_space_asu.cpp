#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace crystal { namespace direct_space_asu {

namespace {

  struct float_cut_plane_wrappers
  {
    typedef float_cut_plane<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_inside_overloads, is_inside, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("direct_space_asu_float_cut_plane", no_init)
        .def(init<fractional<double> const&, double>((arg_("n"), arg_("c"))))
        .add_property("n",
          make_getter(&w_t::n, rbv()),
          make_setter(&w_t::n, dcp()))
        .def_readwrite("c", &w_t::c)
        .def("evaluate", &w_t::evaluate, (arg_("point")))
        .def("is_inside", &w_t::is_inside, is_inside_overloads(
          (arg_("point"), arg_("epsilon"))))
        .def("get_point_in_plane", &w_t::get_point_in_plane)
        .def("add_buffer", &w_t::add_buffer,
          (arg_("unit_cell"), arg_("thickness")))
      ;
    }
  };

  struct float_asu_wrappers
  {
    typedef float_asu<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      volume_vertices_overloads, volume_vertices, 0, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      box_min_overloads, box_min, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      box_max_overloads, box_max, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("direct_space_asu_float_asu", no_init)
        .def(init<uctbx::unit_cell const&,
                  w_t::facets_t const&,
                  optional<double const&> >(
          (arg_("unit_cell"), arg_("facets"), arg_("is_inside_epsilon"))))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("facets", &w_t::facets, ccr())
        .def("is_inside_epsilon", &w_t::is_inside_epsilon)
        .def("is_inside", &w_t::is_inside, (arg_("point")))
        .def("_add_buffer", &w_t::add_buffer)
        .def("volume_vertices", &w_t::volume_vertices,
          volume_vertices_overloads((arg_("cartesian"), arg_("epsilon"))))
        .def("box_min", &w_t::box_min, box_min_overloads(
          (arg_("cartesian")))[ccr()])
        .def("box_max", &w_t::box_max, box_max_overloads(
          (arg_("cartesian")))[ccr()])
      ;
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping<w_t::facets_t, fixed_capacity_policy>();
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

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_active_overloads, is_active, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("direct_space_asu_asu_mapping_index_pair", no_init)
        .def_readonly("i_seq", &w_t::i_seq)
        .def_readonly("j_seq", &w_t::j_seq)
        .def_readonly("j_sym", &w_t::j_sym)
        .def("is_active", &w_t::is_active, is_active_overloads(
          (arg_("minimal"))))
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

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_overloads, process, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_frac_overloads, process_sites_frac, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_sites_cart_overloads, process_sites_cart, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t, boost::shared_ptr<asu_mappings<> > >(
          "direct_space_asu_asu_mappings", no_init)
        .def(init<sgtbx::space_group const&,
                  float_asu<> const&,
                  double const&>((
          arg_("space_group"),
          arg_("asu"),
          arg_("buffer_thickness"))))
        .def("reserve", &w_t::reserve, (arg_("n_sites_final")))
        .def("space_group", &w_t::space_group, rir())
        .def("asu", &w_t::asu, rir())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("buffer_thickness", &w_t::buffer_thickness)
        .def("asu_buffer", &w_t::asu_buffer, rir())
        .def("buffer_covering_sphere", &w_t::buffer_covering_sphere, rir())
        .def("process",
          (void(w_t::*)(fractional<> const&, double const&))
            &w_t::process, process_overloads((
          arg_("original_site"), arg_("min_distance_sym_equiv"))))
        .def("process",
          (void(w_t::*)(fractional<> const&, sgtbx::site_symmetry_ops const&))
            &w_t::process, (
          arg_("original_site"), arg_("site_symmetry_ops")))
        .def("process_sites_frac",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            double const&))
              &w_t::process_sites_frac, process_sites_frac_overloads((
          arg_("original_sites"), arg_("min_distance_sym_equiv"))))
        .def("process_sites_frac",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_frac, (
          arg_("original_sites"), arg_("site_symmetry_table")))
        .def("process_sites_cart",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            double const&))
              &w_t::process_sites_cart, process_sites_cart_overloads((
          arg_("original_sites"), arg_("min_distance_sym_equiv"))))
        .def("process_sites_cart",
          (void(w_t::*)(
            af::const_ref<scitbx::vec3<double> > const&,
            sgtbx::site_symmetry_table const&))
              &w_t::process_sites_cart, (
          arg_("original_sites"), arg_("site_symmetry_table")))
        .def("n_sites_in_asu_and_buffer", &w_t::n_sites_in_asu_and_buffer)
        .def("n_sites_in_asu_and_buffer", &w_t::n_sites_in_asu_and_buffer)
        .def("lock", &w_t::lock)
        .def("is_locked", &w_t::is_locked)
        .def("mappings", &w_t::mappings, ccr())
        .def("mapped_sites_min", &w_t::mapped_sites_min, ccr())
        .def("mapped_sites_max", &w_t::mapped_sites_max, ccr())
        .def("mapped_sites_span", &w_t::mapped_sites_span)
        .def("special_op", &w_t::special_op, (arg_("i_seq")), ccr())
        .def("site_symmetry_table", &w_t::site_symmetry_table, rir())
        .def("get_rt_mx", &w_t::get_rt_mx, (arg_("i_seq"), arg_("i_sym")))
        .def("diff_vec", &w_t::diff_vec, (arg_("pair")))
        .def("map_moved_site_to_asu", &w_t::map_moved_site_to_asu, (
          arg_("moved_original_site"), arg_("i_seq"), arg_("i_sym")))
        .def("r_inv_cart", &w_t::r_inv_cart, (arg_("i_seq"), arg_("i_sym")))
        .def("is_simple_interaction", &w_t::is_simple_interaction, (
          arg_("pair")))
        .def("interaction_type_id", &w_t::interaction_type_id, (
          arg_("pair")))
        .def("make_trial_pair", &w_t::make_trial_pair, (
          arg_("i_seq"), arg_("j_seq"), arg_("j_sym")))
        .def("make_pair", &w_t::make_pair, (
          arg_("i_seq"), arg_("j_seq"), arg_("j_sym")))
        .def("find_i_sym", &w_t::find_i_sym, (arg_("i_seq"), arg_("rt_mx")))
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
