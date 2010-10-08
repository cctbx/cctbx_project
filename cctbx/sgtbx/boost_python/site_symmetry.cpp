#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/array_family/versa_matrix.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct site_constraints_wrappers
  {
    typedef site_constraints<> w_t;

    static
    af::versa<int, af::c_grid<2> >
    row_echelon_form_as_versa(w_t const& self)
    {
      return af::mat_const_ref_as_versa(self.row_echelon_form());
    }

    static
    af::versa<double, af::c_grid<2> >
    gradient_sum_matrix_as_versa(w_t const& self)
    {
      return af::mat_const_ref_as_versa(self.gradient_sum_matrix());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("site_constraints", no_init)
        .def_readonly("row_echelon_lcm", &w_t::row_echelon_lcm)
        .def("row_echelon_form", row_echelon_form_as_versa)
        .add_property("row_echelon_constants",
          make_getter(&w_t::row_echelon_constants, rbv()))
        .add_property("independent_indices",
          make_getter(&w_t::independent_indices, rbv()))
        .def("n_independent_params", &w_t::n_independent_params)
        .def("n_dependent_params", &w_t::n_dependent_params)
        .def("independent_params", &w_t::independent_params, (
          arg("all_params")))
        .def("all_params", &w_t::all_params, (
          arg("independent_params")))
        .def("all_shifts", &w_t::all_shifts, (
          arg("independent_shifts")))
        .def("gradient_sum_matrix", gradient_sum_matrix_as_versa)
        .def("independent_gradients",
             (af::small<double, 3> (w_t::*)(af::const_ref<double> const&) const)
             &w_t::independent_gradients,
             (arg("all_gradients")))
        .def("independent_gradients",
             (af::small<double, 3> (w_t::*)(scitbx::vec3<double> const&) const)
             &w_t::independent_gradients,
             (arg("all_gradients")))
        .def("independent_curvatures", &w_t::independent_curvatures, (
          arg("all_curvatures")))
      ;
    }
  };

  struct site_symmetry_ops_wrappers
  {
    typedef site_symmetry_ops w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("site_symmetry_ops", no_init)
        .enable_pickling()
        .def(init<int, rt_mx const&, af::shared<rt_mx> const&>(( // for pickle
          arg("multiplicity"),
          arg("special_op"),
          arg("matrices"))))
        .def("multiplicity", &w_t::multiplicity)
        .def("special_op", &w_t::special_op, ccr())
        .def("matrices", &w_t::matrices, ccr())
        .def("n_matrices", &w_t::n_matrices)
        .def("is_point_group_1", &w_t::is_point_group_1)
        .def("__contains__", &w_t::contains)
        .def("__eq__", &w_t::operator==)
        .def("is_compatible_u_star",
          (bool(w_t::*)(scitbx::sym_mat3<double> const&, double) const)
            &w_t::is_compatible_u_star, (
            arg("u_star"), arg("tolerance")=1e-6))
        .def("average_u_star",
          (scitbx::sym_mat3<double>
            (w_t::*)(scitbx::sym_mat3<double> const&) const)
          &w_t::average_u_star, (arg("u_star")))
        .def("change_basis", &w_t::change_basis, (arg("cb_op")))
        .def("site_constraints", &w_t::site_constraints, rir())
        .def("adp_constraints", &w_t::adp_constraints, rir())
        .def("cartesian_adp_constraints", &w_t::cartesian_adp_constraints, (
          arg("unit_cell"), arg("unit_cell_has_changed")=false), rir())
      ;
    }
  };

  struct site_symmetry_wrappers
  {
    typedef site_symmetry w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t, bases<site_symmetry_ops> >("site_symmetry", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  fractional<double> const&,
                  optional<double, bool> >(
          (arg("unit_cell"),
           arg("space_group"),
           arg("original_site"),
           arg("min_distance_sym_equiv")=0.5,
           arg("assert_min_distance_sym_equiv")=true)))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("original_site", &w_t::original_site, ccr())
        .def("min_distance_sym_equiv", &w_t::min_distance_sym_equiv)
        .def("exact_site", &w_t::exact_site, ccr())
        .def("distance_moved", &w_t::distance_moved)
        .def("shortest_distance", &w_t::shortest_distance)
        .def("check_min_distance_sym_equiv",&w_t::check_min_distance_sym_equiv)
        .def("multiplicity", &w_t::multiplicity)
        .def("point_group_type", &w_t::point_group_type)
        .def("unique_ops", &w_t::unique_ops)
      ;
    }
  };

  struct site_symmetry_table_wrappers
  {
    typedef site_symmetry_table w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      object none;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("site_symmetry_table")
        .enable_pickling()
        .def(init<
          af::shared<std::size_t> const&,
          af::shared<site_symmetry_ops> const&,
          af::shared<std::size_t> const&>(( // for pickle
            arg("indices"),
            arg("table"),
            arg("special_position_indices"))))
        .def("process",
          (void(w_t::*)(std::size_t, site_symmetry_ops const&))
            &w_t::process, (
              arg("insert_at_index"),
              arg("site_symmetry_ops")))
        .def("process",
          (void(w_t::*)(site_symmetry_ops const&))
            &w_t::process,
          (arg("site_symmetry_ops")))
        .def("process",
          (void(w_t::*)(
            uctbx::unit_cell const&,
            space_group const&,
            af::const_ref<scitbx::vec3<double> > const&,
            af::const_ref<bool> const&,
            double,
            bool)) &w_t::process, (
          arg("unit_cell"),
          arg("space_group"),
          arg("original_sites_frac"),
          arg("unconditional_general_position_flags")=none,
          arg("min_distance_sym_equiv")=0.5,
          arg("assert_min_distance_sym_equiv")=true))
        .def("is_special_position", &w_t::is_special_position, (arg("i_seq")))
        .def("get", &w_t::get, (arg("i_seq")), rir())
        .def("n_special_positions", &w_t::n_special_positions)
        .def("special_position_indices",&w_t::special_position_indices,ccr())
        .def("n_unique", &w_t::n_unique)
        .def("indices", &w_t::indices, ccr())
        .def("table", &w_t::table, ccr())
        .def("reserve", &w_t::reserve)
        .def("deep_copy", &w_t::deep_copy)
        .def("change_basis", &w_t::change_basis, (arg("cb_op")))
        .def("select",
          (site_symmetry_table(w_t::*)(
            af::const_ref<std::size_t> const& ) const)
              &w_t::select,
          (arg("selection")))
        .def("select",
          (site_symmetry_table(w_t::*)(
            af::const_ref<bool> const& ) const)
              &w_t::select,
          (arg("selection")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_site_symmetry()
  {
    site_constraints_wrappers::wrap();
    site_symmetry_ops_wrappers::wrap();
    site_symmetry_wrappers::wrap();
    site_symmetry_table_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
