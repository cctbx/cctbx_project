#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/tensor_rank_2.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct tensor_rank_2_constraints_wrappers
  {
    typedef tensor_rank_2::constraints<> w_t;

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
      class_<w_t>("tensor_rank_2_constraints", no_init)
        .def(init<sgtbx::space_group const&, bool>((
          arg("space_group"),
          arg("reciprocal_space"))))
        .def(init<af::shared<rt_mx> const&, std::size_t, bool>((
          arg("symmetry_matrices"),
          arg("i_first_matrix_to_use"),
          arg("reciprocal_space"))))
        .def("row_echelon_form", row_echelon_form_as_versa)
        .add_property("independent_indices",
          make_getter(&w_t::independent_indices, rbv()))
        .def("gradient_sum_matrix", gradient_sum_matrix_as_versa)
        .def("n_independent_params", &w_t::n_independent_params)
        .def("n_dependent_params", &w_t::n_dependent_params)
        .def("independent_params", &w_t::independent_params,
          (arg("all_params")))
        .def("all_params", &w_t::all_params, (arg("independent_params")))
        .def("independent_gradients", &w_t::independent_gradients,
          (arg("all_gradients")))
        .def("independent_curvatures", &w_t::independent_curvatures,
          (arg("all_curvatures")))
      ;
    }
  };


  struct tensor_rank_2_cartesian_constraints_wrapper
  {
    typedef tensor_rank_2::cartesian_constraints<> wt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt>("tensor_rank_2_cartesian_constraints", no_init)
        .def(init<uctbx::unit_cell const&, sgtbx::space_group const&>((
          arg("unit_cell"), arg("space_group"))))
        .def("n_independent_params", &wt::n_independent_params)
        .def("all_params", &wt::all_params,
          arg("independent_params"))
        .def("independent_params", &wt::independent_params,
             arg("all_params"))
        .def("independent_gradients", &wt::independent_gradients,
          arg("all_gradients"))
        .def("jacobian", &wt::jacobian, rir)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_tensor_rank_2()
  {
    tensor_rank_2_constraints_wrappers::wrap();
    tensor_rank_2_cartesian_constraints_wrapper::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
