#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/tensor_rank_2.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct tensor_rank_2_constraints_wrappers
  {
    typedef tensor_rank_2::constraints<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("tensor_rank_2_constraints", no_init)
        .def(init<sgtbx::space_group const&, bool, optional<bool> >((
          arg_("space_group"),
          arg_("reciprocal_space"),
          arg_("initialize_gradient_handling")=false)))
        .def(init<sgtbx::site_symmetry_ops const&, bool, optional<bool> >((
          arg_("site_symmetry_ops"),
          arg_("reciprocal_space"),
          arg_("initialize_gradient_handling")=false)))
        .def(init<
            af::shared<rt_mx> const&, std::size_t, bool, optional<bool> >((
          arg_("symmetry_matrices"),
          arg_("i_first_matrix_to_use"),
          arg_("reciprocal_space"),
          arg_("initialize_gradient_handling")=false)))
        .add_property("row_echelon_form",
          make_getter(&w_t::row_echelon_form_memory, rbv()))
        .add_property("independent_indices",
          make_getter(&w_t::independent_indices, rbv()))
        .add_property("gradient_sum_coeffs",
          make_getter(&w_t::gradient_sum_coeffs, rbv()))
        .def("n_independent_params", &w_t::n_independent_params)
        .def("n_dependent_params", &w_t::n_dependent_params)
        .def("independent_params", &w_t::independent_params,
          (arg_("all_params")))
        .def("all_params", &w_t::all_params, (arg_("independent_params")))
        .def("sym_gradients", &w_t::sym_gradients, (arg_("asu_gradients")))
        .def("independent_gradients", &w_t::independent_gradients,
          (arg_("all_gradients")))
        .def("independent_curvatures", &w_t::independent_curvatures,
          (arg_("all_curvatures")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_tensor_rank_2()
  {
    tensor_rank_2_constraints_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
