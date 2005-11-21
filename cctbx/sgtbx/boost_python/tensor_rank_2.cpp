#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/tensor_rank_2.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

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
          arg_("space_group"),
          arg_("reciprocal_space"))))
        .def(init<sgtbx::site_symmetry_ops const&, bool>((
          arg_("site_symmetry_ops"),
          arg_("reciprocal_space"))))
        .def(init<af::shared<rt_mx> const&, std::size_t, bool>((
          arg_("symmetry_matrices"),
          arg_("i_first_matrix_to_use"),
          arg_("reciprocal_space"))))
        .def("row_echelon_form", row_echelon_form_as_versa)
        .add_property("independent_indices",
          make_getter(&w_t::independent_indices, rbv()))
        .def("gradient_sum_matrix", gradient_sum_matrix_as_versa)
        .def("n_independent_params", &w_t::n_independent_params)
        .def("n_dependent_params", &w_t::n_dependent_params)
        .def("independent_params", &w_t::independent_params,
          (arg_("all_params")))
        .def("all_params", &w_t::all_params, (arg_("independent_params")))
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
