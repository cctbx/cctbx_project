#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/lbfgsb.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace lbfgsb { namespace {

  template <typename FloatType>
  struct minimizer_wrappers
  {
    typedef minimizer<FloatType> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      process_overloads, process, 3, 4)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("minimizer", no_init)
        .def(init<int const&,
                  int const&,
                  af::shared<FloatType>,
                  af::shared<FloatType>,
                  af::shared<int>,
                  FloatType const&,
                  FloatType const&,
                  int const&>())
        .def("process", &w_t::process, process_overloads())
        .def("requests_f_and_g", &w_t::requests_f_and_g)
        .def("is_terminated", &w_t::is_terminated)
        .def("task", &w_t::task)
        .def("f_list", &w_t::f_list)
        .def("f", &w_t::f)
        .def("request_restart", &w_t::request_restart)
        .def("request_stop", &w_t::request_stop)
        .def("request_stop_with_restore", &w_t::request_stop_with_restore)
        .def("n", &w_t::n)
        .def("m", &w_t::m)
        .def("l", &w_t::l)
        .def("u", &w_t::u)
        .def("nbd", &w_t::nbd)
        .def("factr", &w_t::factr)
        .def("pgtol", &w_t::pgtol)
        .def("iprint", &w_t::iprint)
        .def(  "initial_x_replaced_by_projection",
          &w_t::initial_x_replaced_by_projection)
        .def("is_constrained", &w_t::is_constrained)
        .def("is_fully_constrained", &w_t::is_fully_constrained)
        .def("n_iteration", &w_t::n_iteration)
        .def("n_fg_evaluations_total", &w_t::n_fg_evaluations_total)
        .def("n_fg_evaluations_iter", &w_t::n_fg_evaluations_iter)
        .def(  "n_intervals_explored_cauchy_search_total",
          &w_t::n_intervals_explored_cauchy_search_total)
        .def(  "n_intervals_explored_cauchy_search_iter",
          &w_t::n_intervals_explored_cauchy_search_iter)
        .def(  "n_skipped_bfgs_updates_total",
          &w_t::n_skipped_bfgs_updates_total)
        .def("n_bfgs_updates_total", &w_t::n_bfgs_updates_total)
        .def(  "subspace_argmin_is_within_box",
          &w_t::subspace_argmin_is_within_box)
        .def("n_free_variables", &w_t::n_free_variables)
        .def("n_active_constraints", &w_t::n_active_constraints)
        .def(  "n_variables_leaving_active_set",
          &w_t::n_variables_leaving_active_set)
        .def(  "n_variables_entering_active_set",
          &w_t::n_variables_entering_active_set)
        .def("theta_bfgs_matrix_current", &w_t::theta_bfgs_matrix_current)
        .def("f_previous_iteration", &w_t::f_previous_iteration)
        .def("floating_point_epsilon", &w_t::floating_point_epsilon)
        .def(  "factr_times_floating_point_epsilon",
          &w_t::factr_times_floating_point_epsilon)
        .def(  "two_norm_line_search_direction_vector",
          &w_t::two_norm_line_search_direction_vector)
        .def(  "two_norm_line_search_direction_vector_sq",
          &w_t::two_norm_line_search_direction_vector_sq)
        .def(  "accumulated_time_cauchy_search",
          &w_t::accumulated_time_cauchy_search)
        .def(  "accumulated_time_subspace_minimization",
          &w_t::accumulated_time_subspace_minimization)
        .def(  "accumulated_time_line_search",
          &w_t::accumulated_time_line_search)
        .def(  "slope_line_search_function_current",
          &w_t::slope_line_search_function_current)
        .def(  "slope_line_search_function_start",
          &w_t::slope_line_search_function_start)
        .def(  "maximum_relative_step_length",
          &w_t::maximum_relative_step_length)
        .def(  "relative_step_length_line_search",
          &w_t::relative_step_length_line_search)
        .def(  "infinity_norm_projected_gradient",
          &w_t::infinity_norm_projected_gradient)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    minimizer_wrappers<double>::wrap();
  }

}}} // namespace scitbx::lbfgs::<anonymous>

BOOST_PYTHON_MODULE(scitbx_lbfgsb_ext)
{
  scitbx::lbfgsb::init_module();
}
