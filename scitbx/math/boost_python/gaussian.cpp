#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/tuple.hpp>
#include <scitbx/math/gaussian/fit.h>
#include <boost_adaptbx/optional_conversions.h>
#include <boost/optional.hpp>

namespace scitbx { namespace math { namespace gaussian {

namespace {

  struct term_wrappers
  {
    typedef term<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("gaussian_term", no_init)
        .def(init<double const&, double const&>((arg("a"), arg("b"))))
        .def_readwrite("a", &w_t::a)
        .def_readwrite("b", &w_t::b)
        .def("at_x_sq", &w_t::at_x_sq, (arg("x_sq")))
        .def("at_x", &w_t::at_x, (arg("x")))
        .def("gradient_dx_at_x", &w_t::gradient_dx_at_x, (arg("x")))
        .def("integral_dx_at_x", &w_t::integral_dx_at_x, (
          arg("x"),
          arg("b_min_for_erf_based_algorithm")=1e-3))
        .def("gradients_d_ab_at_x_sq", &w_t::gradients_d_ab_at_x_sq, (
          arg("x_sq")))
      ;
    }
  };

  struct sum_wrappers : boost::python::pickle_suite
  {
    typedef sum<> w_t;

    static
    boost::python::tuple
    getinitargs(w_t const& w)
    {
      return boost::python::make_tuple(
        w.array_of_a(), w.array_of_b(), w.c(), w.use_c());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("gaussian_sum", no_init)
        .def(init<double const&, optional<bool> >((
          arg("c"), arg("use_c")=true)))
        .def(init<
          af::small<double, w_t::max_n_terms> const&,
          af::small<double, w_t::max_n_terms> const&,
          optional<double const&, bool> >((
            arg("a"), arg("b"), arg("c")=0, arg("use_c")=true)))
        .def(init<
          af::const_ref<double> const&,
          optional<double const&, bool> >((
            arg("ab"), arg("c")=0, arg("use_c")=true)))
        .def(init<sum<> const&>())
        .def("n_terms", &w_t::n_terms)
        .def("array_of_a", &w_t::array_of_a)
        .def("array_of_b", &w_t::array_of_b)
        .def("c", &w_t::c, ccr())
        .def("use_c", &w_t::use_c)
        .def("n_parameters", &w_t::n_parameters)
        .def("parameters", &w_t::parameters)
        .def("at_x_sq", (double(w_t::*)(double const&) const) &w_t::at_x_sq, (
          arg("x_sq")))
        .def("at_x_sq",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
            &w_t::at_x_sq, (
              arg("x_sq")))
        .def("at_x", (double(w_t::*)(double const&) const) &w_t::at_x, (
          arg("x")))
        .def("at_x",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
            &w_t::at_x, (
              arg("x")))
        .def("gradient_dx_at_x", &w_t::gradient_dx_at_x, (arg("x")))
        .def("integral_dx_at_x", &w_t::integral_dx_at_x, (
          arg("x"),
          arg("b_min_for_erf_based_algorithm")=1e-3))
        .def_pickle(sum_wrappers())
      ;
      boost_adaptbx::optional_conversions::to_and_from_python<
        boost::optional<sum<> > >();
    }
  };

  struct fit_wrappers
  {
    typedef fit<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<sum<> > >("gaussian_fit", no_init)
        .def(init<
          af::shared<double> const&,
          af::shared<double> const&,
          af::shared<double> const&,
          sum<double> const&>((
            arg("table_x"),
            arg("table_y"),
            arg("table_sigmas"),
            arg("start"))))
        .def(init<
          af::shared<double> const&,
          sum<double> const&,
          af::shared<double> const&,
          sum<double> const&>((
            arg("table_x"),
            arg("reference"),
            arg("table_sigmas"),
            arg("start"))))
        .def("table_x", &w_t::table_x)
        .def("table_y", &w_t::table_y)
        .def("table_sigmas", &w_t::table_sigmas)
        .def("fitted_values", &w_t::fitted_values)
        .def("differences", &w_t::differences)
        .def("significant_relative_errors", &w_t::significant_relative_errors)
        .def("bound_flags", &w_t::bound_flags, (
          arg("a_bounded"), arg("b_bounded")))
        .def("apply_shifts", &w_t::apply_shifts, (
          arg("shifts"), arg("enforce_positive_b")))
        .def("target_function", &w_t::target_function, (
          arg("power"), arg("use_sigmas"), arg("differences")))
        .def("gradients_d_abc", &w_t::gradients_d_abc, (
          arg("power"), arg("use_sigmas"), arg("differences")))
        .def("gradients_d_shifts", &w_t::gradients_d_shifts, (
          arg("shifts"), arg("gradients_d_abc")))
        .def("least_squares_jacobian_abc", &w_t::least_squares_jacobian_abc)
        .def("least_squares_hessian_abc_as_packed_u",
          &w_t::least_squares_hessian_abc_as_packed_u)
      ;
    }
  };

}} // gaussian::namespace <anonymous>

namespace boost_python {

  void wrap_gaussian()
  {
    gaussian::term_wrappers::wrap();
    gaussian::sum_wrappers::wrap();
    gaussian::fit_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
