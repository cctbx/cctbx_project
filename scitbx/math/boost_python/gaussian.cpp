#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/math/gaussian/fit.h>

namespace scitbx { namespace math { namespace gaussian {

namespace {

  struct term_wrappers
  {
    typedef term<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      integral_dx_at_x_overloads, integral_dx_at_x, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("gaussian_term", no_init)
        .def(init<double const&, double const&>())
        .def_readwrite("a", &w_t::a)
        .def_readwrite("b", &w_t::b)
        .def("at_x_sq", &w_t::at_x_sq)
        .def("at_x", &w_t::at_x)
        .def("gradient_dx_at_x", &w_t::gradient_dx_at_x)
        .def("integral_dx_at_x", &w_t::integral_dx_at_x,
          integral_dx_at_x_overloads())
        .def("gradients_d_ab_at_x_sq", &w_t::gradients_d_ab_at_x_sq)
      ;
    }
  };

  struct sum_wrappers : boost::python::pickle_suite
  {
    typedef sum<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      integral_dx_at_x_overloads, integral_dx_at_x, 1, 2)

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
        .def(init<double const&, optional<bool> >())
        .def(init<af::small<double, w_t::max_n_terms> const&,
                  af::small<double, w_t::max_n_terms> const&,
                  optional<double const&, bool> >())
        .def(init<af::const_ref<double> const&,
                  optional<double const&, bool> >())
        .def(init<sum<> const&>())
        .def("n_terms", &w_t::n_terms)
        .def("array_of_a", &w_t::array_of_a)
        .def("array_of_b", &w_t::array_of_b)
        .def("c", &w_t::c, ccr())
        .def("use_c", &w_t::use_c)
        .def("n_parameters", &w_t::n_parameters)
        .def("at_x_sq", (double(w_t::*)(double const&) const) &w_t::at_x_sq)
        .def("at_x_sq",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
          &w_t::at_x_sq)
        .def("at_x", (double(w_t::*)(double const&) const) &w_t::at_x)
        .def("at_x",
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const)
          &w_t::at_x)
        .def("gradient_dx_at_x", &w_t::gradient_dx_at_x)
        .def("integral_dx_at_x", &w_t::integral_dx_at_x,
          integral_dx_at_x_overloads())
        .def_pickle(sum_wrappers())
      ;
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
        .def(init<af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  sum<double> const&>())
        .def(init<af::shared<double> const&,
                  sum<double> const&,
                  af::shared<double> const&,
                  sum<double> const&>())
        .def("table_x", &w_t::table_x)
        .def("table_y", &w_t::table_y)
        .def("table_sigmas", &w_t::table_sigmas)
        .def("fitted_values", &w_t::fitted_values)
        .def("differences", &w_t::differences)
        .def("significant_relative_errors", &w_t::significant_relative_errors)
        .def("bound_flags", &w_t::bound_flags)
        .def("apply_shifts", &w_t::apply_shifts)
        .def("target_function", &w_t::target_function)
        .def("gradients_d_abc", &w_t::gradients_d_abc)
        .def("gradients_d_shifts", &w_t::gradients_d_shifts)
      ;
    }
  };

}} // gaussian::namespace <anoymous>

namespace boost_python {

  void wrap_gaussian()
  {
    gaussian::term_wrappers::wrap();
    gaussian::sum_wrappers::wrap();
    gaussian::fit_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
