#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/math/gaussian/sum.h>

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
      ;
    }
  };

  struct sum_wrappers
  {
    typedef sum<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      integral_dx_at_x_overloads, integral_dx_at_x, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("gaussian_sum", no_init)
        .def(init<double const&>())
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  optional<double const&, bool> >())
        .def("n_terms", &w_t::n_terms)
        .def("array_of_a", &w_t::array_of_a)
        .def("array_of_b", &w_t::array_of_b)
        .def("c", &w_t::c, ccr())
        .def("use_c", &w_t::use_c)
        .def("all_zero", &w_t::all_zero)
        .def("at_x_sq", &w_t::at_x_sq)
        .def("at_x", &w_t::at_x)
        .def("gradient_dx_at_x", &w_t::gradient_dx_at_x)
        .def("integral_dx_at_x", &w_t::integral_dx_at_x,
          integral_dx_at_x_overloads())
      ;
    }
  };

}} // gaussian::namespace <anoymous>

namespace boost_python {

  void wrap_gaussian()
  {
    gaussian::term_wrappers::wrap();
    gaussian::sum_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
