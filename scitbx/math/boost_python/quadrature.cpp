#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <scitbx/math/quadrature.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace scitbx { namespace math {

namespace {

  struct gauss_hermite_engine_wrappers
  {
    typedef scitbx::math::quadrature::gauss_hermite_engine<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("gauss_hermite_engine",no_init)
        .def(init<int const& > ((arg_("n_points"))))
        .def("f", &w_t::f)
        .def("refine", &w_t::refine)
        .def("x", &w_t::x)
        .def("w", &w_t::w)
        .def("w_exp_x_squared", &w_t::w_exp_x_squared)
        ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_gh()
  {
    gauss_hermite_engine_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
