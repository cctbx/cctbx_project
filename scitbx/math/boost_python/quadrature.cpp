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

  struct seven_twelve_0120_wrappers
  {
    typedef scitbx::math::quadrature::seven_twelve_0120<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("seven_twelve_0120",no_init)
        .def(init<> () )
        .def("coord", &w_t::coord)
        .def("weight", &w_t::weight)
        ;
    }
  };


  struct five_nine_1001_wrappers
  {
    typedef scitbx::math::quadrature::five_nine_1001<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("five_nine_1001",no_init)
        .def(init<> () )
        .def("coord", &w_t::coord)
        .def("weight", &w_t::weight)
        ;
    }
  };


  struct five_nine_1110_wrappers
  {
    typedef scitbx::math::quadrature::five_nine_1110<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("five_nine_1110",no_init)
        .def(init<> () )
        .def("coord", &w_t::coord)
        .def("weight", &w_t::weight)
        ;
    }
  };


  struct nine_twentyone_1012_wrappers
  {
    typedef scitbx::math::quadrature::nine_twentyone_1012<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("nine_twentyone_1012",no_init)
        .def(init<> () )
        .def("coord", &w_t::coord)
        .def("weight", &w_t::weight)
        ;
    }
  };







} // namespace <anonymous>

namespace boost_python {

  void wrap_quadrature()
  {
    gauss_hermite_engine_wrappers::wrap();
    seven_twelve_0120_wrappers::wrap();
    five_nine_1001_wrappers::wrap();
    five_nine_1110_wrappers::wrap();
    nine_twentyone_1012_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
