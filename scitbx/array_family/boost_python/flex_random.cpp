#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/random.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct mersenne_twister_wrappers
  {
    typedef random::mersenne_twister w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      using boost::python::arg;
      class_<w_t>("mersenne_twister", no_init)
        .def(init<unsigned>((arg("seed")=0)))
        .def(init<boost_random::mt19937 &>())
        .def("random_size_t_min", &w_t::random_size_t_min)
        .def("random_size_t_max", &w_t::random_size_t_max)
        .def("seed", &w_t::seed, (arg("value")=0))
        .def("random_size_t", (std::size_t(w_t::*)()) &w_t::random_size_t)
        .def("random_size_t",
          (af::shared<std::size_t>(w_t::*)(std::size_t))
            &w_t::random_size_t, (arg("size")))
        .def("random_size_t",
          (af::shared<std::size_t>(w_t::*)(std::size_t, std::size_t))
            &w_t::random_size_t, (arg("size"), arg("modulus")))
        .def("random_double", (double(w_t::*)()) &w_t::random_double)
        .def("random_double", (af::shared<double>(w_t::*)(std::size_t))
          &w_t::random_double, (arg("size")))
        .def("random_double",
          (af::shared<double>(w_t::*)(std::size_t, double))
            &w_t::random_double, (arg("size"), arg("factor")))
        .def("random_bool", &w_t::random_bool, (
          arg("size"), arg("threshold")))
        .def("random_permutation", &w_t::random_permutation, (arg("size")))
        .def("random_double_point_on_sphere",
          &w_t::random_double_point_on_sphere)
        .def("random_double_unit_quaternion",
          &w_t::random_double_unit_quaternion)
        .def("random_double_r3_rotation_matrix",
          &w_t::random_double_r3_rotation_matrix)
        .def("random_double_r3_rotation_matrix_arvo_1992",
          &w_t::random_double_r3_rotation_matrix_arvo_1992)
        .def("random_int_gaussian_distribution",
          (af::shared<int>(w_t::*)(std::size_t,const double&,const double&))
          &w_t::random_int_gaussian_distribution,
          (arg("size"), arg("mu"), arg("sigma"))
          )
        .def("getstate", &w_t::getstate)
        .def("setstate", &w_t::setstate, (arg("state")))
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_random()
  {
    mersenne_twister_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
