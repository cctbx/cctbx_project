#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/random.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct mersenne_twister_wrappers
  {
    typedef random::mersenne_twister w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      seed_overloads, seed, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("mersenne_twister", no_init)
        .def(init<optional<unsigned> >((arg_("seed")=0)))
        .def("random_size_t_min", &w_t::random_size_t_min)
        .def("random_size_t_max", &w_t::random_size_t_max)
        .def("seed", &w_t::seed, seed_overloads((arg_("value")=0)))
        .def("random_size_t", (std::size_t(w_t::*)()) &w_t::random_size_t)
        .def("random_size_t",
          (af::shared<std::size_t>(w_t::*)(std::size_t))
            &w_t::random_size_t, (arg_("size")))
        .def("random_size_t",
          (af::shared<std::size_t>(w_t::*)(std::size_t, std::size_t))
            &w_t::random_size_t, (arg_("size"), arg_("modulus")))
        .def("random_double", (double(w_t::*)()) &w_t::random_double)
        .def("random_double", (af::shared<double>(w_t::*)(std::size_t))
          &w_t::random_double, (arg_("size")))
        .def("random_double",
          (af::shared<double>(w_t::*)(std::size_t, double))
            &w_t::random_double, (arg_("size"), arg_("factor")))
        .def("random_permutation", &w_t::random_permutation, (arg_("size")))
        .def("random_double_point_on_sphere",
          &w_t::random_double_point_on_sphere)
        .def("getstate", &w_t::getstate)
        .def("setstate", &w_t::setstate, (arg_("state")))
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_random()
  {
    mersenne_twister_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
