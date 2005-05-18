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
        .def("random_size_t", &w_t::random_size_t, (arg_("size")))
        .def("random_double", &w_t::random_double, (arg_("size")))
        .def("random_permutation", &w_t::random_permutation, (arg_("size")))
        .def("random_integer", &w_t::random_integer,
              ((arg_("size"),
               arg_("limit")) ))
        .def("seed", &w_t::seed, seed_overloads((arg_("value")=0)))
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_random()
  {
    mersenne_twister_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
