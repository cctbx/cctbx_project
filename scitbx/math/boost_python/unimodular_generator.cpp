#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/unimodular_generator.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace math {

namespace {

  struct unimodular_generator_wrappers
  {
    typedef unimodular_generator<int> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("unimodular_generator", no_init)
        .def(init<int const&>((arg("range"))))
        .def("at_end", &w_t::at_end)
        .def("next", &w_t::next)
        .def("__next__", &w_t::next)
        .def("count", &w_t::count)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_unimodular_generator()
  {
    unimodular_generator_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
