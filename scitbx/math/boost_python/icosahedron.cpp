#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <scitbx/math/icosahedron.h>

namespace scitbx { namespace math {

namespace {

  struct icosahedron_wrappers
  {
    typedef icosahedron<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("icosahedron", no_init)
        .def(init<int>((arg_("level"))))
        .def_readonly("level", &w_t::level)
        .add_property("sites", make_getter(&w_t::sites, rbv()))
      ;
    }
  };

} // namespace <anoymous>

namespace boost_python {

  void wrap_icosahedron()
  {
    icosahedron_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
