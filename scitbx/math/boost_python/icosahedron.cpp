#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/def.hpp>


#include <scitbx/math/icosahedron.h>

namespace scitbx { namespace math {

namespace {

  struct icosahedron_wrappers
  {
    typedef icosahedron<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("icosahedron", no_init)
        .def(init<int>())
        .add_property("level",make_getter(&w_t::level))
        .add_property("sites",make_getter(&w_t::sites,rbv()))
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
