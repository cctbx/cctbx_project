#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <iotbx/mtz/crystal.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace iotbx { namespace mtz {
namespace {

  struct crystal_wrappers
  {
    typedef crystal w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("crystal", no_init)
        .def(init<object const&, int>((arg_("object"), arg_("i_crystal"))))
        .def("object", &w_t::object)
        .def("i_crystal", &w_t::i_crystal)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_crystal");
      }
    }
  };

  void
  wrap_all()
  {
    crystal_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_crystal() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
