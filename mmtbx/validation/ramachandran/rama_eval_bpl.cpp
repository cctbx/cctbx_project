#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/validation/ramachandran/rama_eval.h>

namespace mmtbx { namespace validation { namespace ramachandran {
namespace {

  void init_module()
  {
    using namespace boost::python;

    class_<rama_eval>("rama_eval", no_init)

      .def(init<>())
      .def("get_value", &rama_eval::get_value,
        (arg("rama_class"), arg("phi"), arg("psi")))
    ;
  }

} // namespace <anonymous>
}}} // namespace mmtbx::validation::ramachandran

BOOST_PYTHON_MODULE(mmtbx_validation_ramachandran_ext)
{
  mmtbx::validation::ramachandran::init_module();
}
