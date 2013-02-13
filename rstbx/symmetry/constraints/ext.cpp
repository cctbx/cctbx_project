#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>

#include <rstbx/symmetry/constraints/a_g_conversion.h>

#include <vector>
#include <map>

using namespace boost::python;
namespace rstbx{
namespace boost_python { namespace {

  void
  symmetry_constraints_init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    class_<rstbx::symmetry::AG>("AGconvert", init<>())
      .def("forward", &rstbx::symmetry::AG::forward)
      .def("setAngles", &rstbx::symmetry::AG::setAngles,
           (arg_("phi"),arg_("psi"),arg_("theta"))
          )
      .def("validate_and_setG", &rstbx::symmetry::AG::validate_and_setG)
      .def("back_as_orientation", &rstbx::symmetry::AG::back_as_orientation)
      .def("back", &rstbx::symmetry::AG::back)
      .add_property("G",make_getter(&rstbx::symmetry::AG::G, rbv()),
                        make_setter(&rstbx::symmetry::AG::G, dcp())
          )
      .add_property("phi",make_getter(&rstbx::symmetry::AG::phi, rbv()))
      .add_property("psi",make_getter(&rstbx::symmetry::AG::psi, rbv()))
      .add_property("theta",make_getter(&rstbx::symmetry::AG::theta, rbv()))
    ;

  }

}
}} // namespace rstbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(rstbx_symmetry_constraints_ext)
{
  rstbx::boost_python::symmetry_constraints_init_module();

}
