#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/dynamics/dynamics.h>

namespace mmtbx { namespace dynamics {
namespace {

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_;
    class_<kinetic_energy_and_temperature>("kinetic_energy_and_temperature",
                                           init<af::shared<vec3<double> > const&,
                                                af::shared<double> const&>())

      .def("kinetic_energy", &kinetic_energy_and_temperature::kinetic_energy)
      .def("temperature", &kinetic_energy_and_temperature::temperature)
    ;
    def("vxyz_at_t_plus_dt_over_2",vxyz_at_t_plus_dt_over_2)
   ;
  }

} // namespace <anonymous>
}} // namespace mmtbx::dynamics

BOOST_PYTHON_MODULE(mmtbx_dynamics_ext)
{
  mmtbx::dynamics::init_module();
}
