#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/rsr/rsr.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace rsr {
namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;
    class_<manager<> >("manager",
      init<int const&,
           int const&,
           int const&,
           cctbx::xray::scattering_type_registry const&,
           cctbx::uctbx::unit_cell const&,
           af::const_ref<cctbx::xray::scatterer<> > const&,
           optional<double const&,
                    double const&> >((
                                  arg("nx"),
                                  arg("ny"),
                                  arg("nz"),
                                  arg("scattering_type_registry"),
                                  arg("unit_cell"),
                                  arg("scatterers"),
                                  arg("exp_table_one_over_step_size")=-100,
                                  arg("wing_cutoff")=1.e-3)))
      .add_property("density_array", make_getter(&manager<>::density_array, rbv()))
    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::rsr

BOOST_PYTHON_MODULE(mmtbx_rsr_ext)
{
  mmtbx::rsr::init_module();
}
