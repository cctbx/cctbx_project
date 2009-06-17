#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/curvatures_simple.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace curvatures_simple {
namespace boost_python {

namespace {

  struct grads_and_curvs_wrappers
  {
    typedef grads_and_curvs<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>(
          "structure_factors_curvatures_simple_grads_and_curvs", no_init)
        .def(init<cctbx::miller::index<> const&>((arg_("hkl"))))
        .def("compute", &w_t::compute<scatterer<> >, (
          arg_("space_group"),
          arg_("hkl"),
          arg_("d_star_sq"),
          arg_("scatterers"),
          arg_("scattering_type_registry"),
          arg_("site_symmetry_table"),
          arg_("i_scatterer")))
        .def("copy_gradients", &w_t::copy_gradients)
        .def("copy_curvatures", &w_t::copy_curvatures)
      ;
    }
  };

} // namespace <anoymous>

}}} // namespace structure_factors::curvatures_simple::boost_python

namespace boost_python {

  void wrap_curvatures_simple()
  {
    structure_factors::curvatures_simple::boost_python
      ::grads_and_curvs_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
