#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/xray/minimization.h>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_minimization()
  {
    using namespace boost::python;

    def("minimization_apply_shifts",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        af::const_ref<scatterer<> > const&,
        gradient_flags const&,
        af::const_ref<double> const&))
          minimization::apply_shifts, (
      arg_("unit_cell"),
      arg_("scatterers"),
      arg_("gradient_flags"),
      arg_("shifts")));

    def("minimization_add_site_gradients",
      (void(*)(
        af::const_ref<scatterer<> > const&,
        gradient_flags const&,
        af::ref<double> const&,
        af::const_ref<scitbx::vec3<double> > const&))
          minimization::add_site_gradients, (
      arg_("scatterers"),
      arg_("gradient_flags"),
      arg_("xray_gradients"),
      arg_("site_gradients")));
    def("minimization_extract_site_gradients",
      (af::shared<scitbx::vec3<double> >(*)(
        af::const_ref<scatterer<> > const&,
        gradient_flags const&,
        af::const_ref<double> const&))
          minimization::extract_site_gradients, (
      arg_("scatterers"),
      arg_("gradient_flags"),
      arg_("xray_gradients")));
    def("minimization_add_u_iso_gradients",
      (void(*)(
        af::const_ref<scatterer<> > const&,
        gradient_flags const&,
        af::ref<double> const&,
        af::const_ref<double> const&))
          minimization::add_u_iso_gradients, (
      arg_("scatterers"),
      arg_("gradient_flags"),
      arg_("xray_gradients"),
      arg_("u_iso_gradients")));
  }

}}} // namespace cctbx::xray::boost_python
