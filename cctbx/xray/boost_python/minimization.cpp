#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/xray/minimization.h>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_minimization()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround

    def("minimization_apply_shifts",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group_type const&,
        af::const_ref<scatterer<> > const&,
        scattering_dictionary const&,
        gradient_flags const&,
        af::const_ref<double> const&,
        double const&))
          minimization::apply_shifts, (
      arg_("unit_cell"),
      arg_("space_group_type"),
      arg_("scatterers"),
      arg_("scattering_dict"),
      arg_("gradient_flags"),
      arg_("shifts"),
      arg_("d_min")));

    def("minimization_add_geometry_restraints_site_gradients",
      (void(*)(
        af::const_ref<scatterer<> > const&,
        gradient_flags const&,
        af::ref<double> const&,
        af::const_ref<scitbx::vec3<double> > const&))
          minimization::add_geometry_restraints_site_gradients, (
      arg_("scatterers"),
      arg_("gradient_flags"),
      arg_("xray_gradients"),
      arg_("geometry_restraints_site_gradients")));
  }

}}} // namespace cctbx::xray::boost_python
