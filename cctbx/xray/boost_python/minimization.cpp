#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <cctbx/xray/minimization.h>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_minimization()
  {
    using namespace boost::python;

    def("minimization_apply_shifts",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group_type const&,
        af::const_ref<scatterer<> > const&,
        scattering_dictionary const&,
        gradient_flags const&,
        af::const_ref<double> const&,
        double const&)) minimization::apply_shifts);
  }

}}} // namespace cctbx::xray::boost_python
