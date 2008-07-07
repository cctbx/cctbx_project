#include <smtbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>

#include <cctbx/xray/scatterer.h>

#include <smtbx/refinement/constraints/special_positions.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

template<typename FloatType, class XrayScattererType>
struct special_positions_wrapper
{
  typedef special_positions<FloatType, XrayScattererType,
                            af::boost_python::flex_1d>
          wt;

  static void wrap(char const *name) {
    using namespace boost::python;

    class_<wt>(name, no_init)
      .def(init<uctbx::unit_cell const &, // 2nd
                sgtbx::site_symmetry_table const &, // 3rd
                af::shared<typename wt::xray_scatterer_type>, // 4th
                typename wt::parameter_map_type const &, // 5th
                af::ref<xray::scatterer_flags> >( // 6th
           (arg_("unit_cell"),
            arg_("site_symmetry_table"),
            arg_("scatterers"),
            arg_("crystallographic_parameter_map"),
            arg_("constraint_flags")))
           [with_custodian_and_ward<1,2,
            with_custodian_and_ward<1,3,
            with_custodian_and_ward<1,5> > >()])
      .add_property("already_constrained", &wt::already_constrained)
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      ;

  }
};

void wrap_special_positions() {
  special_positions_wrapper<double, xray::scatterer<> >::wrap(
    "special_positions");
}

}}}} // smtbx::refinement::constraints::boost_python
