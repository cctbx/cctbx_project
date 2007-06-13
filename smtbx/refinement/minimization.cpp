#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <smtbx/refinement/minimization.h>

namespace smtbx { namespace refinement { namespace boost_python {

  struct apply_special_position_constrained_shifts_wrappers
  {
    typedef minimization::apply_special_position_constrained_shifts<
              scatterer<>, double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(
        "apply_special_position_constrained_shifts",
        no_init
      )
        .def(init<
          uctbx::unit_cell const&,
          const sgtbx::site_symmetry_table&,
          af::const_ref<scatterer<> > const&,
          af::const_ref<double> const& >((
            arg_("unit_cell"),
            arg_("site_symmetry_table"),
            arg_("scatterers"),
            arg_("shifts"))))
        .add_property("shifted_scatterers",
          make_getter(&w_t::shifted_scatterers, rbv()))
        .add_property("u_iso_refinable_params",
          make_getter(&w_t::u_iso_refinable_params, rbv()))
      ;
    }
  };

  struct special_position_constrained_gradients_wrappers
  {
    typedef minimization::special_position_constrained_gradients<
              scatterer<>, double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(
        "special_position_constrained_gradients",
        no_init
      )
        .def(init<
          uctbx::unit_cell const&,
          const sgtbx::site_symmetry_table&,
          af::const_ref<scatterer<> > const&,
          af::ref<double> const& >((
            arg_("unit_cell"),
            arg_("site_symmetry_table"),
            arg_("scatterers"),
            arg_("xray_gradients"))))
        .add_property("reduced_gradients",
          make_getter(&w_t::reduced_gradients, rbv()))
      ;
    }
  };

  void wrap_minimization()
  {
    using namespace boost::python;

    apply_special_position_constrained_shifts_wrappers::wrap();
    special_position_constrained_gradients_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
