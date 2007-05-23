#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <smtbx/refinement/minimization.h>

namespace smtbx { namespace refinement { namespace boost_python {

  struct apply_shifts_wrappers
  {
    typedef minimization::apply_shifts<scatterer<>, double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("minimization_apply_shifts", no_init)
        .def(init<
          uctbx::unit_cell const&,
          af::const_ref<scatterer<> > const&,
          af::const_ref<double> const& >((
            arg_("unit_cell"),
            arg_("scatterers"),
            arg_("shifts"))))
        .add_property("shifted_scatterers",
          make_getter(&w_t::shifted_scatterers, rbv()))
        .add_property("u_iso_refinable_params",
          make_getter(&w_t::u_iso_refinable_params, rbv()))
      ;
    }
  };
  
  struct reduce_gradient_as_per_special_position_constraints_wrappers
  {
    typedef minimization::reduce_gradient_as_per_special_position_constraints<
              scatterer<>, double> w_t;
    
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("reduce_gradient_as_per_special_position_constraints_wrappers",
                  no_init)
        .def(init<
          uctbx::unit_cell const&,
          af::const_ref<scatterer<> > const&,
          af::ref<double> const& >((
            arg_("unit_cell"),
            arg_("scatterers"),
            arg_("shifts"))))
        .add_property("reduced_gradients",
          make_getter(&w_t::reduced_gradients, rbv()))
      ;
    }
  };

  void wrap_minimization()
  {
    using namespace boost::python;

    apply_shifts_wrappers::wrap();
    reduce_gradient_as_per_special_position_constraints_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
