#include <cctbx/boost_python/flex_fwd.h>

#include <smtbx/structure_factors/direct/standard_xray.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace smtbx { namespace structure_factors { namespace direct {
  namespace boost_python {

    template <typename FloatType,
              template<typename, bool> class ObservableType>
    struct linearisation_wrapper
    {
      typedef one_h_linearisation::base<FloatType, true, ObservableType> wt;
      typedef typename wt::float_type float_type;

      static void wrap(char const *name) {
        using namespace boost::python;
        return_value_policy<return_by_value> rbv;

        typedef void (wt::*compute_with_exact_exp_i_2pi_ptr)(
                  miller::index<> const &);
        typedef void (wt::*compute_with_tabulated_exp_i_2pi_ptr)(
                  miller::index<> const &,
                  cctbx::math::cos_sin_table<float_type> const &);
        typedef void (wt::*compute_with_interpolated_tabulated_exp_i_2pi_ptr)(
                  miller::index<> const &,
                  cctbx::math::cos_sin_lin_interp_table<float_type> const &);

        class_<wt>(name, no_init)
          .def(init<std::size_t,
                    uctbx::unit_cell const &,
                    sgtbx::space_group const &,
                    af::shared< xray::scatterer<float_type> > const &,
                    xray::scattering_type_registry const &>
               ((arg("n_parameters"),
                 arg("unit_cell"),
                 arg("space_group"),
                 arg("scatterers"),
                 arg("scattering_type_registry"))))
          .def("compute",
               (compute_with_exact_exp_i_2pi_ptr) &wt::compute,
               args("miller_index"))
          .def("compute",
               (compute_with_tabulated_exp_i_2pi_ptr) &wt::compute,
               args("miller_index", "tabulated_exp_i_2pi"))
          .def("compute",
               (compute_with_interpolated_tabulated_exp_i_2pi_ptr) &wt::compute,
               args("miller_index", "interpolated_tabulated_exp_i_2pi"))
          .def_readonly("f_calc", &wt::f_calc)
          .add_property("grad_f_calc",
                        make_getter(&wt::grad_f_calc, rbv))
          .def_readonly("observable", &wt::observable)
          .add_property("grad_observable",
                        make_getter(&wt::grad_observable, rbv))
          ;
      }
    };


    void wrap_standard_xray() {
      using namespace one_h_linearisation;
      linearisation_wrapper<double, modulus_squared>::wrap(
        "linearisation_of_f_calc_modulus_squared");

      linearisation_wrapper<double, modulus>::wrap(
        "linearisation_of_f_calc_modulus");
    }
  }
}}}
