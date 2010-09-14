#include <cctbx/boost_python/flex_fwd.h>

#include <smtbx/structure_factors/direct/standard_xray.h>

#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

namespace smtbx { namespace structure_factors { namespace direct {
  namespace boost_python {

    template <class wt>
    struct linearisation_class_ : boost::python::class_<wt>
    {
      typedef typename wt::complex_type complex_type;
      typedef typename wt::float_type float_type;

      static void compute(wt &self, miller::index<> const &h) {
        self.compute(h);
      }

      static complex_type f_calc(wt const &self) {
        return self.f_calc;
      }

      static float_type observable(wt const &self) {
        return self.observable;
      }

      static af::shared<complex_type> grad_f_calc(wt const &self) {
        return self.grad_f_calc.array();
      }

      static af::shared<float_type> grad_observable(wt const &self) {
        return self.grad_observable.array();
      }

      linearisation_class_(std::string const &name)
        : boost::python::class_<wt>(name.c_str(), boost::python::no_init)
      {
        using namespace boost::python;
        (*this)
          .def("compute", compute, args("miller_index"))
          .add_property("f_calc", f_calc)
          .add_property("grad_f_calc", grad_f_calc)
          .add_property("observable", observable)
          .add_property("grad_observable", grad_observable)
          ;
      }
    };

    template <typename FloatType,
              template<typename, bool> class ObservableType,
              template<typename> class ExpI2PiFunctor>
    struct linearisation_wrapper
    {
      typedef FloatType float_type;
      typedef ExpI2PiFunctor<float_type> exp_i_2pi_functor;

      static void wrap_custom_trigo(char const *core_name) {
        using namespace boost::python;
        typedef one_h_linearisation::custom_trigonometry<FloatType,
                                                            true,
                                                            ObservableType,
                                                            ExpI2PiFunctor>
                wt;
        std::string name(core_name);
        name += std::string("_with_custom_trigonometry");
        linearisation_class_<wt>(name)
          .def(init<uctbx::unit_cell const &,
                    sgtbx::space_group const &,
                    af::shared< xray::scatterer<float_type> > const &,
                    xray::scattering_type_registry const &,
                    exp_i_2pi_functor const &>
               ((arg("unit_cell"),
                 arg("space_group"),
                 arg("scatterers"),
                 arg("scattering_type_registry"),
                 arg("exp_i_2pi_functor")))
                [with_custodian_and_ward<1, 2,
                 with_custodian_and_ward<1, 3,
                 with_custodian_and_ward<1, 4,
                 with_custodian_and_ward<1, 5,
                 with_custodian_and_ward<1, 6> > > > >()])
          ;
        }

      static void wrap_std_trigo(char const *core_name) {
        using namespace boost::python;
        typedef one_h_linearisation::std_trigonometry<FloatType,
                                                        true,
                                                        ObservableType>
                wt;
        std::string name(core_name);
        name += std::string("_with_std_trigonometry");
        linearisation_class_<wt>(name)
          .def(init<uctbx::unit_cell const &,
                    sgtbx::space_group const &,
                    af::shared< xray::scatterer<float_type> > const &,
                    xray::scattering_type_registry const &>
               ((arg("unit_cell"),
                 arg("space_group"),
                 arg("scatterers"),
                 arg("scattering_type_registry")))
                [with_custodian_and_ward<1, 2,
                 with_custodian_and_ward<1, 3,
                 with_custodian_and_ward<1, 4,
                 with_custodian_and_ward<1, 5> > > >()])
          ;
      }

      static void wrap(char const *core_name) {
        wrap_custom_trigo(core_name);
        wrap_std_trigo(core_name);
      }
    };


    void wrap_standard_xray() {
      using namespace one_h_linearisation;

      linearisation_wrapper<double, modulus_squared,
                            cctbx::math::cos_sin_table>
        ::wrap("linearisation_of_f_calc_modulus_squared");

      linearisation_wrapper<double, modulus,
                            cctbx::math::cos_sin_table>
        ::wrap("linearisation_of_f_calc_modulus");
    }
  }
}}}
