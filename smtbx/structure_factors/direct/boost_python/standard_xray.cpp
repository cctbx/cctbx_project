#include <cctbx/boost_python/flex_fwd.h>

#include <smtbx/structure_factors/direct/standard_xray.h>

#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

namespace smtbx { namespace structure_factors { namespace direct {
  namespace boost_python {

    template <class wt>
    struct fc_for_one_h_class : boost::python::class_<wt>
    {
      typedef typename wt::complex_type complex_type;
      typedef typename wt::float_type float_type;

      static void evaluate(wt &self, miller::index<> const &h) {
        self.evaluate(h);
      }

      static void linearise(wt &self, miller::index<> const &h) {
        self.linearise(h);
      }

      static complex_type f_calc(wt const &self) {
        return self.f_calc;
      }

      static float_type observable(wt const &self) {
        return self.observable;
      }

      static boost::python::object grad_f_calc(wt const &self) {
        using namespace boost::python;
        return self.is_linearisation() ? object(self.grad_f_calc.array())
                                       : object();
      }

      static boost::python::object grad_observable(wt const &self) {
        using namespace boost::python;
        return self.is_linearisation() ? object(self.grad_observable.array())
                                       : object();
      }

      fc_for_one_h_class(std::string const &name)
        : boost::python::class_<wt>(name.c_str(), boost::python::no_init)
      {
        using namespace boost::python;
        (*this)
          .def("evaluate" , evaluate , args("miller_index"))
          .def("linearise", linearise, args("miller_index"))
          .add_property("f_calc", f_calc)
          .add_property("observable", observable)
          .add_property("grad_f_calc", grad_f_calc)
          .add_property("grad_observable", grad_observable)
          ;
      }
    };

    template <template<typename> class ObservableType>
    struct observable_traits
    {};

    template <>
    struct observable_traits<one_h::modulus_squared>
    {
      static char const *name() { return "modulus_squared"; }
    };

    template <>
    struct observable_traits<one_h::modulus>
    {
      static char const *name() { return "modulus"; }
    };

    template <typename FloatType,
              template<typename> class ObservableType,
              template<typename> class ExpI2PiFunctor>
    struct fc_for_one_h_wrapper
    {
      typedef FloatType float_type;
      typedef ExpI2PiFunctor<float_type> exp_i_2pi_functor;

      static void wrap_custom_trigo(char const *core_name) {
        using namespace boost::python;
        typedef one_h::custom_trigonometry<FloatType,
                                           ObservableType,
                                           ExpI2PiFunctor>
                wt;
        std::string name(core_name);
        name += std::string("_with_custom_trigonometry");
        fc_for_one_h_class<wt>(name)
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
        typedef one_h::std_trigonometry<FloatType,
                                        ObservableType>
                wt;
        std::string name(core_name);
        name += std::string("_with_std_trigonometry");
        fc_for_one_h_class<wt>(name)
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

      static void wrap() {
        std::string core_name = "f_calc_";
        core_name += observable_traits<ObservableType>::name();
        wrap_custom_trigo(core_name.c_str());
        wrap_std_trigo(core_name.c_str());
      }
    };


    void wrap_standard_xray() {
      fc_for_one_h_wrapper<double, one_h::modulus_squared,
                           cctbx::math::cos_sin_table>::wrap();

      fc_for_one_h_wrapper<double, one_h::modulus,
                           cctbx::math::cos_sin_table>::wrap();
    }
  }
}}}
