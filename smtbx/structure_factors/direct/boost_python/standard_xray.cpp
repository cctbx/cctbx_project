#include <cctbx/boost_python/flex_fwd.h>

#include <smtbx/structure_factors/direct/standard_xray.h>
#include <smtbx/structure_factors/direct/table_based.h>

#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/manage_new_object.hpp>

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
      typedef one_scatterer_one_h::scatterer_contribution<float_type>
        scatterer_contribution_type;
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
                    exp_i_2pi_functor const &,
                    scatterer_contribution_type *,
                    bool>
               ((arg("unit_cell"),
                 arg("space_group"),
                 arg("scatterers"),
                 arg("exp_i_2pi_functor"),
                 arg("scatter_contribution"),
                 arg("own_scatterer_contribution") = false))
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
                    scatterer_contribution_type *,
                    bool>
               ((arg("unit_cell"),
                 arg("space_group"),
                 arg("scatterers"),
                 arg("scatter_contribution"),
                 arg("own_scatterer_contribution") = false))
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

    template <typename FloatType>
    struct scatterer_contribution_wrapper {
      typedef one_scatterer_one_h::scatterer_contribution<FloatType> wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt, boost::noncopyable>("scatterer_contribution", no_init)
          .def("get", &wt::get,
            (arg("scatterer_index"),
             arg("h")))
          .def("at_d_star_sq", &wt::at_d_star_sq,
            (arg("d_start_sq")), return_internal_reference<>())
          ;
      }
    };

    template <typename FloatType>
    struct isotropic_scatterer_contribution_wrapper {
      typedef one_scatterer_one_h::isotropic_scatterer_contribution<FloatType>
        wt;
      typedef one_scatterer_one_h::scatterer_contribution<FloatType>
        scatterer_contribution_type;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
          bases<scatterer_contribution_type> >
          ("isotropic_scatterer_contribution", no_init)
          .def(init<af::shared< xray::scatterer<FloatType> > const &,
            xray::scattering_type_registry const &>(
            (arg("scatterers"),
             arg("scattering_type_registry"))))
          .def(init<af::shared< xray::scatterer<FloatType> > const &,
            xray::scattering_type_registry const &,
            uctbx::unit_cell const &,
            cctbx::xray::observations<FloatType> const &>(
            (arg("scatterers"),
              arg("scattering_type_registry"),
              arg("unit_cell"),
              arg("reflections"))))
          ;
      }
    };

    template <typename FloatType>
    struct table_based_wrapper {
      typedef table_based::builder<FloatType> wt;
      typedef one_scatterer_one_h::scatterer_contribution<FloatType>
        scatterer_contribution_type;

      static void wrap() {
        using namespace boost::python;
        class_<wt, boost::noncopyable>("table_based_scatterer_contribution", no_init)
          .def("build", &wt::build,
            (arg("unit_cell"),
              arg("scatterers"),
              arg("file_name"),
              arg("space_group"),
              arg("anomalous_flag")),
            return_value_policy<manage_new_object>())
          .staticmethod("build")
          .def("build_lookup_based_for_tests", &wt::build_lookup_based_for_tests,
          (arg("unit_cell"),
            arg("space_group"),
            arg("scatterers"),
            arg("scattering_type_registry"),
            arg("indices")),
            return_value_policy<manage_new_object>())
          .staticmethod("build_lookup_based_for_tests")
          ;
      }
    };

    void wrap_standard_xray() {
      fc_for_one_h_wrapper<double, one_h::modulus_squared,
                           cctbx::math::cos_sin_table>::wrap();

      fc_for_one_h_wrapper<double, one_h::modulus,
                           cctbx::math::cos_sin_table>::wrap();

      scatterer_contribution_wrapper<double>::wrap();
      isotropic_scatterer_contribution_wrapper<double>::wrap();
      table_based_wrapper<double>::wrap();
    }
  }
}}}
