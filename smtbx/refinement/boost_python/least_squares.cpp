#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/weighting_schemes.h>


namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

  template <typename FloatType,
            class NormalEquations,
            class OneMillerIndexLinearisation>
  struct normal_equation_building
  {
    typedef build_normal_equations<FloatType> wt;

    template <template<typename> class WeightingSchemeType>
    static void def_init(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass.def(init<NormalEquations &,                    // normal_equations
                cctbx::xray::observations<FloatType> const &,         // miller_indices+data+sigmas
                af::const_ref<std::complex<FloatType> > const &,      // f_mask
                WeightingSchemeType<FloatType> const &,               // weighting_scheme
                boost::optional<FloatType>,                           // scale_factor
                OneMillerIndexLinearisation &,                        // f_calc_function
                scitbx::sparse::matrix<FloatType> const &,            // jacobian_transpose_matching_grad_fc
                cctbx::xray::extinction_correction<FloatType> &,      // extinction
                optional<bool> >                                      // objective_only=false
                ());
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt> klass(name, no_init);
        def_init<mainstream_shelx_weighting>(klass);
        def_init<unit_weighting            >(klass);
        def_init<sigma_weighting           >(klass);
      klass
        .def("observables", &wt::observables)
        .def("f_calc", &wt::f_calc)
        .def("weights", &wt::weights)
        ;
    }
  };

  void wrap_least_squares() {
    using namespace boost::python;

    normal_equation_building<
      double,
      lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<double>,
      structure_factors::direct::one_h::std_trigonometry<
        double,
        structure_factors::direct::one_h::modulus_squared
      >
    >::wrap("build_normal_equations");
  }


}}}}
