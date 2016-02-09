#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/weighting_schemes.h>


namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

  template <typename FloatType,
            class OneMillerIndexLinearisation>
  struct normal_equation_building
  {
    typedef build_normal_equations<FloatType> wt;

    template <class NormalEquations,
              template<typename> class WeightingSchemeType>
    static void def_init(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      // below the comment after each type is the argument name
      // (c.f. smtbx::refinement::build_normal_equations)
      klass.def(
        init<
          NormalEquations &, // normal_equations
          cctbx::xray::observations<FloatType> const &, // reflections
          af::const_ref<std::complex<FloatType> > const &, // f_mask
          WeightingSchemeType<FloatType> const &, // weighting_scheme
          boost::optional<FloatType>, // scale_factor
          OneMillerIndexLinearisation &, // f_calc_function
          scitbx::sparse::matrix<FloatType> const &,
            // jacobian_transpose_matching_grad_fc
          cctbx::xray::extinction_correction<FloatType> const &, // exti
          optional<bool, bool> // objective_only=false, may_parallelise_=false
        >((arg("normal_equations"), arg("reflections"), arg("f_mask"),
           arg("weighting_scheme"), arg("scale_factor"),
           arg("f_calc_function"), arg("jacobian_transpose_matching_grad_fc"),
           arg("extinction"), arg("objective_only")=false,
           arg("may_parallelise")=false)));
    }


    static void wrap(char const *name) {
      using namespace boost::python;

      class_<wt> klass(name, no_init);

      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
          FloatType,
          scitbx::matrix::sum_of_symmetric_rank_1_updates>
        NormalEquations_BLAS2;
      def_init<NormalEquations_BLAS2, mainstream_shelx_weighting>(klass);
      def_init<NormalEquations_BLAS2, unit_weighting            >(klass);
      def_init<NormalEquations_BLAS2, sigma_weighting           >(klass);

#ifdef CCTBX_HAS_LAPACKE
      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
          FloatType,
          scitbx::matrix::rank_n_update>
        NormalEquations_BLAS3;
      def_init<NormalEquations_BLAS3, mainstream_shelx_weighting>(klass);
      def_init<NormalEquations_BLAS3, unit_weighting            >(klass);
      def_init<NormalEquations_BLAS3, sigma_weighting           >(klass);
#endif

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
      structure_factors::direct::one_h::std_trigonometry<
        double,
        structure_factors::direct::one_h::modulus_squared
      >
    >::wrap("build_normal_equations");
  }


}}}}
