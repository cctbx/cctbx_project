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
  struct wrapper {
    template <class ObjectType,
      class NormalEquations,
      template<typename> class WeightingSchemeType>
    static void def_init_(boost::python::class_<ObjectType> &klass) {
      using namespace boost::python;
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
        // objective_only=false, may_parallelise=false, use_openpm=false
        optional<bool, bool, bool>
        >((arg("normal_equations"), arg("reflections"), arg("f_mask"),
          arg("weighting_scheme"), arg("scale_factor"),
          arg("f_calc_function"), arg("jacobian_transpose_matching_grad_fc"),
          arg("extinction"), arg("objective_only") = false,
          arg("may_parallelise") = false, arg("use_openmp") = false)));
    }

    template <class ObjectType>
    static void wrap_init(char const *name,
      boost::python::class_<ObjectType> &klass)
    {
      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
        FloatType,
        scitbx::matrix::sum_of_symmetric_rank_1_updates>
        NormalEquations_BLAS2;
      def_init_<ObjectType, NormalEquations_BLAS2, mainstream_shelx_weighting>(klass);
      def_init_<ObjectType, NormalEquations_BLAS2, unit_weighting            >(klass);
      def_init_<ObjectType, NormalEquations_BLAS2, sigma_weighting           >(klass);

      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
        FloatType,
        scitbx::matrix::rank_n_update>
        NormalEquations_BLAS3;
      def_init_<ObjectType, NormalEquations_BLAS3, mainstream_shelx_weighting>(klass);
      def_init_<ObjectType, NormalEquations_BLAS3, unit_weighting            >(klass);
      def_init_<ObjectType, NormalEquations_BLAS3, sigma_weighting           >(klass);
    }

    struct normal_equation_building {
      typedef build_normal_equations<FloatType> wt;

      static void wrap(char const *name) {
        using namespace boost::python;
        class_<wt> klass(name, no_init);
        wrap_init<wt>(name, klass);
        klass
          .def("observables", &wt::observables)
          .def("f_calc", &wt::f_calc)
          .def("weights", &wt::weights)
          .add_static_property("available_threads",
            &wt::get_available_threads,
            &wt::set_available_threads)
          .def("hasOpenMP", &wt::has_openmp)
            ;
      }
    };

    struct design_matrix_building {
      typedef build_design_matrix<FloatType> wt;

      static void wrap(char const *name) {
        using namespace boost::python;
        class_<wt> klass(name, no_init);
        wrap_init<wt>(name, klass);
        klass
          .def("observables", &wt::observables)
          .def("f_calc", &wt::f_calc)
          .def("weights", &wt::weights)
          .add_static_property("available_threads",
            &wt::get_available_threads,
            &wt::set_available_threads)
          .def("hasOpenMP", &wt::has_openmp)
          .def("design_matrix", &wt::design_matrix)
          ;
      }
    };
  };

  void wrap_least_squares() {
    using namespace boost::python;
    typedef wrapper<double,
      structure_factors::direct::one_h::std_trigonometry<
        double,
        structure_factors::direct::one_h::modulus_squared>
    > wrapper_t;
    wrapper_t::normal_equation_building::wrap("build_normal_equations");
    wrapper_t::design_matrix_building::wrap("build_design_matrix");
  }


}}}}
