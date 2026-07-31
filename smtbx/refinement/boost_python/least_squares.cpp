#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/least_squares_matrix_free.h>
#include <smtbx/refinement/weighting_schemes.h>
#include <smtbx/refinement/least_squares_fc_ed_N_beam.h>
#include <cctbx/math/cos_sin_table.h>

void SetRefinementProgressListener(ProgressListener l) {
  smtbx::refinement::least_squares::GetRefinementProgressListener() = l;
}

namespace smtbx { namespace refinement { namespace least_squares {

  ProgressListener& GetRefinementProgressListener() {
    static ProgressListener l = 0;
    return l;
  }

  namespace boost_python {

    using namespace boost::python;
  template <typename FloatType>
  struct wrapper {

    template <class ObjectType, class NormalEquations>
    static void def_init_(class_<ObjectType, bases<builder_base<FloatType> > >& klass) {
      using namespace boost::python;
      typedef void (ObjectType::* build_t)(NormalEquations&,
        IWeightingScheme<FloatType> const&);

      klass.def(
        init<
        NormalEquations&, // normal_equations
        cctbx::xray::observations<FloatType> const&, // reflections
        MaskData<FloatType> const&, // f_mask_data
        IWeightingScheme<FloatType> const&, // weighting_scheme
        boost::optional<FloatType>, // scale_factor
        f_calc_function_base<FloatType> &, // f_calc_function
        scitbx::sparse::matrix<FloatType> const&,
        // jacobian_transpose_matching_grad_fc
        cctbx::xray::fc_correction<FloatType> const&, // exti, swat
        // objective_only=false, may_parallelise_=false, use_openmp=false, max_mem=300
        optional<bool, bool, bool, int>
        >((arg("normal_equations"), arg("reflections"), arg("f_mask_data"),
          arg("weighting_scheme"), arg("scale_factor"),
          arg("f_calc_function"), arg("jacobian_transpose_matching_grad_fc"),
          arg("fc_correction"), arg("objective_only") = false,
          arg("may_parallelise") = false,
          arg("use_openmp") = false,
          arg("max_memory") = 300)))
        .def(
          init<
          cctbx::xray::observations<FloatType> const&, // reflections
          MaskData<FloatType> const&, // f_mask_data
          boost::optional<FloatType>, // scale_factor
          f_calc_function_base<FloatType>&, // f_calc_function
          scitbx::sparse::matrix<FloatType> const&,
          // jacobian_transpose_matching_grad_fc
          cctbx::xray::fc_correction<FloatType> const&, // exti, swat
          // objective_only=false, may_parallelise_=false, use_openmp
          optional<bool, bool, bool>
          >((arg("reflections"), arg("f_mask_data"), arg("scale_factor"),
            arg("f_calc_function"), arg("jacobian_transpose_matching_grad_fc"),
            arg("fc_correction"), arg("objective_only") = false,
            arg("may_parallelise") = false,
            arg("use_openmp") = false,
            arg("max_memory") = 300)))
        .def("build", (build_t)&ObjectType::build)
        ;
    }

    template <class ObjectType>
    static void wrap_init(char const* name,
      boost::python::class_<ObjectType, bases<builder_base<FloatType> > >& klass)
    {
      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
        FloatType,
        scitbx::matrix::sum_of_symmetric_rank_1_updates>
        NormalEquations_BLAS2;

      typedef
        lstbx::normal_equations::non_linear_ls_with_separable_scale_factor<
        FloatType,
        scitbx::matrix::rank_n_update>
        NormalEquations_BLAS3;
      def_init_<ObjectType, NormalEquations_BLAS2>(klass);
      def_init_<ObjectType, NormalEquations_BLAS3>(klass);
      // the accumulators which never form a normal matrix, c.f.
      // smtbx/refinement/least_squares_matrix_free.h
      def_init_<ObjectType, separable_scale_factor_summary<FloatType> >(klass);
      def_init_<ObjectType, separable_scale_factor_product<FloatType> >(klass);
    }

    /// The matrix free accumulators, as targets for the builders above.
    struct matrix_free_accumulators {
      static void wrap() {
        using namespace boost::python;
        {
          typedef separable_scale_factor_summary<FloatType> wt;
          return_value_policy<return_by_value> rbv;
          class_<wt>("separable_scale_factor_summary", no_init)
            .def(init<int, af::const_ref<int> const&, af::const_ref<int> const&>(
              (arg("n_parameters"), arg("block_parameters"), arg("block_sizes"))))
            .def("n_parameters", &wt::n_parameters)
            .def("scale_factor", &wt::scale_factor)
            .def("objective", &wt::objective)
            .def("sum_w_yo_sq", &wt::sum_w_yo_sq_value)
            .def("sum_w_yc_sq", &wt::sum_w_yc_sq_value)
            .def("n_equations", &wt::n_equations_value)
            .def("grad_scale_factor", &wt::grad_scale_factor, rbv)
            .def("right_hand_side", &wt::right_hand_side, rbv)
            .def("blocks", &wt::blocks, rbv)
            ;
        }
        {
          typedef separable_scale_factor_product<FloatType> wt;
          return_value_policy<return_by_value> rbv;
          class_<wt>("separable_scale_factor_product", no_init)
            .def(init<int, FloatType, af::const_ref<FloatType> const&,
                      af::const_ref<FloatType> const&, FloatType>(
              (arg("n_parameters"), arg("scale_factor"),
               arg("grad_scale_factor"), arg("p"), arg("sum_w_yo_sq"))))
            .def("n_parameters", &wt::n_parameters)
            .def("result", &wt::result, rbv)
            ;
        }
      }
    };

    struct normal_equation_building {
      static void wrap_base() {
        typedef builder_base<FloatType> wt;
        return_value_policy<return_by_value> rbv;
        class_<wt, boost::noncopyable>("builder_base", no_init)
          .def("design_matrix", &wt::design_matrix, rbv)
          .def("observables", &wt::observables, rbv)
          //.def("reflections", &wt::reflections)
          .def("f_calc", &wt::f_calc, rbv)
          .def("weights", &wt::weights, rbv)
          .add_static_property("available_threads",
            &wt::get_available_threads,
            &wt::set_available_threads)
          .def("hasOpenMP", &wt::has_openmp)
          .def("has_design_matrix", &wt::has_design_matrix)
          ;
      }

      static void wrap(char const* name) {
        wrap_base();
        typedef build_normal_equations<FloatType> wt;
        class_<wt, bases<builder_base<FloatType> > > klass(name, no_init);
        wrap_init<wt>(name, klass);
      }
    };

    /** StoreType is what the design matrix is held in. Wrapped twice: once at
        FloatType, and once at float, which halves a matrix that is large on
        a protein. The narrow one still accumulates its products in FloatType --
        see build_design_matrix::times -- so it is mixed precision and not
        single precision, and the difference is 2.6e-08 against 3.0e-06.
     */
    template <typename StoreType>
    struct design_matrix_building_ {

      static void wrap(char const* name) {
        using namespace boost::python;
        typedef build_design_matrix<FloatType, StoreType> wt;
        class_<wt, bases<builder_base<FloatType> > > klass(name, no_init);
        wrap_init<wt>(name, klass);
        // the products against the stored matrix, so that it never has to be
        // copied out to be multiplied by
        klass
          .def("times", &wt::times, (arg("x")))
          .def("transpose_times", &wt::transpose_times, (arg("x")))
          .def("is_mixed_precision", &wt::is_mixed_precision)
          .staticmethod("is_mixed_precision")
          ;
      }
    };

    typedef design_matrix_building_<FloatType> design_matrix_building;
    typedef design_matrix_building_<float> design_matrix_building_single;

    struct f_calc_function_wrapper {
      static void wrap_base() {
        using namespace boost::python;
        typedef f_calc_function_base<FloatType> wt;
        typedef void (wt::* evaluate_t1) (miller::index<> const&);
        typedef void (wt::* evaluate_t2) (miller::index<> const&, std::complex<FloatType> const&);
        typedef void (wt::* linearise_t1) (miller::index<> const&);
        typedef void (wt::* linearise_t2) (miller::index<> const&, std::complex<FloatType> const&);
        typedef void (wt::* compute_t1)(miller::index<> const&,
          boost::optional<std::complex<FloatType> > const&, twin_fraction<FloatType> const*, bool);
        typedef void (wt::* compute_t2)(miller::index<> const&, twin_fraction<FloatType> const&,
          std::complex<FloatType> const&, bool);
        class_<wt, boost::noncopyable>("f_calc_function_base", no_init)
          .def("compute", (compute_t1)&wt::compute,
            (arg("index"), arg("f_mask"), arg("fraction"), arg("cumpute_grad") = false))
          .def("compute", (compute_t2)&wt::compute,
            (arg("index"), arg("f_mask"), arg("fraction"), arg("cumpute_grad") = false))
          .def("evaluate", (evaluate_t1)&wt::evaluate,
            (arg("index")))
          .def("evaluate", (evaluate_t2)&wt::evaluate,
            (arg("index"), arg("f_mask")))
          .def("linearise", (linearise_t1)&wt::linearise,
            (arg("index")))
          .def("linearise", (linearise_t2)&wt::linearise,
            (arg("index"), arg("f_mask")))
          .add_property("f_calc", &wt::get_f_calc)
          .add_property("observable", &wt::get_observable)
          ;
      }

      static void wrap_default() {
        using namespace boost::python;
        typedef structure_factors::direct::one_h::std_trigonometry<double,
          structure_factors::direct::one_h::modulus_squared> default_f_calc_func_t;
        typedef f_calc_function_default<FloatType, default_f_calc_func_t> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          boost::shared_ptr<wt> >("f_calc_function_default", no_init)
          .def(init<boost::shared_ptr<default_f_calc_func_t> >(
            (arg("f_calc_function"))))
          ;
      }

      /// The same, over a trigonometry the caller supplies
      /** The default computes exp(i 2pi h.x) with std::cos and std::sin, once
          per reflection per scatterer per symmetry operation, which is the
          bulk of a reflection pass. cctbx carries a lookup table as an
          alternative and smtbx::structure_factors accepts one through
          custom_trigonometry; this is what lets a refinement be handed it.
       */
      static void wrap_default_custom_trigonometry() {
        using namespace boost::python;
        typedef structure_factors::direct::one_h::custom_trigonometry<double,
          structure_factors::direct::one_h::modulus_squared,
          cctbx::math::cos_sin_table> default_f_calc_func_t;
        typedef f_calc_function_default<FloatType, default_f_calc_func_t> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          boost::shared_ptr<wt> >(
            "f_calc_function_default_with_custom_trigonometry", no_init)
          .def(init<boost::shared_ptr<default_f_calc_func_t> >(
            (arg("f_calc_function"))))
          ;
      }

      static void wrap_default_fc() {
        using namespace boost::python;
        typedef structure_factors::direct::one_h::std_trigonometry_fc<double>
          default_f_calc_func_t;
        typedef f_calc_function_default<FloatType, default_f_calc_func_t> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          boost::shared_ptr<wt> >("f_calc_function_default_fc", no_init)
          .def(init<boost::shared_ptr<default_f_calc_func_t> >(
            (arg("f_calc_function"))))
          ;
      }

      static void wrap_caching() {
        using namespace boost::python;
        typedef f_calc_function_with_cache<FloatType> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          boost::shared_ptr<wt> >("f_calc_function_with_cache", no_init)
          .def(init<boost::shared_ptr<f_calc_function_base<FloatType> >, bool>(
            (arg("f_calc_function"), arg("use_cache")=false)))
          ;
      }

      static void wrap_ed_N_beam() {
        using namespace boost::python;
        typedef f_calc_function_ed_N_beam<FloatType> wt;
        typedef f_calc_function_base<FloatType> f_calc_f_t;

        class_<wt, bases<f_calc_function_base<FloatType> >,
          boost::shared_ptr<wt> >("f_calc_function_ed_N_beam", no_init)
          .def(init<N_beam_shared_data<FloatType> const&>(
              (arg("data"),
                arg("params"))))
          ;
      }

      static void wrap_beam_width_cache() {
        using namespace boost::python;
        typedef beam_width_cache<FloatType> wt;

        class_<wt, boost::shared_ptr<wt> >("beam_width_cache", no_init)
          .def(init<beam_width_cache<FloatType> const&>(
            (arg("cache"))))
          .def("find_width", &wt::find_width)
          ;
      }

      static void wrap_mask_data() {
        using namespace boost::python;
        typedef MaskData<FloatType> wt;
        return_value_policy<return_by_value> rbv;
        class_<wt>("MaskData", no_init)
          .def(init<af::const_ref<typename wt::complex_type> const&>((arg("f_mask"))))
          .def(init<cctbx::xray::observations<FloatType> const&,
            sgtbx::space_group const&,
            bool,
            af::const_ref<miller::index<> > const&,
            af::const_ref<typename wt::complex_type> const&>(
              (arg("observations"),
                arg("space_group"), arg("anomalous_flag"),
                arg("indices"), arg("f_mask"))))
          .def("size", &wt::size)
          .def("__len__", &wt::size)
          .def("get", &wt::get, rbv)
          .def("__getitem__", &wt::get, rbv)
          ;
      }

      static void wrap() {
        wrap_base();
        wrap_default();
        wrap_default_custom_trigonometry();
        wrap_default_fc();
        wrap_caching();
        wrap_beam_width_cache();
        wrap_ed_N_beam();
        wrap_mask_data();
      }
    };
  };

  void wrap_least_squares() {
    using namespace boost::python;
    typedef wrapper<double> wrapper_t;

    wrapper_t::matrix_free_accumulators::wrap();
    wrapper_t::normal_equation_building::wrap("build_normal_equations");
    wrapper_t::design_matrix_building::wrap("build_design_matrix");
    wrapper_t::design_matrix_building_single::wrap("build_design_matrix_single");
    wrapper_t::f_calc_function_wrapper::wrap();
  }


}}}}
