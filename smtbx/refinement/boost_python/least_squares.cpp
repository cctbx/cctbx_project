#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/weighting_schemes.h>
#include <smtbx/refinement/least_squares_fc_ed.h>
#include <smtbx/refinement/least_squares_fc_ed_n.h>
#include <smtbx/refinement/least_squares_fc_ed_two_beam.h>


namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

    using namespace boost::python;
  template <typename FloatType>
  struct wrapper {

    template <class ObjectType,
      class NormalEquations,
      template<typename> class WeightingSchemeType>
    static void def_init_(class_<ObjectType, bases<builder_base<FloatType> > >& klass) {
      using namespace boost::python;
      typedef void (ObjectType::* build_t)(NormalEquations&,
        WeightingSchemeType<FloatType> const&);

      klass.def(
        init<
        NormalEquations&, // normal_equations
        cctbx::xray::observations<FloatType> const&, // reflections
        MaskData<FloatType> const&, // f_mask_data
        WeightingSchemeType<FloatType> const&, // weighting_scheme
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
      def_init_<ObjectType, NormalEquations_BLAS2, mainstream_shelx_weighting>(klass);
      def_init_<ObjectType, NormalEquations_BLAS2, unit_weighting            >(klass);
      def_init_<ObjectType, NormalEquations_BLAS2, sigma_weighting           >(klass);
      def_init_<ObjectType, NormalEquations_BLAS3, mainstream_shelx_weighting>(klass);
      def_init_<ObjectType, NormalEquations_BLAS3, unit_weighting            >(klass);
      def_init_<ObjectType, NormalEquations_BLAS3, sigma_weighting           >(klass);
    }

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

    struct design_matrix_building {

      static void wrap(char const* name) {
        using namespace boost::python;
        typedef build_design_matrix<FloatType> wt;
        class_<wt, bases<builder_base<FloatType> > > klass(name, no_init);
        wrap_init<wt>(name, klass);
      }
    };

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
          std::auto_ptr<wt> >("f_calc_function_default", no_init)
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
          std::auto_ptr<wt> >("f_calc_function_default_fc", no_init)
          .def(init<boost::shared_ptr<default_f_calc_func_t> >(
            (arg("f_calc_function"))))
          ;
      }

      static void wrap_caching() {
        using namespace boost::python;
        typedef f_calc_function_with_cache<FloatType> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          std::auto_ptr<wt> >("f_calc_function_with_cache", no_init)
          .def(init<boost::shared_ptr<f_calc_function_base<FloatType> >, bool>(
            (arg("f_calc_function"), arg("use_cache")=false)))
          ;
      }

      static void wrap_ed_shared_data() {
        using namespace boost::python;
        typedef ed_shared_data<FloatType> wt;
        return_value_policy<return_by_value> rbv;

        class_<wt, std::auto_ptr<wt> >("ed_shared_data", no_init)
          .def(init<const scitbx::sparse::matrix<FloatType>&,
            f_calc_function_base<FloatType>&,
            sgtbx::space_group const&,
            bool,
            af::shared<FrameInfo<FloatType> >,
            cctbx::xray::thickness<FloatType> const&,
            RefinementParams<FloatType> const&,
            bool, bool>(
              (arg("Jt_matching_grad_fc"),
                arg("f_calc_function"),
                arg("space_group"), arg("anomalous_flag"),
                arg("frames"), arg("thickness"),
                arg("params"), arg("compute_grad"), arg("build") = true)))
          .def("build", &wt::build)
          .add_property("frames", make_getter(&wt::frames, rbv))
          ;
      }

      static void wrap_ed() {
        using namespace boost::python;
        typedef f_calc_function_ed<FloatType> wt;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          std::auto_ptr<wt> >("f_calc_function_ed", no_init)
          .def(init<ed_shared_data<FloatType> const&>(
              (arg("data"))))
          ;
      }

      static void wrap_two_beam_shared_data() {
        using namespace boost::python;
        typedef two_beam_shared_data<FloatType> wt;
        return_value_policy<return_by_value> rbv;

        class_<wt, std::auto_ptr<wt> >("two_beam_shared_data", no_init)
          .def(init<const scitbx::sparse::matrix<FloatType>&,
            f_calc_function_base<FloatType>&,
            sgtbx::space_group const&,
            bool,
            af::shared<FrameInfo<FloatType> >,
            cctbx::xray::thickness<FloatType> const&,
            RefinementParams<FloatType> const&,
            bool, bool>(
              (arg("Jt_matching_grad_fc"),
                arg("f_calc_function"),
                arg("space_group"), arg("anomalous_flag"),
                arg("frames"), arg("thickness"),
                arg("params"), arg("compute_grad"), arg("build") = true)))
          .def("build", &wt::build)
          .add_property("frames", make_getter(&wt::frames, rbv))
          ;
      }

      static void wrap_ed_two_beam() {
        using namespace boost::python;
        typedef f_calc_function_ed_two_beam<FloatType> wt;
        typedef f_calc_function_base<FloatType> f_calc_f_t;

        class_<wt, bases<f_calc_function_base<FloatType> >,
          std::auto_ptr<wt> >("f_calc_function_ed_two_beam", no_init)
          .def(init<two_beam_shared_data<FloatType> const&>(
              (arg("data"),
                arg("params"))))
          ;
      }

      static void wrap_ed_n_shared_data() {
        using namespace boost::python;
        typedef ed_n_shared_data<FloatType> wt;
        typedef f_calc_function_base<FloatType> f_calc_f_t;
        return_value_policy<return_by_value> rbv;

        class_<wt, std::auto_ptr<wt> >("ed_n_shared_data", no_init)
          .def(init<reparametrisation const&,
            f_calc_f_t&,
            cctbx::xray::fc_correction<FloatType> const&,
            sgtbx::space_group const&,
            bool,
            af::shared<FrameInfo<FloatType> >,
            cctbx::xray::thickness<FloatType> const&,
            RefinementParams<FloatType> const&,
            bool, bool>(
              (arg("reparametrisation"),
                arg("f_calc_function"), arg("fc_correction"),
                arg("space_group"), arg("anomalous_flag"),
                arg("frames"), arg("thickness"),
                arg("params"), arg("compute_grad"), arg("build") = true)))
          .def("build", &wt::build)
          .def("process_frame_id", &wt::process_frame_id)
          .add_property("Fcs_k", make_getter(&wt::Fcs_k, rbv))
          ;
      }

      static void wrap_ed_n() {
        using namespace boost::python;
        typedef f_calc_function_ed_n<FloatType> wt;
        typedef f_calc_function_base<FloatType> at;
        class_<wt, bases<f_calc_function_base<FloatType> >,
          std::auto_ptr<wt> >("f_calc_function_ed_n", no_init)
          .def(init<ed_n_shared_data<FloatType> const&>(
              (arg("data"))))
          //.add_property("ratio", &wt::get_ratio)
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
        wrap_default_fc();
        wrap_caching();
        wrap_ed_shared_data();
        wrap_ed();
        wrap_ed_n_shared_data();
        wrap_ed_n();
        wrap_two_beam_shared_data();
        wrap_ed_two_beam();
        wrap_mask_data();
      }
    };
  };
  
  void wrap_least_squares() {
    using namespace boost::python;
    typedef wrapper<double> wrapper_t;

    wrapper_t::normal_equation_building::wrap("build_normal_equations");
    wrapper_t::design_matrix_building::wrap("build_design_matrix");
    wrapper_t::f_calc_function_wrapper::wrap();
  }


}}}}
