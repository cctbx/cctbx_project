#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <mmtbx/scaling/twinning.h>

namespace mmtbx { namespace scaling {



namespace{


  struct twin_r_wrapper
  {
    typedef twinning::twin_r<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("twin_r", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const& >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("symop") )))
        .def("r_abs_value", &w_t::r_abs_value)
        .def("r_sq_value", &w_t::r_sq_value)
        .def("r_abs_pair", &w_t::r_abs_pair)
        .def("r_sq_pair", &w_t::r_abs_pair)
        ;
    }

  };




  struct l_test_wrapper
  {
    typedef twinning::l_test<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("l_test", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             long const&,
             long const&,
             long const&,
             std::size_t const& >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("parity_h"),
               arg_("parity_k"),
               arg_("parity_l"),
               arg_("max_delta_h") )))
        .def("mean_l", &w_t::mean_l)
        .def("mean_l2", &w_t::mean_l2)
        .def("ml_alpha", &w_t::ml_alpha)
        .def("cumul", &w_t::cumul)
        ;
    }

  };


  struct detwin_wrapper
  {
    typedef twinning::detwin<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("detwin", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&
             >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("sigma"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("twin_law") )))
        .def("detwin_with_alpha", &w_t::detwin_with_alpha)
        .def("detwinned_i", &w_t::detwinned_i)
        .def("detwinned_sigi", &w_t::detwinned_sigi)
        .def("detwinned_hkl", &w_t::detwinned_hkl)
        .def("completeness", &w_t::completeness)
        .def("location", &w_t::location)
        ;
    }

  };


  struct h_test_wrapper
  {
    typedef twinning::h_test<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("h_test", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&,
             double const&
             >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("sigma"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("twin_law"),
               arg_("fraction")) ))
        .def("h_array", &w_t::h_array)
        .def("h_values", &w_t::h_values)
        .def("h_cumul_obs", &w_t::h_cumul_obs)
        .def("h_cumul_fit", &w_t::h_cumul_fit)
        .def("distance", &w_t::distance)
        .def("alpha", &w_t::alpha)
        .def("mean_h", &w_t::mean_h)
        .def("mean_h2", &w_t::mean_h2)
        ;
    }

  };


  struct ml_murray_rust_wrapper
  {
    typedef twinning::ml_murray_rust<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("ml_murray_rust", no_init)
        .def(init<
             scitbx::af::const_ref<double> const&,
             scitbx::af::const_ref<double> const&,
             scitbx::af::const_ref< cctbx::miller::index<> > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&,
             long const&
             > (( arg_("z"),
                  arg_("sig_z"),
                  arg_("indices"),
                  arg_("space_group"),
                  arg_("anomalous_flag"),
                  arg_("twin_law"),
                  arg_("n_hermite"))
                ))
        .def("num_int", &w_t::num_int,
             (( arg_("ito1"),
                arg_("sito1"),
                arg_("ito2"),
                arg_("sito2"),
                arg_("low_sig"),
                arg_("high_sig"),
                arg_("twin_fraction"),
                arg_("n")
                ))
            )
        .def("p_raw", &w_t::p_raw)
        .def("log_p_given_t", &w_t::log_p_given_t)
        .def("fast_log_p_given_t", &w_t::fast_log_p_given_t)

        ;

    }


  };


  struct ml_twin_with_ncs
  {
    typedef twinning::ml_twin_with_ncs<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("ml_twin_with_ncs", no_init)
        .def(init<
             scitbx::af::const_ref<double> const&,
             scitbx::af::const_ref<double> const&,
             scitbx::af::const_ref< cctbx::miller::index<> > const&,
             scitbx::af::const_ref< std::size_t > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&,
             cctbx::uctbx::unit_cell const&,
             long const&
             > (( arg_("z"),
                  arg_("sig_z"),
                  arg_("indices"),
                  arg_("bins"),
                  arg_("space_group"),
                  arg_("anomalous_flag"),
                  arg_("twin_law"),
                  arg_("unit_cell"),
                  arg_("n_quadrature"))
                ))
        .def("p_raw", &w_t::p_raw)
        .def("num_int", &w_t::num_int )
        .def("p_tot_given_t_and_coeff", &w_t::p_tot_given_t_and_coeff)
        ;

    }


  };



  struct very_quick_erf_wrapper
  {
    typedef twinning::very_quick_erf<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("very_quick_erf", no_init)
        .def(init<double const& > ((arg_("step_size"))))
        .def("erf", &w_t::erf)
        .def("loop_for_timings", &w_t::loop_for_timings, (
          arg_("number_of_iterations"), arg_("optimized")))
        ;
    }
  };



  struct quick_ei0_wrapper
  {
    typedef twinning::quick_ei0<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("quick_ei0", no_init)
        .def(init<int const& > ((arg_("n_points"))))
        .def("ei0", &w_t::ei0)
        .def("loop_for_timings", &w_t::loop_for_timings, (
          arg_("number_of_iterations"), arg_("optimized")))
        ;
    }
  };

}  // namespace <anonymous>

namespace boost_python {

  void wrap_twinning()
  {
    twin_r_wrapper::wrap();
    l_test_wrapper::wrap();
    detwin_wrapper::wrap();
    h_test_wrapper::wrap();
    ml_murray_rust_wrapper::wrap();
    ml_twin_with_ncs::wrap();
    very_quick_erf_wrapper::wrap();
    quick_ei0_wrapper::wrap();
  }

}}}
