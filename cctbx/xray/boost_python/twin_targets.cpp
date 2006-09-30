#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/twin_targets.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/boost_python/is_polymorphic_workaround.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  cctbx::xray::twin_targets::single_twin_likelihood<double>)

namespace cctbx { namespace xray { namespace twin_targets { namespace boost_python {
  namespace {

    struct least_squares_hemihedral_twinning_on_i_wrappers
    {
      typedef least_squares_hemihedral_twinning_on_i<double> w_t;

      // A thin wrapper that wraps up the derivatives
      static boost::python::tuple
      d_target_d_ab(w_t const& self, scitbx::af::const_ref<std::complex<double> > const& f_model)
      {
        scitbx::af::tiny<scitbx::af::shared<double>, 2> result;
        result = self.d_target_d_ab( f_model );
        return boost::python::make_tuple( result[0], result[1] );
      }

      static void
      wrap()
      {
        using namespace boost::python;


        class_<w_t>("least_squares_hemihedral_twinning_on_i", no_init)
          .def(init<
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 1 indices
               scitbx::af::const_ref< double > const&,                  // 2 i_obs
               scitbx::af::const_ref< double > const&,                  // 3 w_obs
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 4 hl_obs
               sgtbx::space_group const&,                               // 5 space group
               bool const&,                                             // 6 anomalous flag
               double const&,                                           // 7 alpha
               scitbx::mat3<double> const&                              // 8 twin law)
               >
               ((arg_("hkl_obs"),
                 arg_("i_obs"),
                 arg_("w_obs"),
                 arg_("hkl_calc"),
                 arg_("space_group"),
                 arg_("anomalous_flag"),
                 arg_("alpha"),
                 arg_("twin_law")
                 )))
          .def("target", &w_t::target)
          .def("d_target_d_ab", d_target_d_ab )
          .def("d_target_d_fmodel", &w_t::d_target_d_fmodel)
          .def("d_target_d_alpha", &w_t::d_target_d_alpha )
          .def("alpha", (void(w_t::*)(double)) &w_t::alpha)
          .def("alpha", (double(w_t::*)()) &w_t::alpha)
          ;
      }
    };


    struct least_squares_hemihedral_twinning_on_f_wrappers
    {
      typedef least_squares_hemihedral_twinning_on_f<double> w_t;

      // A thin wrapper that wraps up the derivatives
      static boost::python::tuple
      d_target_d_ab(w_t const& self, scitbx::af::const_ref<std::complex<double> > const& f_model)
      {
        scitbx::af::tiny<scitbx::af::shared<double>, 2> result;
        result = self.d_target_d_ab( f_model );
        return boost::python::make_tuple( result[0], result[1] );
      }

      static void
      wrap()
      {
        using namespace boost::python;


        class_<w_t>("least_squares_hemihedral_twinning_on_f", no_init)
          .def(init<
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 1 indices
               scitbx::af::const_ref< double > const&,                  // 2 i_obs
               scitbx::af::const_ref< double > const&,                  // 3 w_obs
               scitbx::af::const_ref< cctbx::miller::index<> > const&,  // 4 hl_obs
               sgtbx::space_group const&,                               // 5 space group
               bool const&,                                             // 6 anomalous flag
               double const&,                                           // 7 alpha
               scitbx::mat3<double> const&                              // 8 twin law)
               >
               ((arg_("hkl_obs"),
                 arg_("f_obs"),
                 arg_("w_obs"),
                 arg_("hkl_calc"),
                 arg_("space_group"),
                 arg_("anomalous_flag"),
                 arg_("alpha"),
                 arg_("twin_law")
                 )))
          .def("target", &w_t::target)
          .def("d_target_d_ab", d_target_d_ab )
          .def("d_target_d_fmodel", &w_t::d_target_d_fmodel)
          .def("d_target_d_alpha", &w_t::d_target_d_alpha )
          .def("alpha", (void(w_t::*)(double)) &w_t::alpha)
          .def("alpha", (double(w_t::*)()) &w_t::alpha)
          ;
      }
    };



    struct hemihedral_r_values_wrappers
    {
      typedef hemihedral_r_values<double> w_t;

      static void
      wrap()
      {
        using namespace boost::python;
        class_<w_t>("hemihedral_r_values", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> > const&,   // 1 obs indices
             scitbx::af::const_ref< cctbx::miller::index<> > const&,   // 2 calc indices
             sgtbx::space_group const&,                                // 3 space group
             bool const&,                                              // 4 anomalous flag
             scitbx::mat3<double> const&                               // 5 twin law
             >
             ((arg_("hkl_obs"),
               arg_("hkl_calc"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("twin_law")
             )))
        .def("r_intensity_abs", &w_t::r_intensity_abs,
             ( arg_("f_obs"),
               arg_("f_model"),
               arg_("selection"),
               arg_("twin_fraction")
             )
            )
        .def("r_intensity_sq", &w_t::r_intensity_sq,
             ( arg_("f_obs"),
               arg_("f_model"),
               arg_("selection"),
               arg_("twin_fraction")
             )
            )
        .def("r_amplitude_abs", &w_t::r_amplitude_abs,
             ( arg_("f_obs"),
               arg_("f_model"),
               arg_("selection"),
               arg_("twin_fraction")
             )
            )
        .def("r_amplitude_sq", &w_t::r_amplitude_sq,
             ( arg_("f_obs"),
               arg_("f_model"),
               arg_("selection"),
               arg_("twin_fraction")
             )
            )




        ;
      }
    };








   struct hemihedral_detwinner_wrappers
   {
     typedef hemihedral_detwinner<double> w_t;

      static boost::python::tuple
      detwin_with_twin_fraction(w_t const& self,
                                scitbx::af::const_ref<double> const& i_obs,
                                scitbx::af::const_ref<double> const& sig_obs,
                                double const& twin_fraction)
      {
        scitbx::af::tiny<scitbx::af::shared<double>, 2> result;
        result = self.detwin_with_twin_fraction(i_obs,sig_obs,twin_fraction);
        return boost::python::make_tuple( result[0], result[1] );
      }


      static boost::python::tuple
      detwin_with_model_data(w_t const& self,
                             scitbx::af::const_ref<double> const& i_obs,
                             scitbx::af::const_ref<double> const& sig_obs,
                             scitbx::af::const_ref<std::complex<double> > const& f_model,
                             double const& twin_fraction)
      {
        scitbx::af::tiny<scitbx::af::shared<double>, 2> result;
        result = self.detwin_with_model_data(i_obs,sig_obs,f_model,twin_fraction);
        return boost::python::make_tuple( result[0], result[1] );
      }


     static void
     wrap()
     {
       using namespace boost::python;
       class_<w_t>("hemihedral_detwinner", no_init)
       .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,   // 1 obs indices
                  scitbx::af::const_ref< cctbx::miller::index<> > const&,   // 2 calc indices
                  sgtbx::space_group const&,                                // 3 space group
                  bool const&,                                              // 4 anomalous flag
                  scitbx::mat3<double> const&                               // 5 twin law
                >
                ((arg_("hkl_obs"),
                  arg_("hkl_calc"),
                  arg_("space_group"),
                  arg_("anomalous_flag"),
                  arg_("twin_law")
                )))
       .def("detwin_with_twin_fraction", detwin_with_twin_fraction,
             ( arg_("i_obs"),
               arg_("sigma_obs"),
               arg_("twin_fraction")
             ) )
       .def("detwin_with_model_data", detwin_with_model_data,
             ( arg_("i_obs"),
               arg_("sigma_obs"),
               arg_("f_model"),
               arg_("twin_fraction")
             ) )
         ;
     }
   };



   struct single_twin_likelihood_wrappers
   {
     typedef single_twin_likelihood<double> w_t;
     static void wrap()
     {
       using namespace boost::python;
       class_<w_t>("single_twin_likelihood", no_init)
         .def(init<
              double const&, double const&,
              double const&, double const&,
              double const&, double const&,
              double const&, double const&,
              bool   const&, bool   const&,
              double const&, double const&,
              double const&, int    const&
              >
              ((arg_("i_obs1"),  arg_("s_obs1"),
                arg_("i_obs2"),  arg_("s_obs2"),
                arg_("f_calc1"), arg_("f_calc2"),
                arg_("eps1"),    arg_("eps2"),
                arg_("centric1"),arg_("centric2"),
                arg_("alpha"),   arg_("beta"),
                arg_("twin_fraction"), arg_("n_quad")
                ))
              )
         .def("log_p" , &w_t::log_p)
         .def("d_log_p_d_f", &w_t::d_log_p_d_f)
         .def("dd_log_p_dd_f", &w_t::dd_log_p_dd_f)
         .def("num_integrate", &w_t::num_integrate)
         .def("laplace_integrate", &w_t::laplace_integrate)
         ;
     }


   };









  }}} // namespace cctbx::xray

  namespace boost_python {

    void wrap_twin_targets()
    {
      twin_targets::boost_python::least_squares_hemihedral_twinning_on_i_wrappers::wrap();
      twin_targets::boost_python::least_squares_hemihedral_twinning_on_f_wrappers::wrap();
      twin_targets::boost_python::hemihedral_r_values_wrappers::wrap();
      twin_targets::boost_python::hemihedral_detwinner_wrappers::wrap();
      twin_targets::boost_python::single_twin_likelihood_wrappers::wrap();
    }

  }

}}
