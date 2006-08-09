#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <mmtbx/scaling/scaling.h>
#include <mmtbx/scaling/absolute_scaling.h>
#include <mmtbx/scaling/relative_scaling.h>
#include <mmtbx/scaling/twinning.h>
#include <mmtbx/scaling/outlier.h>


namespace mmtbx { namespace scaling {
namespace boost_python{

  // outlier detection  
  void wrap_outlier();
  // twinnning related stuff
  void wrap_twinning();
  // relative scaling
  void wrap_local_scaling_moment_based();
  void wrap_local_scaling_ls_based();
  void wrap_local_scaling_nikonov();
  void wrap_least_squares_on_i();
  void wrap_least_squares_on_f();
  void wrap_least_squares_on_i_wt();
  void wrap_least_squares_on_f_wt();

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
     wilson_single_nll_overloads,
     absolute_scaling::wilson_single_nll, 9, 10)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
     wilson_total_nll_overloads,
     absolute_scaling::wilson_total_nll, 9, 10)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
     ml_normalise_overloads,
     absolute_scaling::ml_normalise, 8, 9)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
     ml_normalise_aniso_overloads,
     absolute_scaling::ml_normalise_aniso, 5, 6)

  void init_module()
  {
    using namespace boost::python;

    // Wilson parameters and other handy stuff
    def("get_d_star_sq_low_limit", get_d_star_sq_low_limit);
    def("get_d_star_sq_high_limit", get_d_star_sq_high_limit);

    def("get_gamma_prot", ( scitbx::af::shared<double>(*)(
        scitbx::af::const_ref<double> const&))
        get_gamma_prot);
    def("gamma_prot", (double(*)(double const&)) gamma_prot);

    def("sigma_prot_sq", (double(*)(double const&, double const&))
        sigma_prot_sq);
    def("get_sigma_prot_sq", (scitbx::af::shared<double>(*)(
                          scitbx::af::const_ref<double> const&,
                          double const&)) get_sigma_prot_sq);



    // Isotropic scaling
    def("wilson_single_nll",
        (double(*)(double const&,
                   double const&,
                   double const&,
                   double const&,
                   double const&,
                   double const&,
                   bool const&,
                   double const&,
                   double const&,
                   bool const&
                   ))
        absolute_scaling::wilson_single_nll,
        wilson_single_nll_overloads((arg_("d_star_sq"),
                                     arg_("f_obs"),
                                     arg_("sigma_f_obs"),
                                     arg_("epsilon"),
                                     arg_("sigma_sq"),
                                     arg_("gamma_prot"),
                                     arg_("centric"),
                                     arg_("p_scale"),
                                     arg_("p_B_wilson"),
                                     arg_("transform")=true))
        );

    def("wilson_total_nll",
        (double(*)(scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<bool> const&,
                   double const&,
                   double const&,
                   bool))
        absolute_scaling::wilson_total_nll,
        wilson_total_nll_overloads((arg_("d_star_sq"),
                                    arg_("f_obs"),
                                    arg_("sigma_f_obs"),
                                    arg_("epsilon"),
                                    arg_("sigma_sq"),
                                    arg_("gamma_prot"),
                                    arg_("centric"),
                                    arg_("p_scale"),
                                    arg_("p_B_wilson"),
                                    arg_("transform")=true))
        );


    def("wilson_single_nll_gradient", (scitbx::af::tiny<double,2>(*)(
                                        double const&,
                                        double const&,
                                        double const&,
                                        double const&,
                                        double const&,
                                        double const&,
                                        bool const&,
                                        double const&,
                                        double const&))
        absolute_scaling::wilson_single_nll_gradient,
          (arg_("d_star_sq"),
           arg_("f_obs"),
           arg_("sigma_f_obs"),
           arg_("epsilon"),
           arg_("sigma_sq"),
           arg_("gamma_prot"),
           arg_("centric"),
           arg_("p_scale"),
           arg_("p_B_wilson"))
        );

    def("wilson_total_nll_gradient", (scitbx::af::tiny<double,2>(*)(
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<double> const&,
                                      scitbx::af::const_ref<bool> const&,
                                      double const&,
                                      double const&))
        absolute_scaling::wilson_total_nll_gradient,
          (arg_("d_star_sq"),
           arg_("f_obs"),
           arg_("sigma_f_obs"),
           arg_("epsilon"),
           arg_("sigma_sq"),
           arg_("gamma_prot"),
           arg_("centric"),
           arg_("p_scale"),
           arg_("p_B_wilson"))
        );

    def ("ml_normalise", (scitbx::af::shared<double>(*)(
                          scitbx::af::const_ref<double> const&,
                          scitbx::af::const_ref<double> const&,
                          scitbx::af::const_ref<double> const&,
                          scitbx::af::const_ref<double> const&,
                          scitbx::af::const_ref<double> const&,
                          scitbx::af::const_ref<bool> const&,
                          double const&,
                          double const&,
                          bool const&))
         absolute_scaling::ml_normalise,
         ml_normalise_overloads(
         (arg_("d_star_sq"),
          arg_("f_obs"),
          arg_("epsilon"),
          arg_("sigma_sq"),
          arg_("gamma_prot"),
          arg_("centric"),
          arg_("p_scale"),
          arg_("p_B_wilson"),
          arg_("wiggle")=true))
         );




    // Anisotropic scaling
    def("wilson_get_aniso_scale",
        (double(*)(cctbx::miller::index<> const&,
                   double const&,
                   double const&,
                   scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_get_aniso_scale);

    def("wilson_single_nll_aniso",
        (double(*)(cctbx::miller::index<> const&,
                   double const&,
                   double const&,
                   double const&,
                   double const&,
                   double const&,
                   bool const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_single_nll_aniso);


    def("wilson_total_nll_aniso",
        (double(*)(scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<double> const&,
                   scitbx::af::const_ref<bool> const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_total_nll_aniso);

    def("wilson_single_nll_aniso_gradient",
        (scitbx::af::shared<double>(*)
         ( cctbx::miller::index<> const&,
           double const&,
           double const&,
           double const&,
           double const&,
           double const&,
           bool const&,
           double const&,
           cctbx::uctbx::unit_cell const&,
           scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_single_nll_aniso_gradient);


    def("wilson_total_nll_aniso_gradient_orbit",
        (scitbx::af::shared<double>(*)
         (scitbx::af::const_ref< cctbx::miller::index<> > const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<bool> const&,
          double const&,
          cctbx::uctbx::unit_cell const&,
          cctbx::sgtbx::space_group const&,
          scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_total_nll_aniso_gradient_orbit);

    def("wilson_total_nll_aniso_gradient",
        (scitbx::af::shared<double>(*)
         (scitbx::af::const_ref< cctbx::miller::index<> > const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<double> const&,
          scitbx::af::const_ref<bool> const&,
          double const&,
          cctbx::uctbx::unit_cell const&,
          scitbx::sym_mat3<double> const&))
        absolute_scaling::wilson_total_nll_aniso_gradient);




    def("ml_normalise_aniso",
        (scitbx::af::shared<double>(*)
         (scitbx::af::const_ref< cctbx::miller::index<> >const&,
          scitbx::af::const_ref<double> const&,
          double const&,
          cctbx::uctbx::unit_cell const& ,
          scitbx::sym_mat3<double> const&,
          bool const&))
         absolute_scaling::ml_normalise_aniso,
         ml_normalise_aniso_overloads(
           (arg_("hkl"),
            arg_("f_obs"),
            arg_("p_scale"),
            arg_("unit_cell"),
            arg_("u_star"),
            arg_("volume_correction_rwgk")=false)));


    def("kernel_normalisation",
        (scitbx::af::shared<double>(*)
         (scitbx::af::const_ref< double > const&,
          scitbx::af::const_ref< double > const&,
          scitbx::af::const_ref< double > const&,
          scitbx::af::const_ref< double > const&,
          double const&))
        absolute_scaling::kernel_normalisation,
        (arg_("d_star_sq_hkl"),
         arg_("I_hkl"),
         arg_("epsilon"),
         arg_("d_star_sq_array"),
         arg_("kernel_width")));


    //outlier detection
    wrap_outlier();
    // twinning related stuff
    wrap_twinning();
    // relative scaling
    wrap_local_scaling_moment_based();
    wrap_local_scaling_ls_based();
    wrap_local_scaling_nikonov();
    wrap_least_squares_on_i();
    wrap_least_squares_on_f();
    wrap_least_squares_on_i_wt();
    wrap_least_squares_on_f_wt();
  }

}}}} // mmtbx::scaling::<anonymous>

BOOST_PYTHON_MODULE(mmtbx_scaling_ext)
{
  mmtbx::scaling::boost_python::init_module();
}
