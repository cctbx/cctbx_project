#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/max_lik/max_lik.h>

#include <boost/python/list.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace mmtbx { namespace max_lik {
namespace {

  void init_module()
  {
    using namespace boost::python;
    class_<wat_dist>("wat_dist")
      .def("do_wat_dist", &wat_dist::do_wat_dist,
        (arg_("shell"),
         arg_("xyzf"),  arg_("atmrad"), arg_("element_symbol"),
         arg_("uc"),  arg_("sg"), arg_("nxnynz"), arg_("sel_flag"), arg_("rad"), arg_("nshells")))
      .def("data", &wat_dist::data)
      .def("max_number_of_shells", &wat_dist::max_number_of_shells)
      .def("as_xplor_map", &wat_dist::as_xplor_map,
        (arg_("outputfile")))
    ;
    class_<alpha_beta_est>("alpha_beta_est",
                 init<boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      cctbx::sgtbx::space_group const& >((arg_("fo_test"),
                                                          arg_("fm_test"),
                                                          arg_("indices"),
                                                          arg_("space_group"))))
      .def("alpha", &alpha_beta_est::alpha)
      .def("beta",  &alpha_beta_est::beta)
    ;
    class_<fom_and_phase_error>("fom_and_phase_error",
                 init<af::shared<double> const&,
                      af::shared<double> const&,
                      af::shared<double> const&,
                      af::shared<double> const&,
                      cctbx::sgtbx::space_group const&,
                      af::const_ref<cctbx::miller::index<> > const& >((
                                                          arg_("f_obs"),
                                                          arg_("f_model"),
                                                          arg_("alpha"),
                                                          arg_("beta"),
                                                          arg_("space_group"),
                                                          arg_("miller_indices"))))
      .def("fom",         &fom_and_phase_error::fom)
      .def("phase_error", &fom_and_phase_error::phase_error)
    ;
    class_<f_star_w_star_mu_nu>("f_star_w_star_mu_nu",
                 init<af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      cctbx::sgtbx::space_group const&,
                      af::const_ref<cctbx::miller::index<> > const& >((
                                                          arg_("f_obs"),
                                                          arg_("f_model"),
                                                          arg_("alpha"),
                                                          arg_("beta"),
                                                          arg_("space_group"),
                                                          arg_("miller_indices"))))
      .def("f_star", &f_star_w_star_mu_nu::f_star)
      .def("w_star", &f_star_w_star_mu_nu::w_star)
      .def("mu",     &f_star_w_star_mu_nu::mu)
      .def("nu",     &f_star_w_star_mu_nu::nu)
      .def("k1",     &f_star_w_star_mu_nu::f_star)
      .def("k2",     &f_star_w_star_mu_nu::w_star)
      .def("ls_star_1", &f_star_w_star_mu_nu::mu)
      .def("ls_star_2", &f_star_w_star_mu_nu::nu)
      .def("number_of_f_star_zero", &f_star_w_star_mu_nu::number_of_f_star_zero)
    ;
    class_<sasha_error_calculator>("sasha_error_calculator",
                 init<af::const_ref<vec3<double> > const&,
                      af::const_ref<vec3<double> > const&,
                      af::const_ref<vec3<double> > const&,
                      af::const_ref<vec3<double> > const&,
                      af::const_ref<std::string> const&,
                      af::const_ref<std::string> const&,
                      cctbx::uctbx::unit_cell const&,
                      cctbx::sgtbx::space_group const&,
                      double >((arg_("r1f"),
                                    arg_("r1c"),
                                    arg_("r2f"),
                                    arg_("r2c"),
                                    arg_("lab1"),
                                    arg_("lab2"),
                                    arg_("uc"),
                                    arg_("sg"),
                                    arg_("rad"))))
      .def("doptimal", &sasha_error_calculator::doptm)
      .def("dregular", &sasha_error_calculator::distm)
    ;
    class_<peak_clustering>("peak_clustering",
                 init<af::const_ref<vec3<double> > const&,
                      af::const_ref<vec3<double> > const&,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      cctbx::uctbx::unit_cell const&,
                      double const& >((arg_("r1f"),
                                       arg_("r2f"),
                                       arg_("h1"),
                                       arg_("h2"),
                                       arg_("uc"),
                                       arg_("cutoff"))))
      .def("sites", &peak_clustering::sites)
      .def("heights", &peak_clustering::heights)
    ;
  }

} // namespace <anonymous>
}} // namespace mmtbx::max_lik

BOOST_PYTHON_MODULE(mmtbx_max_lik_ext)
{
  mmtbx::max_lik::init_module();
}
